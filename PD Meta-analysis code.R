# ============================================================================
# META-ANALYSIS FOR PARKINSON'S DISEASE (Blood & Brain)
# Approach A: Within-dataset processing + RRA + weighted Z
# Includes scatter plots of effect sizes and all standard diagnostics
# ============================================================================
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("ggplot2", "reshape2", "corrplot", "dplyr", "limma")
for(p in pkgs) {
  if (!require(p, character.only = TRUE, quietly = TRUE))
    install.packages(p)
  library(p, character.only = TRUE)
}

# ----------------------------- 1. SETUP ------------------------------------
pd_dir <- "C:/Users/User/Desktop/Parkinson's Disease-Feb2026/"

# Final dataset lists (blood and good brain)
pd_blood <- c("GSE75249","GSE6613","GSE49126","GSE161199")
pd_brain <- c("GSE42966","GSE19587","GSE7621","GSE8397A","GSE8397B")
pd_all <- c(pd_blood, pd_brain)

# ----------------------------- 2. PREPROCESSING ----------------------------
preprocess_dataset <- function(gse, base_dir) {
  file_raw <- file.path(base_dir, paste0(gse, "_exp_norm_clean.csv"))
  if(!file.exists(file_raw)) {
    cat("  Missing raw file:", basename(file_raw), "\n")
    return(FALSE)
  }
  dat <- read.csv(file_raw, check.names = FALSE)
  if(colnames(dat)[1] != "GENE") colnames(dat)[1] <- "GENE"
  # Remove rows with invalid gene symbols (empty, NA, or only digits)
  gene <- as.character(dat$GENE)
  invalid <- is.na(gene) | gene == "" | grepl("^[0-9]+$", gene)
  if(any(invalid)) {
    dat <- dat[!invalid, ]
    cat("  ", gse, ": removed", sum(invalid), "invalid gene rows\n")
  }
  expr <- as.matrix(dat[, -1])
  mode(expr) <- "numeric"
  
  vals <- as.vector(expr)
  vals <- vals[!is.na(vals)]
  if(length(vals) > 0 && (median(vals) > 10 || mean(vals > 10, na.rm = TRUE) > 0.9)) {
    expr <- log2(expr + 1)
    cat("  ", gse, ": applied log2(x+1)\n")
  } else {
    cat("  ", gse, ": already log‑transformed (or raw but low median)\n")
  }
  
  out_df <- data.frame(GENE = dat[,1], expr, check.names = FALSE)
  file_pre <- file.path(base_dir, paste0(gse, "_exp_norm_preprocessed.csv"))
  write.csv(out_df, file_pre, row.names = FALSE)
  return(TRUE)
}

# ----------------------------- 3. COLLAPSE DUPLICATES ----------------------
collapse_dataset <- function(gse, base_dir) {
  file_pre <- file.path(base_dir, paste0(gse, "_exp_norm_preprocessed.csv"))
  if(!file.exists(file_pre)) return(FALSE)
  dat <- read.csv(file_pre, check.names = FALSE)
  if(colnames(dat)[1] != "GENE") colnames(dat)[1] <- "GENE"
  gene <- dat[,1]
  expr <- as.matrix(dat[, -1])
  mode(expr) <- "numeric"
  
  if(anyDuplicated(gene)) {
    expr_coll <- aggregate(expr, by = list(GENE = gene), FUN = median, na.rm = TRUE)
    rownames(expr_coll) <- expr_coll$GENE
    expr_coll <- expr_coll[, -1, drop = FALSE]
    expr_gene <- as.matrix(expr_coll)
  } else {
    rownames(expr) <- gene
    expr_gene <- expr
  }
  
  # Remove rows where all expression values are exactly 0
  all_zero <- apply(expr_gene, 1, function(x) all(x == 0, na.rm = TRUE))
  if(sum(all_zero) > 0) {
    expr_gene <- expr_gene[!all_zero, , drop = FALSE]
    cat("  ", gse, ": removed", sum(all_zero), "rows with all zero expression\n")
  }
  
  out_df <- data.frame(GENE = rownames(expr_gene), expr_gene, check.names = FALSE)
  file_coll <- file.path(base_dir, paste0(gse, "_genelevel_norm_collapse.csv"))
  write.csv(out_df, file_coll, row.names = FALSE)
  cat("  ", gse, ": collapsed to", nrow(expr_gene), "genes\n")
  return(TRUE)
}

# ----------------------------- 4. SAFETY LOG2 (POST‑COLLAPSE) --------------
safety_log2 <- function(gse_list, base_dir) {
  cat("\n========== SAFETY LOG2 (max > 50) ==========\n")
  for(gse in gse_list) {
    f <- file.path(base_dir, paste0(gse, "_genelevel_norm_collapse.csv"))
    if(!file.exists(f)) next
    dat <- read.csv(f, check.names = FALSE)
    if(colnames(dat)[1] != "GENE") colnames(dat)[1] <- "GENE"
    rownames(dat) <- dat$GENE
    expr <- as.matrix(dat[, -1])
    mode(expr) <- "numeric"
    max_val <- max(expr, na.rm = TRUE)
    if(max_val > 50) {
      expr_log <- log2(expr + 1)
      new_dat <- data.frame(GENE = rownames(expr_log), expr_log, check.names = FALSE)
      write.csv(new_dat, f, row.names = FALSE)
      cat("  SAFETY: applied log2 to", gse, "(max was", round(max_val,1), ")\n")
    } else {
      cat("  SAFETY: skipped", gse, "(max", round(max_val,1), ")\n")
    }
  }
}

# ----------------------------- 5. FULL COMPATIBILITY -----------------------
full_compatibility <- function(gse_list, base_dir, tissue_name) {
  cat("\n========== FULL COMPATIBILITY: PD", tissue_name, "==========\n")
  expr_list <- list()
  for(gse in gse_list) {
    f <- file.path(base_dir, paste0(gse, "_genelevel_norm_collapse.csv"))
    if(!file.exists(f)) next
    dat <- read.csv(f, check.names = FALSE)
    rownames(dat) <- dat[,1]
    mat <- as.matrix(dat[,-1])
    # Remove rows with any NA/Inf
    finite <- apply(mat, 1, function(x) all(is.finite(x)))
    if(sum(finite) == 0) next
    mat <- mat[finite, , drop=FALSE]
    # Also remove rows with all zero (if not already removed)
    all_zero <- apply(mat, 1, function(x) all(x == 0, na.rm=TRUE))
    if(sum(all_zero) > 0) mat <- mat[!all_zero, , drop=FALSE]
    expr_list[[gse]] <- mat
    cat("  ", gse, "genes:", nrow(mat), "samples:", ncol(mat),
        "median:", round(median(mat),2), "max:", round(max(mat),2), "\n")
  }
  if(length(expr_list) < 2) {
    cat("  Not enough datasets.\n")
    return()
  }
  common <- Reduce(intersect, lapply(expr_list, rownames))
  cat("  Common genes:", length(common), "\n")
  if(length(common) < 100) {
    cat("  Too few common genes – skipping full plots.\n")
    return()
  }
  expr_common <- lapply(expr_list, function(x) x[common, , drop=FALSE])
  
  # Summary statistics
  stats_df <- data.frame(
    Dataset = names(expr_common),
    nSamples = sapply(expr_common, ncol),
    Median = sapply(expr_common, function(x) median(x, na.rm=TRUE)),
    Mean = sapply(expr_common, function(x) mean(x, na.rm=TRUE)),
    SD = sapply(expr_common, function(x) sd(x, na.rm=TRUE)),
    IQR = sapply(expr_common, function(x) IQR(x, na.rm=TRUE))
  )
  write.csv(stats_df, file.path(base_dir, paste0("PD_", tissue_name, "_stats.csv")), row.names=FALSE)
  
  # Prepare long format for plotting
  long_list <- list()
  for(i in seq_along(expr_common)) {
    mat <- expr_common[[i]]
    df <- reshape2::melt(mat, varnames=c("Gene","Sample"), value.name="Expression")
    df <- df[is.finite(df$Expression), ]
    df$Dataset <- names(expr_common)[i]
    long_list[[i]] <- df
  }
  long_all <- do.call(rbind, long_list)
  if(nrow(long_all) == 0) {
    cat("  No finite values after melting – cannot plot.\n")
    return()
  }
  
  # Boxplot (trim outliers to 1st–99th percentile)
  p_box <- ggplot(long_all, aes(x=Dataset, y=Expression, fill=Dataset)) +
    geom_boxplot(outlier.shape=NA) + 
    coord_cartesian(ylim = quantile(long_all$Expression, c(0.01, 0.99), na.rm=TRUE)) +
    theme_minimal() +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    labs(title=paste("PD", tissue_name, "- Expression distributions"), y="Expression")
  ggsave(file.path(base_dir, paste0("PD_", tissue_name, "_boxplot.png")), p_box, width=14, height=6)
  
  # Density
  p_dens <- ggplot(long_all, aes(x=Expression, color=Dataset)) +
    geom_density() + theme_minimal() +
    labs(title=paste("PD", tissue_name, "- Density curves"))
  ggsave(file.path(base_dir, paste0("PD_", tissue_name, "_density.png")), p_dens, width=10, height=6)
  
  # ---------------------------- PCA (IMPROVED VISIBILITY) ---------------------------
  combined <- do.call(cbind, expr_common)
  pca_input <- t(combined)
  pca_input <- pca_input[complete.cases(pca_input), ]
  if(ncol(pca_input) > 1 && nrow(pca_input) > 1) {
    const <- apply(pca_input, 2, function(x) length(unique(x)) == 1)
    pca_input <- pca_input[, !const, drop=FALSE]
    if(ncol(pca_input) >= 2) {
      pca_res <- prcomp(pca_input, center=TRUE, scale=TRUE)
      pca_df <- data.frame(pca_res$x[,1:2],
                           Dataset = rep(names(expr_common), times=sapply(expr_common, ncol)))
      pca_df <- pca_df[rownames(pca_input), ]
      
      # --- Diagnostic print (optional) ---
      cat("PCA coordinates for", tissue_name, ":\n")
      print(table(pca_df$Dataset))
      if(tissue_name == "Blood" && "GSE75249" %in% pca_df$Dataset) {
        sub <- pca_df[pca_df$Dataset == "GSE75249", ]
        cat("GSE75249 PC1 range:", range(sub$PC1), " PC2 range:", range(sub$PC2), "\n")
      }
      
      # --- Dynamic axis limits (add 5% margin) ---
      x_range <- range(pca_df$PC1)
      y_range <- range(pca_df$PC2)
      x_margin <- diff(x_range) * 0.05
      y_margin <- diff(y_range) * 0.05
      
      if (tissue_name == "Blood") {
        # Make GSE75249 stand out: orange circle with black border, larger size
        # Others: semi‑transparent smaller points
        blood_colors <- c("GSE161199" = "red",
                          "GSE49126"  = "green",
                          "GSE6613"   = "blue",
                          "GSE75249"  = "orange")  # used for legend only
        
        p_pca <- ggplot() +
          # Other datasets (smaller, semi‑transparent)
          geom_point(data = subset(pca_df, Dataset != "GSE75249"),
                     aes(x=PC1, y=PC2, color=Dataset),
                     alpha = 0.6, size = 2) +
          # GSE75249 (large, solid orange with black border)
          geom_point(data = subset(pca_df, Dataset == "GSE75249"),
                     aes(x=PC1, y=PC2),
                     color = "black", fill = "orange", shape = 21, size = 4, stroke = 0.8) +
          scale_color_manual(values = blood_colors) +
          coord_cartesian(xlim = c(x_range[1] - x_margin, x_range[2] + x_margin),
                          ylim = c(y_range[1] - y_margin, y_range[2] + y_margin)) +
          theme_minimal() +
          labs(title = paste("PD", tissue_name, "- PCA (GSE75249 as large orange circles)"))
      } else {
        # Brain: custom distinct colours, all points larger and opaque
        brain_colors <- c("GSE42966" = "#E41A1C",   # red
                          "GSE19587" = "#377EB8",   # blue
                          "GSE7621"  = "#4DAF4A",   # green
                          "GSE8397A" = "#984EA3",   # purple
                          "GSE8397B" = "#FF7F00")   # orange
        p_pca <- ggplot(pca_df, aes(x=PC1, y=PC2, color=Dataset)) +
          geom_point(alpha = 0.9, size = 3) +
          scale_color_manual(values = brain_colors) +
          coord_cartesian(xlim = c(x_range[1] - x_margin, x_range[2] + x_margin),
                          ylim = c(y_range[1] - y_margin, y_range[2] + y_margin)) +
          theme_minimal() +
          labs(title = paste("PD", tissue_name, "- PCA"))
      }
      ggsave(file.path(base_dir, paste0("PD_", tissue_name, "_PCA.png")), p_pca, width=10, height=8)
      cat("  PCA plot saved with enhanced visibility.\n")
    }
  }
  
  # Correlation heatmap of mean expression vectors
  mean_expr <- sapply(expr_common, function(x) rowMeans(x, na.rm=TRUE))
  if(ncol(mean_expr) > 1) {
    mean_expr <- na.omit(mean_expr)
    if(nrow(mean_expr) > 1) {
      cormat <- cor(mean_expr, use="pairwise.complete.obs")
      png(file.path(base_dir, paste0("PD_", tissue_name, "_cor_heatmap.png")), width=800, height=800)
      corrplot::corrplot(cormat, method="color", type="upper", order="hclust",
                         tl.col="black", tl.srt=45, title=paste("PD", tissue_name, "- Mean expression correlation"))
      dev.off()
    }
  }
  cat("  Full compatibility completed.\n")
}

# ----------------------------- 6. RUN PIPELINE -----------------------------
cat("\n========== PREPROCESSING ==========\n")
for(gse in pd_all) preprocess_dataset(gse, pd_dir)

cat("\n========== COLLAPSING ==========\n")
for(gse in pd_all) collapse_dataset(gse, pd_dir)

cat("\n========== SAFETY LOG2 ==========\n")
safety_log2(pd_all, pd_dir)

cat("\n========== FULL COMPATIBILITY ==========\n")
full_compatibility(pd_blood, pd_dir, "Blood")
full_compatibility(pd_brain, pd_dir, "Brain")

cat("\n✅ PD pipeline completed. Check the PD directory for outputs.\n")
# ----------------------------- 0. INSTALL & LOAD PACKAGES ------------------
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("limma", "RobustRankAggreg", "pheatmap", "ggplot2", "ggrepel",
          "reshape2", "RColorBrewer", "WGCNA", "UpSetR", "corrplot",
          "ggVennDiagram", "flashClust", "metafor")
for(p in pkgs) {
  if (!require(p, character.only = TRUE, quietly = TRUE))
    BiocManager::install(p, ask = FALSE)
  library(p, character.only = TRUE)
}
enableWGCNAThreads()

# ----------------------------- 1. PATHS & DATASET LISTS --------------------
data_dir <- "C:/Users/User/Desktop/Parkinson's Disease-Feb2026/"
out_dir <- file.path(data_dir, "PD_MetaAnalysis_Results")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

blood_datasets <- c("GSE75249","GSE6613","GSE49126","GSE161199")
brain_datasets <- c("GSE42966","GSE19587","GSE7621","GSE8397A","GSE8397B")

# ----------------------------- 2. PER-DATASET limma (with optional batch) ---
run_limma_with_batch <- function(gse, base_dir) {
  expr_file <- file.path(base_dir, paste0(gse, "_genelevel_norm_collapse.csv"))
  pheno_file <- file.path(base_dir, paste0(gse, "_pheno.csv"))
  if (!file.exists(expr_file) || !file.exists(pheno_file))
    stop(paste("Missing file for", gse))
  
  expr <- read.csv(expr_file, check.names = FALSE)
  pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)
  
  # Align samples
  sample_cols <- intersect(colnames(expr)[-1], pheno$GSM)
  if (length(sample_cols) == 0) stop("No matching GSM IDs")
  expr <- expr[, c("GENE", sample_cols)]
  pheno <- pheno[match(sample_cols, pheno$GSM), ]
  
  # Build expression matrix
  gene_symbols <- as.character(expr$GENE)
  mat <- as.matrix(expr[, -1])
  rownames(mat) <- make.unique(gene_symbols)
  mode(mat) <- "numeric"
  
  # Filter low-variance genes (20% quantile)
  gene_vars <- apply(mat, 1, var, na.rm = TRUE)
  mat <- mat[gene_vars > quantile(gene_vars, 0.20, na.rm = TRUE), , drop = FALSE]
  
  # Optional: batch effect correction within dataset (if Batch column exists)
  if ("Batch" %in% colnames(pheno)) {
    batch <- factor(pheno$Batch)
    design <- model.matrix(~0 + Diagnosis + batch, data = pheno)
  } else {
    design <- model.matrix(~0 + Diagnosis, data = pheno)
  }
  colnames(design) <- gsub("Diagnosis", "", colnames(design))
  
  # Contrast: PD vs Control
  contrast <- makeContrasts(PD - Control, levels = design)
  fit <- lmFit(mat, design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  res <- topTable(fit2, number = Inf, sort.by = "none")
  
  # Correct Z-score (two-tailed, sign = direction of effect)
  p_two <- 2 * pt(-abs(res$t), df = fit2$df.total)
  res$Z <- sign(res$t) * qnorm(1 - p_two/2)
  res$Gene <- rownames(res)
  res <- res[, c("Gene", "logFC", "AveExpr", "t", "P.Value", "adj.P.Val", "Z")]
  return(res)
}

# ----------------------------- 3. META-ANALYSIS FUNCTION -------------------
run_meta_analysis <- function(dataset_list, base_dir, out_dir, tissue_name) {
  cat("\n========== Processing", tissue_name, "==========\n")
  tissue_dir <- file.path(out_dir, tissue_name)
  dir.create(tissue_dir, showWarnings = FALSE)
  
  z_list <- list()
  fc_list <- list()
  sample_sizes <- c()
  
  for (gse in dataset_list) {
    cat("  Running limma on", gse, "... ")
    deg <- tryCatch(run_limma_with_batch(gse, base_dir), error = function(e) NULL)
    if (is.null(deg)) {
      cat("FAILED\n")
      next
    }
    write.csv(deg, file = file.path(tissue_dir, paste0(gse, "_DEG_results.csv")), row.names = FALSE)
    z <- deg$Z; names(z) <- deg$Gene
    fc <- deg$logFC; names(fc) <- deg$Gene
    # Remove NA/empty/duplicated
    valid <- !is.na(names(z)) & names(z) != "" & !duplicated(names(z))
    z_list[[gse]] <- z[valid]
    fc_list[[gse]] <- fc[valid]
    
    pheno <- read.csv(file.path(base_dir, paste0(gse, "_pheno.csv")), stringsAsFactors = FALSE)
    sample_sizes[gse] <- nrow(pheno)
    cat(nrow(deg), "genes\n")
  }
  
  if (length(z_list) < 2) stop("Not enough datasets processed.")
  
  # Build Z-matrix
  all_genes <- sort(unique(unlist(lapply(z_list, names))))
  zmat <- sapply(z_list, function(z) z[match(all_genes, names(z))])
  rownames(zmat) <- all_genes
  write.csv(zmat, file = file.path(tissue_dir, "Z_matrix_raw.csv"))
  
  # RRA for up- and down-regulated genes
  # Up: sort by decreasing Z (higher Z = more up)
  # Down: sort by increasing Z (lower Z = more down)
  genes_keep <- all_genes  # we keep all for RRA; later filter by consistency
  rra_up <- aggregateRanks(
    glist = lapply(z_list, function(z) names(sort(z[names(z) %in% genes_keep], decreasing = TRUE))),
    N = length(genes_keep)
  )
  rra_down <- aggregateRanks(
    glist = lapply(z_list, function(z) names(sort(z[names(z) %in% genes_keep], decreasing = FALSE))),
    N = length(genes_keep)
  )
  
  # Consistency filter: gene must have same sign in >= 60% of datasets
  min_datasets <- ceiling(length(dataset_list) * 0.6)
  dir_up <- rowSums(sign(zmat) == 1, na.rm = TRUE) >= min_datasets
  dir_down <- rowSums(sign(zmat) == -1, na.rm = TRUE) >= min_datasets
  consistent_genes <- rownames(zmat)[dir_up | dir_down]
  
  # Meta-DEGs: RRA p-value < 0.05 AND consistent
  meta_candidates <- intersect(
    unique(c(rra_up$Name[rra_up$Score < 0.05], rra_down$Name[rra_down$Score < 0.05])),
    consistent_genes
  )
  if (length(meta_candidates) == 0) {
    warning("No meta-DEGs found; relaxing RRA threshold to 0.1")
    meta_candidates <- intersect(
      unique(c(rra_up$Name[rra_up$Score < 0.1], rra_down$Name[rra_down$Score < 0.1])),
      consistent_genes
    )
  }
  
  # Weighted meta-Z (using scaled Z-matrix and sqrt(n) weights)
  zmat_scaled <- scale(zmat)
  write.csv(zmat_scaled, file = file.path(tissue_dir, "Z_matrix_scaled.csv"))
  weights <- sqrt(sample_sizes[colnames(zmat_scaled)])
  z_meta <- zmat_scaled[meta_candidates, , drop = FALSE]
  w_mat <- sweep(matrix(1, nrow(z_meta), ncol(z_meta)), 2, weights, "*")
  w_mat[is.na(z_meta)] <- NA
  weighted_Z <- rowSums(z_meta * w_mat, na.rm = TRUE) / rowSums(w_mat, na.rm = TRUE)
  
  # Build final data frame
  meta_df <- data.frame(
    Gene = meta_candidates,
    Weighted_Z = weighted_Z,
    Direction = ifelse(weighted_Z > 0, "Upregulated", "Downregulated"),
    RRA_score = pmin(
      rra_up$Score[match(meta_candidates, rra_up$Name)],
      rra_down$Score[match(meta_candidates, rra_down$Name)],
      na.rm = TRUE
    ),
    Up_datasets = rowSums(sign(zmat[meta_candidates, ]) == 1, na.rm = TRUE),
    Down_datasets = rowSums(sign(zmat[meta_candidates, ]) == -1, na.rm = TRUE)
  )
  meta_df <- meta_df[order(meta_df$RRA_score), ]
  write.csv(meta_df, file = file.path(tissue_dir, "Meta_DEGs_final.csv"), row.names = FALSE)
  
  # Save gene list for GSEA (background = all genes in Z-matrix)
  write.table(rownames(zmat), file = file.path(tissue_dir, "ShinyGO_background.txt"),
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # ========== VISUALISATIONS ==========
  
  
  # 1. Volcano plot (Weighted_Z vs -log10(RRA_score))
  meta_df$negLogRRA <- -log10(meta_df$RRA_score + 1e-10)
  p_vol <- ggplot(meta_df, aes(x = Weighted_Z, y = negLogRRA, color = Direction)) +
    geom_point(alpha = 0.7) + theme_bw() +
    geom_text_repel(data = head(meta_df, 20), aes(label = Gene), size = 3) +
    labs(title = paste(tissue_name, "- Volcano Plot (RRA vs Weighted Z)"),
         x = "Weighted Z", y = "-log10(RRA score)")
  ggsave(file.path(tissue_dir, "Volcano_Plot.png"), p_vol, width = 10, height = 7)
  
  # 2. Boxplot comparing logFC between datasets (for top 50 meta-DEGs)
  if (nrow(meta_df) >= 2 && length(dataset_list) >= 2) {
    # Extract logFC for top 50 genes across all datasets
    top_genes <- meta_df$Gene[1:min(50, nrow(meta_df))]
    logFC_matrix <- sapply(fc_list, function(fc) fc[top_genes])
    rownames(logFC_matrix) <- top_genes
    # Melt for ggplot
    logFC_long <- melt(logFC_matrix, varnames = c("Gene", "Study"), value.name = "logFC")
 
    #Boxplot of absolute logFC across studies for top genes
    p_box <- ggplot(logFC_long, aes(x = Study, y = abs(logFC), fill = Study)) +
      geom_boxplot() + theme_minimal() +
      labs(title = paste(tissue_name, "- Absolute logFC distribution (top 50 meta-DEGs)"))
    ggsave(file.path(tissue_dir, "Boxplot_abslogFC_top50.png"), p_box, width = 8, height = 6)
  }
  
  # 3. Heatmap of scaled Z-scores for top 50 meta-DEGs
  if (nrow(meta_df) > 0) {
    top_heat <- meta_df$Gene[1:min(50, nrow(meta_df))]
    heat_mat <- zmat_scaled[top_heat, , drop = FALSE]
    heat_mat <- heat_mat[complete.cases(heat_mat), ]
    if (nrow(heat_mat) > 1) {
      pheatmap(heat_mat, cluster_cols = FALSE, scale = "none",
               color = colorRampPalette(rev(brewer.pal(7, "RdYlBu")))(100),
               main = paste(tissue_name, "- Top 50 Meta-DEGs (scaled Z)"),
               filename = file.path(tissue_dir, "Heatmap_Top50.png"),
               width = 10, height = 12)
    }
  }
  
  cat("  Meta-analysis for", tissue_name, "completed.\n")
  return(meta_df)
}

# ----------------------------- 4. CROSS-TISSUE COMPARISON ------------------
cross_tissue_analysis <- function(out_dir, blood_meta, brain_meta) {
  common <- merge(blood_meta, brain_meta, by = "Gene", suffixes = c("_Blood", "_Brain"))
  common$Status <- "Discordant"
  common$Status[common$Weighted_Z_Blood > 0 & common$Weighted_Z_Brain > 0] <- "Concordant Up"
  common$Status[common$Weighted_Z_Blood < 0 & common$Weighted_Z_Brain < 0] <- "Concordant Down"
  
  # Save shared signature
  write.csv(common[, c("Gene", "Status", "Weighted_Z_Blood", "Weighted_Z_Brain")],
            file = file.path(out_dir, "Cross_Tissue_Shared_Signature.csv"), row.names = FALSE)
  
  # Scatter plot of weighted Z (blood vs brain)
  p_scatter <- ggplot(common, aes(x = Weighted_Z_Blood, y = Weighted_Z_Brain, color = Status)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_text_repel(data = subset(common, Status != "Discordant")[1:20, ], aes(label = Gene), size = 3) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("Concordant Up" = "#D73027", "Concordant Down" = "#4575B4",
                                  "Discordant" = "grey80")) +
    labs(x = "Weighted Z (Blood)", y = "Weighted Z (Brain)",
         title = "Cross-Tissue Concordance: Blood vs Brain") +
    theme_minimal()
  ggsave(file.path(out_dir, "Cross_Tissue_Scatter.png"), p_scatter, width = 8, height = 7)
  
  
  # Venn diagram
  p_venn <- ggVennDiagram(list(Blood = blood_meta$Gene, Brain = brain_meta$Gene))
  ggsave(file.path(out_dir, "Cross_Tissue_Venn.png"), p_venn, width = 6, height = 6)
}

# ----------------------------- 5. RUN ANALYSIS -----------------------------
blood_results <- run_meta_analysis(blood_datasets, data_dir, out_dir, "Blood")
brain_results <- run_meta_analysis(brain_datasets, data_dir, out_dir, "Brain")

cross_tissue_analysis(out_dir, blood_results, brain_results)

cat("\n✅ Meta-analysis completed. Results saved in:\n", out_dir, "\n")
# ============================================================================
# CONTINUATION: GSEA (weighted meta-Z) + WGCNA module correlation
# Run AFTER your base meta-analysis script
# ============================================================================

# ----------------------------- 0. Load required packages -------------------
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
pkgs <- c("clusterProfiler", "org.Hs.eg.db", "WGCNA", "ggplot2", "reshape2")
for(p in pkgs) {
  if (!require(p, character.only = TRUE, quietly = TRUE))
    BiocManager::install(p, ask = FALSE)
  library(p, character.only = TRUE)
}
enableWGCNAThreads()

# ----------------------------- 1. Paths (same as in base script) ----------
data_dir <- "C:/Users/User/Desktop/Parkinson's Disease-Feb2026/"
out_dir <- file.path(data_dir, "PD_MetaAnalysis_Results")
blood_datasets <- c("GSE75249","GSE6613","GSE49126","GSE161199")
brain_datasets <- c("GSE42966","GSE19587","GSE7621","GSE8397A","GSE8397B")

# Load meta-DEG results
blood_meta <- read.csv(file.path(out_dir, "Blood/Meta_DEGs_final.csv"), stringsAsFactors = FALSE)
brain_meta <- read.csv(file.path(out_dir, "Brain/Meta_DEGs_final.csv"), stringsAsFactors = FALSE)

run_gsea_weighted <- function(tissue_name, dataset_list, seed = 123, nPerm = 10000) {
  set.seed(seed)
  cat("\n========== GSEA (weighted meta-Z) for", tissue_name, "==========\n")
  
  # 1. Set random seed for reproducibility
  
  
  zmat_file <- file.path(out_dir, tissue_name, "Z_matrix_raw.csv")
  if(!file.exists(zmat_file)) {
    cat("Z_matrix_raw.csv not found for", tissue_name, "\n")
    return(NULL)
  }
  zmat <- read.csv(zmat_file, row.names = 1, check.names = FALSE)
  
  # Sample sizes
  sample_sizes <- c()
  for(gse in dataset_list) {
    pheno_file <- file.path(data_dir, paste0(gse, "_pheno.csv"))
    if(file.exists(pheno_file)) {
      pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)
      sample_sizes[gse] <- nrow(pheno)
    }
  }
  sample_sizes <- sample_sizes[colnames(zmat)]
  weights <- sqrt(sample_sizes)
  
  # Scale Z-matrix and compute weighted meta-Z
  zmat_scaled <- scale(zmat)
  w_mat <- sweep(matrix(1, nrow(zmat_scaled), ncol(zmat_scaled)), 2, weights, "*")
  w_mat[is.na(zmat_scaled)] <- NA
  meta_Z <- rowSums(zmat_scaled * w_mat, na.rm = TRUE) / rowSums(w_mat, na.rm = TRUE)
  meta_Z <- meta_Z[!is.na(meta_Z)]
  
  # Optional: break ties by adding tiny noise (makes ranking deterministic)
  # Only needed if many genes have exactly the same meta_Z value
  if(anyDuplicated(meta_Z)) {
    meta_Z <- meta_Z + runif(length(meta_Z), min = 0, max = 1e-10)
  }
  
  gene_rank <- sort(meta_Z, decreasing = TRUE)
  
  # Convert to Entrez IDs
  gene_df <- data.frame(SYMBOL = names(gene_rank), rank = gene_rank)
  entrez_map <- bitr(gene_df$SYMBOL, fromType = "SYMBOL", toType = "ENTREZID",
                     OrgDb = org.Hs.eg.db, drop = TRUE)
  gene_rank_entrez <- gene_df$rank[match(entrez_map$SYMBOL, gene_df$SYMBOL)]
  names(gene_rank_entrez) <- entrez_map$ENTREZID
  
  # 2. Run GSEA with many permutations (default 1000 -> 10000)
  gsea_res <- tryCatch(
    gseKEGG(geneList = gene_rank_entrez, organism = "hsa",
            minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05,
            nPermSimple = nPerm,
            seed = seed),
    error = function(e) NULL
  )
  
  if(is.null(gsea_res) || nrow(gsea_res@result) == 0) {
    cat("No significant GSEA results for", tissue_name, "\n")
    return(NULL)
  }
  write.csv(gsea_res@result, file = file.path(out_dir, tissue_name, "GSEA_KEGG_weighted_results.csv"), row.names = FALSE)
  p <- dotplot(gsea_res, showCategory = 15, title = paste(tissue_name, "GSEA (weighted meta-Z)"))
  ggsave(file.path(out_dir, tissue_name, "GSEA_KEGG_weighted_dotplot.png"), p, width = 10, height = 8)
  cat("GSEA completed for", tissue_name, "\n")
  return(gsea_res)
}

# Run GSEA (using the same seed for both tissues)
gsea_blood <- run_gsea_weighted("Blood", blood_datasets, seed = 123, nPerm = 10000)
gsea_brain <- run_gsea_weighted("Brain", brain_datasets, seed = 123, nPerm = 10000)

# ----------------------------- 3. Correlate GSEA leading edge with meta-DEGs
correlate_gsea_meta <- function(tissue_name, meta_df) {
  gsea_file <- file.path(out_dir, tissue_name, "GSEA_KEGG_weighted_results.csv")
  if(!file.exists(gsea_file)) {
    cat("No GSEA results for", tissue_name, "\n")
    return(NULL)
  }
  gsea_res <- read.csv(gsea_file, stringsAsFactors = FALSE)
  gsea_sig <- gsea_res[gsea_res$p.adjust < 0.05, ]
  if(nrow(gsea_sig) == 0) {
    cat("No significant pathways for", tissue_name, "\n")
    return(NULL)
  }
  
  meta_entrez <- bitr(meta_df$Gene, fromType = "SYMBOL", toType = "ENTREZID",
                      OrgDb = org.Hs.eg.db, drop = TRUE)$ENTREZID
  
  overlap_stats <- list()
  for(i in 1:nrow(gsea_sig)) {
    leading <- gsea_sig$core_enrichment[i]
    if(is.na(leading)) next
    leading_entrez <- as.character(unlist(strsplit(leading, "/")))
    overlap <- intersect(leading_entrez, meta_entrez)
    overlap_stats[[i]] <- data.frame(
      Pathway = gsea_sig$Description[i],
      NES = gsea_sig$NES[i],
      FDR = gsea_sig$p.adjust[i],
      LeadingEdgeSize = length(leading_entrez),
      MetaDEG_Overlap = length(overlap),
      OverlapPercent = round(length(overlap)/length(leading_entrez)*100, 2)
    )
  }
  overlap_df <- do.call(rbind, overlap_stats)
  write.csv(overlap_df, file = file.path(out_dir, tissue_name, "GSEA_MetaDEG_overlap.csv"), row.names = FALSE)
  
  # Barplot
  ov_top <- head(overlap_df[order(overlap_df$OverlapPercent, decreasing = TRUE), ], 10)
  p <- ggplot(ov_top, aes(x = reorder(Pathway, OverlapPercent), y = OverlapPercent)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(x = "", y = "% of leading edge genes that are meta-DEGs",
         title = paste(tissue_name, "GSEA Leading Edge vs Meta-DEGs")) +
    theme_minimal()
  ggsave(file.path(out_dir, tissue_name, "GSEA_MetaDEG_overlap_barplot.png"), p, width = 8, height = 6)
  cat("Correlation saved for", tissue_name, "\n")
}

correlate_gsea_meta("Blood", blood_meta)
correlate_gsea_meta("Brain", brain_meta)