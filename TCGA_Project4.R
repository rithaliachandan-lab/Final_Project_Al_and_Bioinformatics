
#Set working directory
setwd("C:/Users/Chandan/Downloads/TCGA_LUAD")

untar("tcga_data.tar.gz")
ss <- read.delim("gdc_sample_sheet.tsv", stringsAsFactors = FALSE)

ss_sub <- subset(ss,Sample.Type %in% c("Primary Tumor", "Solid Tissue Normal"))
ss_sub$Dir <- ss_sub$File.ID

# Keep only RNA-seq augmented STAR gene count files
rna_rows <- grepl("rna_seq\\.augmented_star_gene_counts\\.tsv$", ss_sub$File.Name)
ss_rna <- ss_sub[rna_rows, ]

#extract NKX2-1 TPM from a single RNA-seq folder
get_tpm_NKX2 <- function(folder_name) {
  file_path <- list.files(
    folder_name,
    pattern = "rna_seq\\.augmented_star_gene_counts\\.tsv$",
    full.names = TRUE )
  if (length(file_path) == 0) return(NA_real_)
  
# Read the raw lines
  lines <- readLines(file_path[1])
  header_idx <- grep("^gene_id", lines)[1]
  if (is.na(header_idx)) return(NA_real_)
  header_fields <- strsplit(lines[header_idx], "\t", fixed = TRUE)[[1]]
  data_lines <- lines[(header_idx + 1):length(lines)]
  
# Split each data line on tab
  split_data <- strsplit(data_lines, "\t", fixed = TRUE)
  max_len <- max(lengths(split_data))
  
  mat <- t(vapply(
    split_data,
    function(x) c(x, rep(NA_character_, max_len - length(x))),
    character(max_len)
  ))
  
  dat <- as.data.frame(mat, stringsAsFactors = FALSE)
  colnames(dat) <- header_fields
  
#need the gene_name column to find NKX2-1
  if (!"gene_name" %in% colnames(dat)) return(NA_real_)
  nkx_row <- dat[dat$gene_name == "NKX2-1", , drop = FALSE]
  if (nrow(nkx_row) == 0) return(NA_real_)
  as.numeric(nkx_row$tpm_unstranded[1])
}

# Generic version: extract TPM for ANY gene symbol
get_tpm_gene <- function(folder_name, gene_symbol) {
  file_path <- list.files(
    folder_name,
    pattern = "rna_seq\\.augmented_star_gene_counts\\.tsv$",
    full.names = TRUE
  )
  if (length(file_path) == 0) return(NA_real_)
  
  # Read the raw lines
  lines <- readLines(file_path[1])
  header_idx <- grep("^gene_id", lines)[1]
  if (is.na(header_idx)) return(NA_real_)
  header_fields <- strsplit(lines[header_idx], "\t", fixed = TRUE)[[1]]
  data_lines <- lines[(header_idx + 1):length(lines)]
  
  split_data <- strsplit(data_lines, "\t", fixed = TRUE)
  max_len <- max(lengths(split_data))
  
  mat <- t(vapply(
    split_data,
    function(x) c(x, rep(NA_character_, max_len - length(x))),
    character(max_len)
  ))
  
  dat <- as.data.frame(mat, stringsAsFactors = FALSE)
  colnames(dat) <- header_fields
  
  if (!"gene_name" %in% colnames(dat)) return(NA_real_)
  gene_row <- dat[dat$gene_name == gene_symbol, , drop = FALSE]
  if (nrow(gene_row) == 0) return(NA_real_)
  
  as.numeric(gene_row$tpm_unstranded[1])
}



## 4. Apply function to all samples and compute log2(TPM + 1)
ss_rna$TPM_NKX2_1 <- sapply(ss_rna$Dir, get_tpm_NKX2)

ss_clean <- subset(ss_rna, !is.na(TPM_NKX2_1))
ss_clean$log2_TPM_plus_1 <- log2(ss_clean$TPM_NKX2_1 + 1)


## 5. Plot tumor vs normal NKX2-1 expression
boxplot(
  log2_TPM_plus_1 ~ Sample.Type,
  data = ss_clean,
  main = "NKX2-1 Expression in LUAD (TCGA)",
  xlab = "Tissue Type",
  ylab = "log2(TPM + 1)",
  col  = c("tomato", "skyblue")
)


# ================== PART B: MULTI-GENE ANALYSIS + PLOTS ==================

genes <- c("EGFR", "TP53")   # extra genes

for (g in genes) {
  # TPM column name: e.g., TPM_EGFR
  tpm_col <- paste0("TPM_", g)
  # log2 column name: e.g., log2_TPM_EGFR
  log_col <- paste0("log2_TPM_", g)
  
  ss_rna[[tpm_col]] <- sapply(ss_rna$Dir, get_tpm_gene, gene_symbol = g)
  ss_rna[[log_col]] <- log2(ss_rna[[tpm_col]] + 1)
}

# NKX2-1 log2 also in ss_rna for multi-gene plot
ss_rna$log2_TPM_NKX2_1 <- log2(ss_rna$TPM_NKX2_1 + 1)

# Keep rows where all three genes are present
ss_multi <- subset(
  ss_rna,
  !is.na(log2_TPM_NKX2_1) &
    !is.na(log2_TPM_EGFR) &
    !is.na(log2_TPM_TP53)
)

par(mfrow = c(1, 3))

boxplot(
  log2_TPM_NKX2_1 ~ Sample.Type,
  data = ss_multi,
  main = "NKX2-1",
  xlab = "Tissue Type",
  ylab = "log2(TPM + 1)",
  col  = c("green", "pink")
  )

boxplot(
  log2_TPM_EGFR ~ Sample.Type,
  data = ss_multi,
  main = "EGFR",
  xlab = "Tissue Type",
  ylab = "log2(TPM + 1)",
  col  = c("purple", "red")
)

boxplot(
  log2_TPM_TP53 ~ Sample.Type,
  data = ss_multi,
  main = "TP53",
  xlab = "Tissue Type",
  ylab = "log2(TPM + 1)",
  col  = c("yellow", "grey")
)

par(mfrow = c(1, 1))

# ========= STEP: CORRELATION ANALYSIS (TUMOR SAMPLES ONLY) =========

# 1. Keep only Primary Tumor samples
tumor_only <- subset(ss_rna, Sample.Type == "Primary Tumor")

# 2. Remove NA values for all three genes
tumor_only <- subset(
  tumor_only,
  !is.na(log2_TPM_NKX2_1) &
    !is.na(log2_TPM_EGFR) &
    !is.na(log2_TPM_TP53)
)

# 3. Check number of tumor samples
nrow(tumor_only)

# 4. Pearson Correlation Tests

# NKX2-1 vs EGFR
cor_NKX_EGFR <- cor.test(
  tumor_only$log2_TPM_NKX2_1,
  tumor_only$log2_TPM_EGFR,
  method = "pearson"
)
print(cor_NKX_EGFR)

# NKX2-1 vs TP53
cor_NKX_TP53 <- cor.test(
  tumor_only$log2_TPM_NKX2_1,
  tumor_only$log2_TPM_TP53,
  method = "pearson"
)
print(cor_NKX_TP53)

# 5. Scatter Plots
par(mfrow = c(1, 2))

plot(
  tumor_only$log2_TPM_NKX2_1,
  tumor_only$log2_TPM_EGFR,
  xlab = "NKX2-1 (log2 TPM+1)",
  ylab = "EGFR (log2 TPM+1)",
  main = "NKX2-1 vs EGFR (Tumor)"
)

plot(
  tumor_only$log2_TPM_NKX2_1,
  tumor_only$log2_TPM_TP53,
  xlab = "NKX2-1 (log2 TPM+1)",
  ylab = "TP53 (log2 TPM+1)",
  main = "NKX2-1 vs TP53 (Tumor)"
)

par(mfrow = c(1, 1))


# ========== CLEAN HEATMAP (PUBLICATION STYLE) ==========

# Tumor only
tumor_only <- subset(ss_rna, Sample.Type == "Primary Tumor")

tumor_only <- subset(
  tumor_only,
  !is.na(log2_TPM_NKX2_1) &
    !is.na(log2_TPM_EGFR) &
    !is.na(log2_TPM_TP53)
)

# Expression matrix
heatmap_data <- as.matrix(
  tumor_only[, c("log2_TPM_NKX2_1", "log2_TPM_EGFR", "log2_TPM_TP53")]
)

colnames(heatmap_data) <- c("NKX2-1", "EGFR", "TP53")

# Do NOT show row names (this removes congestion)
rownames(heatmap_data) <- NULL

# Better color palette
my_colors <- colorRampPalette(c("navy", "white", "firebrick3"))(100)

# Draw clean heatmap
heatmap(
  heatmap_data,
  Rowv = NA,
  Colv = NA,
  scale = "column",
  col = my_colors,
  labRow = FALSE,   # ✅ hides sample IDs
  main = "Heatmap of LUAD Biomarkers (Tumor Samples)",
  xlab = "Genes",
  ylab = "Tumor Samples"
)

# ========= PCA WITH 3-CLUSTER COLORING =========

# Tumor only
tumor_only <- subset(ss_rna, Sample.Type == "Primary Tumor")

tumor_only <- subset(
  tumor_only,
  !is.na(log2_TPM_NKX2_1) &
    !is.na(log2_TPM_EGFR) &
    !is.na(log2_TPM_TP53)
)

# PCA matrix
pca_data <- as.matrix(
  tumor_only[, c("log2_TPM_NKX2_1",
                 "log2_TPM_EGFR",
                 "log2_TPM_TP53")]
)

# Run PCA
pca_res <- prcomp(pca_data, scale. = TRUE)

# Variance %
pc1_var <- round(summary(pca_res)$importance[2,1] * 100, 2)
pc2_var <- round(summary(pca_res)$importance[2,2] * 100, 2)

# ------------------------
# 3-CLUSTER KMEANS ON PCA
# ------------------------
set.seed(123)
clus <- kmeans(pca_res$x[,1:2], centers = 3)

# Define 3 colors
cluster_colors <- c("red", "blue", "darkgreen")

# CLEAN colored PCA plot
plot(
  pca_res$x[,1],
  pca_res$x[,2],
  pch = 19,
  col = cluster_colors[clus$cluster],
  xlab = paste0("PC1 (", pc1_var, "% variance)"),
  ylab = paste0("PC2 (", pc2_var, "% variance)"),
  main = "PCA of LUAD Tumor Samples (3 Molecular Clusters)"
)

# Add legend
legend(
  "topright",
  legend = c("Cluster 1", "Cluster 2", "Cluster 3"),
  col = cluster_colors,
  pch = 19,
  cex = 0.9
)



saveRDS(ss_rna, file = "ss_rna_for_shiny.rds")
