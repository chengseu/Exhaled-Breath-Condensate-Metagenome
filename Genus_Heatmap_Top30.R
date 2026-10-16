## ============================================================
## 0. Load required R packages
## ============================================================
pkgs <- c("readxl", "pheatmap", "grDevices")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

## ============================================================
## 1. Read the Z-score standardized matrix (Excel file)
## ============================================================
input_file <- "D:/Desktop/Article/Codes/Result/Genus_Zscore_Matrix.xlsx"
zdf <- readxl::read_excel(input_file)

genus_colname <- colnames(zdf)[1]
genus_vec <- zdf[[genus_colname]]
zdf_num <- zdf[, -1]

zmat <- as.matrix(zdf_num)
mode(zmat) <- "numeric"
rownames(zmat) <- genus_vec

## ============================================================
## 2. Select the top 30 genera (ranked by row variance)
## ============================================================
row_var <- apply(zmat, 1, var, na.rm = TRUE)
row_var <- row_var[!is.na(row_var)]

top_n <- min(30, length(row_var))   # top 30 genera
top_genus <- names(sort(row_var, decreasing = TRUE))[1:top_n]

zmat_top30 <- zmat[top_genus, , drop = FALSE]

## ============================================================
## 3. Build the group annotation bar (HSC vs MBE)
## ============================================================
sample_names <- colnames(zmat_top30)
group_vec <- ifelse(grepl("HSC", sample_names), "HSC", "MBE")

annotation_col <- data.frame(Group = factor(group_vec))
rownames(annotation_col) <- sample_names

ann_colors <- list(Group = c(HSC = "#ADE8F4", MBE = "#C9E4CA"))

## ============================================================
## 4. Build a blue-white-red gradient palette with a fixed scale of [-2, 2]
## ============================================================
custom_palette <- colorRampPalette(c("#0000FF", "#FFFFFF", "#FF0000"))(100)

z_limit <- 2
zmat_plot <- zmat_top30
zmat_plot[zmat_plot >  z_limit] <-  z_limit
zmat_plot[zmat_plot < -z_limit] <- -z_limit

## ============================================================
## 5. Output directory
## ============================================================
outdir <- "D:/Desktop/Article/Figures"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

## ============================================================
## 6. Screen preview (recommended to check first)
## ============================================================
pheatmap::pheatmap(
  zmat_plot,
  color = custom_palette,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  fontsize_row = 6,
  fontsize_col = 8,
  border_color = NA,
  show_rownames = TRUE,
  show_colnames = TRUE,
  breaks = seq(-2, 2, length.out = 101)
)

## ============================================================
## 7. Save as TIFF / PDF / SVG
## ============================================================

## TIFF
tiff(file.path(outdir, "Genus_Heatmap_Top30.tiff"),
     width = 2400, height = 2000, res = 300)
pheatmap::pheatmap(
  zmat_plot,
  color = custom_palette,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  fontsize_row = 6,
  fontsize_col = 8,
  border_color = NA,
  show_rownames = TRUE,
  show_colnames = TRUE,
  breaks = seq(-2, 2, length.out = 101)
)
dev.off()

## PDF
pdf(file.path(outdir, "Genus_Heatmap_Top30.pdf"),
    width = 10, height = 8)
pheatmap::pheatmap(
  zmat_plot,
  color = custom_palette,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  fontsize_row = 6,
  fontsize_col = 8,
  border_color = NA,
  show_rownames = TRUE,
  show_colnames = TRUE,
  breaks = seq(-2, 2, length.out = 101)
)
dev.off()

## SVG
svg(file.path(outdir, "Genus_Heatmap_Top30.svg"),
    width = 10, height = 8)
pheatmap::pheatmap(
  zmat_plot,
  color = custom_palette,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  scale = "none",
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  fontsize_row = 6,
  fontsize_col = 8,
  border_color = NA,
  show_rownames = TRUE,
  show_colnames = TRUE,
  breaks = seq(-2, 2, length.out = 101)
)
dev.off()

cat("✅ The top 30 heat map has been successfully generated and saved to:, outdir, "\n")
