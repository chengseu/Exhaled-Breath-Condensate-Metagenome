# ---- Analysis-specific configuration ----
ANALYSIS_NAME  <- "3D-W"
K_DIM          <- 3
INCLUDE_MBEM09 <- TRUE
SEED           <- 20261007
TRYMAX         <- 200
MAXIT          <- 1000
FILE_PREFIX    <- "NMDS_3D_W"

IN_FILE  <- "D:/Desktop/Article1/02.Result-X101SC25075182-Z01-J003-meta/result/04.Annotation/MicroNR/Relative/Unigenes.relative.p.xls"
REF_FILE <- "D:/Desktop/Article1/Codes/Data/NMDS.xlsx"
BASE_DIR <- "D:/Desktop/Article1/Figures/Revised/Fig3c"
OUT_DIR  <- file.path(BASE_DIR, ANALYSIS_NAME)
for (d in c("Figure", "Logs", "R", "Tables"))
  dir.create(file.path(OUT_DIR, d), recursive = TRUE, showWarnings = FALSE)

# ---- Packages ----
suppressPackageStartupMessages({
  library(readr); library(dplyr); library(stringr); library(ggplot2)
  library(ggpubr); library(vegan); library(svglite); library(ragg)
  library(readxl); library(scatterplot3d)
})

# ---- Logging ----
log_lines <- character()
log <- function(...) { msg <- paste0(...); log_lines <<- c(log_lines, msg); cat(msg, "\n") }
log("=== NMDS analysis: ", ANALYSIS_NAME, " ===")
log("Analysis start time: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
log("R version: ", R.version.string)
log("Input file: ", IN_FILE)
log("")

# ---- Read input (tab-delimited .xls text; trailing Description column dropped) ----
raw <- readLines(IN_FILE, warn = FALSE)
hdr <- str_trim(str_split(raw[1], "\t")[[1]]); hdr <- hdr[hdr != ""]
taxa_col <- hdr[1]; samples <- hdr[-1]
parsed <- lapply(raw[-1], function(l) { x <- str_trim(str_split(l, "\t")[[1]]); x[seq_along(hdr)] })
mraw <- do.call(rbind, parsed); colnames(mraw) <- hdr
taxon <- mraw[, taxa_col]
mat <- as.data.frame(lapply(as.data.frame(mraw[, samples, drop = FALSE], stringsAsFactors = FALSE),
                            function(z) suppressWarnings(as.numeric(z))))
rownames(mat) <- taxon
comm <- t(mat)   # sample x taxon

# ---- Metadata ----
meta <- data.frame(
  Sample = samples,
  Group  = ifelse(str_detect(samples, "^HSC"), "HSC", "MBE"),
  Gender = ifelse(str_sub(samples, 4, 4) == "M", "Male", "Female"),
  ID     = str_extract(samples, "\\d+"),
  stringsAsFactors = FALSE
)

log("Data loaded: ", ncol(comm), " taxa x ", nrow(comm), " samples.")
log("Distance method: Bray-Curtis (vegdist 'bray'), autotransform = FALSE.")
log("NMDS dimension k = ", K_DIM, ". Seed = ", SEED, ". trymax = ", TRYMAX, ". maxit = ", MAXIT, ".")
log("")

# ---- MBEM09 handling ----
if (!INCLUDE_MBEM09) {
  keep <- rownames(comm) != "MBEM09"
  comm <- comm[keep, , drop = FALSE]
  meta <- meta[keep, , drop = FALSE]
  log("MBEM09 EXCLUDED from the raw abundance matrix BEFORE distance + NMDS.")
} else {
  log("MBEM09 INCLUDED (all samples participate in NMDS).")
}
n_before <- ncol(comm)
comm <- comm[, colSums(comm) > 0, drop = FALSE]   # drop taxa that became all-zero
if (ncol(comm) < n_before) log("Dropped ", n_before - ncol(comm), " all-zero taxa after sample filtering.")
log("Samples used (N = ", nrow(comm), "): ", paste(rownames(comm), collapse = ", "))
log("")

# ---- Distance + NMDS (k = 3) ----
log("Computing Bray-Curtis distance + 3D NMDS ...")
dist_bray <- vegdist(comm, method = "bray")
warn_msgs <- character()
nmds <- withCallingHandlers(
  { set.seed(SEED); metaMDS(dist_bray, k = K_DIM, trymax = TRYMAX, maxit = MAXIT, trace = 0) },
  warning = function(w) { warn_msgs <<- c(warn_msgs, conditionMessage(w)); invokeRestart("muffleWarning") }
)

stress_val     <- nmds$stress
iters_val      <- nmds$iters
repeats_val    <- nmds$converged   # vegan 2.7.5: count of best-solution repeats
converged_flag <- repeats_val >= 1
log("3D NMDS stress:     ", format(stress_val, digits = 8))
log("Converged:           ", if (converged_flag) paste0("Yes (best solution repeated ", repeats_val, " times)") else "No")
log("Iterations (best):   ", iters_val)
log("Random starts:       ", nmds$tries)
if (length(warn_msgs) > 0) {
  log("WARNINGS: ", paste(warn_msgs, collapse = " | "))
} else {
  log("WARNINGS: none")
}
log("")

# ---- Extract scores (3 columns: NMDS1, NMDS2, NMDS3) ----
sc <- as.data.frame(scores(nmds, display = "sites"))
sc$Sample <- rownames(sc)
sc <- left_join(sc, meta, by = "Sample")

# ---- Visual style (original Figure 3c) ----
fill_colors  <- c("HSC" = "#ADE8F4", "MBE" = "#C9E4CA")
point_colors <- c("HSC" = "#48CAE4", "MBE" = "#87BBA2")

stress_label <- sprintf("3D stress = %.3f", stress_val)
ann_x <- max(sc$NMDS1)
ann_y <- min(sc$NMDS2) - 0.04 * diff(range(sc$NMDS2))

# ---- 2D PROJECTION plot (NMDS1 vs NMDS2), explicitly labelled ----
p <- ggplot(sc, aes(x = NMDS1, y = NMDS2, color = Group, fill = Group)) +
  geom_point(size = 3, alpha = 0.9) +
  stat_ellipse(type = "norm", level = 0.95, geom = "polygon", alpha = 0.4, color = NA) +
  stat_ellipse(type = "norm", level = 0.95, geom = "path",
               linetype = "dashed", color = "grey70", linewidth = 0.6) +
  scale_color_manual(values = point_colors) +
  scale_fill_manual(values = fill_colors) +
  labs(x = "NMDS1", y = "NMDS2",
       subtitle = "3D NMDS projection (NMDS1\u2013NMDS2) \u2014 not a re-computed 2D NMDS") +
  annotate("text", x = ann_x, y = ann_y, label = stress_label,
           size = 4.2, hjust = 1, vjust = 0, color = "black") +
  theme_minimal(base_size = 14, base_family = "Arial") +
  theme(
    legend.position   = "right",
    plot.subtitle     = element_text(size = 10, face = "italic"),
    axis.title.x      = element_text(size = 14),
    axis.title.y      = element_text(size = 14),
    axis.text         = element_text(size = 11, color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks        = element_line(color = "black"),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.75),
    panel.grid.major  = element_line(color = "grey80", linetype = "dashed", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    plot.margin       = margin(t = 10, r = 10, b = 10, l = 10)
  )

# ---- Export projection: TIFF (600 dpi) / PDF / SVG ----
figdir <- file.path(OUT_DIR, "Figure")
tiff_path <- file.path(figdir, paste0(FILE_PREFIX, "_projection.tiff"))
pdf_path  <- file.path(figdir, paste0(FILE_PREFIX, "_projection.pdf"))
svg_path  <- file.path(figdir, paste0(FILE_PREFIX, "_projection.svg"))

ragg::agg_tiff(tiff_path, width = 6, height = 5, units = "in", res = 600)
print(p); dev.off()
cairo_pdf(pdf_path, width = 6, height = 5, family = "Arial")
print(p); dev.off()
svglite::svglite(svg_path, width = 6, height = 5)
print(p); dev.off()

log("Projection figures saved: ", tiff_path, " / ", pdf_path, " / ", svg_path)

# ---- Supplementary TRUE 3D scatter (scatterplot3d) ----
col3d <- ifelse(sc$Group == "HSC", "#48CAE4", "#87BBA2")
s3d_pdf  <- file.path(figdir, paste0(FILE_PREFIX, "_scatter3d.pdf"))
s3d_tiff <- file.path(figdir, paste0(FILE_PREFIX, "_scatter3d.tiff"))

cairo_pdf(s3d_pdf, width = 6, height = 5, family = "Arial")
s3d <- scatterplot3d(sc$NMDS1, sc$NMDS2, sc$NMDS3,
                     color = col3d, pch = 16, cex.symbols = 1.1,
                     xlab = "NMDS1", ylab = "NMDS2", zlab = "NMDS3",
                     angle = 45, main = paste0("3D NMDS (", stress_label, ")"))
legend("topright", legend = c("HSC", "MBE"), col = c("#48CAE4", "#87BBA2"),
       pch = 16, bty = "n", cex = 0.9)
dev.off()

ragg::agg_tiff(s3d_tiff, width = 6, height = 5, units = "in", res = 600)
s3d <- scatterplot3d(sc$NMDS1, sc$NMDS2, sc$NMDS3,
                     color = col3d, pch = 16, cex.symbols = 1.1,
                     xlab = "NMDS1", ylab = "NMDS2", zlab = "NMDS3",
                     angle = 45, main = paste0("3D NMDS (", stress_label, ")"))
legend("topright", legend = c("HSC", "MBE"), col = c("#48CAE4", "#87BBA2"),
       pch = 16, bty = "n", cex = 0.9)
dev.off()

log("3D scatter figures saved: ", s3d_pdf, " / ", s3d_tiff)

# ---- Save 3D coordinates + stress stats ----
write.csv(sc, file.path(OUT_DIR, "Tables", paste0(FILE_PREFIX, "_scores.csv")), row.names = FALSE)
stats <- data.frame(
  Analysis = ANALYSIS_NAME, Dimension = K_DIM,
  MBEM09 = if (INCLUDE_MBEM09) "Included" else "Excluded",
  N_samples = nrow(comm), N_taxa = ncol(comm),
  Stress = stress_val, Distance = "Bray-Curtis", Seed = SEED,
  Iterations = iters_val, RandomStarts = nmds$tries, BestRepeated = repeats_val,
  Convergence = converged_flag,
  Transformation = "None (relative abundance, autotransform=FALSE)",
  stringsAsFactors = FALSE
)
write.csv(stats, file.path(OUT_DIR, "Tables", paste0(FILE_PREFIX, "_stress.csv")), row.names = FALSE)
log("Saved: ", FILE_PREFIX, "_scores.csv and ", FILE_PREFIX, "_stress.csv")

# ---- Write log ----
log("")
log("Output directory: ", OUT_DIR)
log("Analysis finished: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
writeLines(log_lines, file.path(OUT_DIR, "Logs", paste0(FILE_PREFIX, ".log")))
cat("DONE ", ANALYSIS_NAME, " | 3D stress = ", stress_val, " | converged = ", converged_flag, " (repeated ", repeats_val, "x)\n")
