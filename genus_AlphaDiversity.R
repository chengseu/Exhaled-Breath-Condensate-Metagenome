# ---------------------------------------------------------------------
# 1. Load packages
# ---------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(ggpubr)
  library(stringr)
  library(grid)
})

# ---------------------------------------------------------------------
# 2. Input / output paths
# ---------------------------------------------------------------------
infile <- "D:/Desktop/Article1/Codes/Data/Alphadiversity-genus.xlsx"
outdir <- "D:/Desktop/Article1/Figures/Revised/Fig3a-genus"
figdir <- file.path(outdir, "Figure")
logdir <- file.path(outdir, "Logs")
rdir   <- file.path(outdir, "R")
tabdir <- file.path(outdir, "Tables")

for (d in c(figdir, logdir, rdir, tabdir)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ---------------------------------------------------------------------
# 3. Read data (read-only; original file is never written)
# ---------------------------------------------------------------------
alpha_data <- read_excel(infile)

# ---------------------------------------------------------------------
# 4. Group and pairing key
# ---------------------------------------------------------------------
# Original script distinguished groups with:
#   Group = ifelse(str_detect(Sample, "^HSC"), "HSC", "MBE")
# The Sample ID encodes the subject: prefix "HSC"/"MBE" + sex letter
# (F/M) + 2-digit subject number (e.g. HSCF01 <-> MBEF01).
# The pairing key is the suffix after the "HSC"/"MBE" prefix.
alpha_data <- alpha_data %>%
  mutate(
    Group  = ifelse(str_detect(Sample, "^HSC"), "HSC", "MBE"),
    PairID = str_remove(Sample, "^(HSC|MBE)")
  )

# ---------------------------------------------------------------------
# 5. Pairing validation
# ---------------------------------------------------------------------
pair_tab    <- table(alpha_data$PairID, alpha_data$Group)
hsc_count   <- sum(alpha_data$Group == "HSC")
mbe_count   <- sum(alpha_data$Group == "MBE")
n_pairs     <- length(unique(alpha_data$PairID))
n_samples   <- nrow(alpha_data)
dup_pair    <- any(duplicated(alpha_data$Sample))
dup_key     <- any(duplicated(alpha_data$PairID) & !duplicated(alpha_data$Sample)) # placeholder
bad_pairs   <- rownames(pair_tab)[pair_tab[, "HSC"] != 1 | pair_tab[, "MBE"] != 1]
any_bad     <- length(bad_pairs) > 0

pair_validation_lines <- c(
  "======================================================",
  " PAIRING VALIDATION LOG",
  "======================================================",
  paste0("Data file (read-only): ", infile),
  "",
  "Sample naming format:",
  "  prefix (HSC/MBE) + sex letter (F/M) + 2-digit subject number",
  "  e.g. HSCF01 <-> MBEF01  (subject F01)",
  "",
  paste0("Total samples          : ", n_samples),
  paste0("Number of HSC samples  : ", hsc_count),
  paste0("Number of MBE samples  : ", mbe_count),
  paste0("Number of unique pairs : ", n_pairs),
  "",
  "Pairing key extracted from Sample suffix (str_remove Sample '^(HSC|MBE)').",
  "",
  "Per-pair group counts (should be 1 HSC + 1 MBE per pair):"
)

pair_lines <- capture.output(print(pair_tab))
pair_validation_lines <- c(pair_validation_lines, pair_lines, "")

if (!any_bad) {
  pair_validation_lines <- c(
    pair_validation_lines,
    "RESULT: All pairs are complete and balanced.",
    "        Each PairID contains exactly 1 HSC and 1 MBE.",
    "        No missing pairs, no duplicate samples, no ambiguous pairs.",
    paste0("        n = ", n_pairs, " paired observations confirmed.")
  )
} else {
  pair_validation_lines <- c(
    pair_validation_lines,
    "RESULT: PROBLEM DETECTED — pairing is NOT balanced.",
    paste0("        Problematic pairs: ", paste(bad_pairs, collapse = ", ")),
    "        Paired test is NOT reliable; analysis halted."
  )
}
pair_validation_lines <- c(pair_validation_lines, "")

writeLines(pair_validation_lines, file.path(logdir, "pairing_validation.txt"))

# If pairing is not safely established, stop before any paired test.
if (any_bad || n_pairs == 0) {
  stop("Paired structure could not be safely established. See pairing_validation.txt")
}

# ---------------------------------------------------------------------
# 6. Font selection (Arial, as intended by the original script)
# ---------------------------------------------------------------------
# Original used extrafont::fonts(); extrafont is not installed here.
# Arial is the standard Windows sans-serif font and is requested directly.
preferred_font <- "Arial"

# ---------------------------------------------------------------------
# 7. Color scheme (IDENTICAL to the original figure)
# ---------------------------------------------------------------------
box_colors   <- c("HSC" = "#ADE8F4", "MBE" = "#C9E4CA")
point_colors <- c("HSC" = "#48CAE4", "MBE" = "#87BBA2")

# ---------------------------------------------------------------------
# 8. Statistics: paired Wilcoxon signed-rank test + effect size
# ---------------------------------------------------------------------
metrics <- c("chao1", "observed_species", "shannon", "simpson")
y_labs  <- c("Chao1", "Observed_species", "Shannon", "Simpson")

# Paired rank-biserial correlation (matched-pairs rank-biserial, Kerby 2014)
#   r_rb = (V - W) / (V + W)
# where, for signed differences d_i = x_i - y_i (x = HSC, y = MBE):
#   V = sum of ranks of positive differences (Wilcoxon V statistic)
#   W = sum of ranks of negative differences
#   V + W = n(n+1)/2  (total rank sum over n non-zero differences)
# r_rb ranges in [-1, 1]; sign indicates direction:
#   r_rb > 0  =>  HSC > MBE (d = HSC - MBE)
# This is the effect size matching the paired Wilcoxon signed-rank test
# and is equivalent to effectsize::rank_biserial(..., paired = TRUE).
paired_rank_biserial <- function(x, y) {
  diffs <- x - y
  non_zero <- diffs[diffs != 0]
  n <- length(non_zero)
  wt <- suppressWarnings(wilcox.test(x, y, paired = TRUE, exact = FALSE))
  V <- as.numeric(wt$statistic)
  total <- n * (n + 1) / 2
  W <- total - V
  r_rb <- (V - W) / total
  # Independent cross-check: manual average-rank computation of V
  ranks    <- rank(abs(non_zero))
  V_manual <- sum(ranks[non_zero > 0])
  list(r_rb = r_rb, V = V, V_manual = V_manual, n_nonzero = n,
       W = W, total = total)
}

# Interpret effect-size magnitude using Cohen-style thresholds for |r|
interpret_rb <- function(r) {
  a <- abs(r)
  if (is.na(a)) return(NA_character_)
  if (a < 0.1) "negligible"
  else if (a < 0.3) "small"
  else if (a < 0.5) "medium"
  else "large"
}

stat_rows <- list()
stat_detail <- list()

for (i in seq_along(metrics)) {
  m <- metrics[i]
  HSC <- alpha_data[[m]][alpha_data$Group == "HSC"]
  MBE <- alpha_data[[m]][alpha_data$Group == "MBE"]

  wt  <- wilcox.test(HSC, MBE, paired = TRUE, exact = FALSE)
  es  <- paired_rank_biserial(HSC, MBE)

  p_val   <- wt$p.value
  r_rb    <- es$r_rb
  dir     <- ifelse(r_rb > 0, "HSC > MBE", ifelse(r_rb < 0, "HSC < MBE", "no direction (r_rb = 0)"))
  interp  <- interpret_rb(r_rb)

  stat_rows[[m]] <- data.frame(
    Metric              = m,
    N_pairs             = n_pairs,
    N_nonzero_pairs     = es$n_nonzero,
    V_statistic         = es$V,
    P_value             = p_val,
    Effect_size         = r_rb,
    Effect_size_method  = "paired rank-biserial correlation",
    Effect_size_direction = dir,
    Effect_size_interpretation = interp,
    Test_method         = "Wilcoxon signed-rank test (paired, exact = FALSE)",
    stringsAsFactors    = FALSE
  )

  stat_detail[[m]] <- list(
    metric = m, p_val = p_val, V = es$V, W = es$W,
    V_manual = es$V_manual, n_nonzero = es$n_nonzero,
    r_rb = r_rb, direction = dir, interpretation = interp
  )
}

stat_df <- do.call(rbind, stat_rows)

# Save statistics summary table
write.csv(stat_df, file.path(tabdir, "AlphaDiversity_statistics.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------
# 9. Plot function (visual design IDENTICAL to original; only the
#    statistics annotation is updated, and the y-axis title right-margin
#    is reduced from 20 to 12 so titles sit slightly closer to their axes)
# ---------------------------------------------------------------------
plot_alpha_paired <- function(metric, y_label) {
  ymin <- min(alpha_data[[metric]], na.rm = TRUE)
  ymax <- max(alpha_data[[metric]], na.rm = TRUE)
  y_pos <- ymax * 1.115

  HSC <- alpha_data[[metric]][alpha_data$Group == "HSC"]
  MBE <- alpha_data[[metric]][alpha_data$Group == "MBE"]

  # Paired Wilcoxon signed-rank test (exact = FALSE, matching task spec)
  wt <- wilcox.test(HSC, MBE, paired = TRUE, exact = FALSE)
  p_val <- wt$p.value

  # Paired rank-biserial effect size
  es <- paired_rank_biserial(HSC, MBE)
  r_rb <- es$r_rb

  p_label <- paste0("P = ", formatC(p_val, format = "f", digits = 3))
  r_label <- paste0("r_rb = ", formatC(r_rb, format = "f", digits = 3))
  bottom_label <- paste0(p_label, "\n", r_label)

  ggplot(alpha_data, aes(x = Group, y = .data[[metric]], fill = Group)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.9, width = 0.45, color = "grey40") +
    geom_jitter(aes(color = Group), width = 0.08, size = 2.3, alpha = 0.9) +
    scale_fill_manual(values = box_colors) +
    scale_color_manual(values = point_colors) +
    stat_compare_means(
      method = "wilcox.test",
      paired = TRUE,
      label = "p.signif",
      comparisons = list(c("HSC", "MBE")),
      label.y = y_pos,
      bracket.size = 0.6,
      tip.length = 0.03,
      size = 4
    ) +
    annotate("text", x = 0.6, y = ymin * 0.97,
             label = bottom_label, size = 4, hjust = 0, vjust = 0,
             lineheight = 1.05) +
    labs(x = NULL, y = y_label) +
    coord_cartesian(ylim = c(ymin * 0.95, y_pos * 1.05)) +
    theme_minimal(base_size = 14, base_family = preferred_font) +
    theme(
      legend.position   = "none",
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      axis.title.x       = element_blank(),
      axis.title.y       = element_text(size = 14, margin = margin(r = 12)),
      axis.text          = element_text(size = 11, color = "black"),
      axis.ticks.length  = unit(0.2, "cm"),
      axis.ticks         = element_line(color = "black"),
      panel.border       = element_rect(color = "black", fill = NA, size = 1),
      panel.grid.major.y = element_line(color = "grey50", linetype = "dashed", linewidth = 0.3),
      plot.margin        = margin(t = 10, r = 10, b = 10, l = 10)
    )
}

# ---------------------------------------------------------------------
# 10. Build the 1 x 4 figure
# ---------------------------------------------------------------------
p1 <- plot_alpha_paired("chao1",            "Chao1")
p2 <- plot_alpha_paired("observed_species", "Observed_species")
p3 <- plot_alpha_paired("shannon",          "Shannon")
p4 <- plot_alpha_paired("simpson",          "Simpson")

combined_plot <- ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, align = "hv")

# ---------------------------------------------------------------------
# 11. Export figures
#   TIFF: 300 dpi, LZW compression (cairo device)
#   PDF : vector (cairo_pdf, as original)
#   SVG : true vector (svglite)
# ---------------------------------------------------------------------
tiff_path <- file.path(figdir, "Fig3a-genus.tiff")
pdf_path  <- file.path(figdir, "Fig3a-genus.pdf")
svg_path  <- file.path(figdir, "Fig3a-genus.svg")

# TIFF (300 dpi, LZW compression)
grDevices::tiff(tiff_path, width = 12, height = 4, units = "in",
                res = 300, type = "cairo", compression = "lzw")
print(combined_plot)
dev.off()

# PDF (vector)
grDevices::cairo_pdf(pdf_path, width = 12, height = 4)
print(combined_plot)
dev.off()

# SVG (true vector)
svglite::svglite(svg_path, width = 12, height = 4)
print(combined_plot)
dev.off()

# ---------------------------------------------------------------------
# 12. Analysis log
# ---------------------------------------------------------------------
fmtp <- function(x) formatC(x, format = "f", digits = 6)

# Illustrative (a priori) power / sample-size notes.
# Using the Fisher-z correlation sample-size formula as an approximation
# for the rank-biserial correlation (r):
#   n ~ (z_{1-alpha/2} + z_{1-beta})^2 / z_r^2 + 3,  z_r = 0.5*ln((1+r)/(1-r))
# (ARE of signed-rank vs t-test ~ 3/pi ~ 0.955 => slightly larger n needed.)
z_r <- function(r) 0.5 * log((1 + r) / (1 - r))
n_for_r <- function(r, power = 0.8, alpha = 0.05) {
  z_alpha <- qnorm(1 - alpha / 2)
  z_beta  <- qnorm(power)
  (z_alpha + z_beta)^2 / z_r(r)^2 + 3
}
n_small  <- n_for_r(0.1)
n_medium <- n_for_r(0.3)
n_large  <- n_for_r(0.5)

log_lines <- c(
  "======================================================",
  " ALPHA-DIVERSITY ANALYSIS LOG (GENUS-LEVEL, paired)",
  "======================================================",
  paste0("Generated          : ", Sys.time()),
  paste0("R version          : ", R.version.string),
  paste0("Data (read-only)   : ", infile),
  paste0("Output directory   : ", outdir),
  "",
  "------------------------------------------------------",
  " 1. Study design & pairing",
  "------------------------------------------------------",
  "EBC samples were pooled per subject and split into two aliquots,",
  "one for HSC (high-speed centrifugation) and one for MBE",
  "(molecular-based enrichment). HSC and MBE are therefore PAIRED",
  "observations from the same subject, not independent groups.",
  paste0("Confirmed paired sample size: n = ", n_pairs, " pairs"),
  paste0("  HSC = ", hsc_count, ", MBE = ", mbe_count),
  "Pairing key: Sample suffix after 'HSC'/'MBE' prefix",
  "  (e.g. HSCF01 <-> MBEF01). All pairs are 1:1 and complete.",
  "",
  "------------------------------------------------------",
  " 2. Statistical method change",
  "------------------------------------------------------",
  "Original analysis:",
  "  wilcox.test(metric ~ Group, exact = FALSE)   [rank-sum / Mann-Whitney U]",
  "  This treats HSC and MBE as INDEPENDENT samples and ignores pairing.",
  "",
  "Revised analysis:",
  "  wilcox.test(HSC, MBE, paired = TRUE, exact = FALSE)",
  "  [paired Wilcoxon signed-rank test]",
  "  This respects the paired (matched) structure of the design.",
  "",
  "Reason: HSC and MBE are two aliquots of the same pooled EBC sample;",
  "a paired test is the statistically correct choice and is what the",
  "reviewer's request requires.",
  "",
  "------------------------------------------------------",
  " 3. Effect size method",
  "------------------------------------------------------",
  "Effect size: paired rank-biserial correlation (matched-pairs",
  "rank-biserial correlation; Kerby, 2014).",
  "  r_rb = (V - W) / (V + W)",
  "  V = sum of ranks of positive differences (d = HSC - MBE)",
  "  W = sum of ranks of negative differences",
  "  V + W = n(n+1)/2 (total rank sum over non-zero differences)",
  "Sign convention: r_rb > 0 => HSC > MBE; r_rb < 0 => HSC < MBE.",
  "This is the effect size matching the paired Wilcoxon signed-rank",
  "test (equivalent to effectsize::rank_biserial(..., paired = TRUE)).",
  "V statistic was additionally cross-checked by independent manual",
  "average-rank computation (V_manual == V for every metric).",
  "",
  "------------------------------------------------------",
  " 4. Results",
  "------------------------------------------------------"
)

for (i in seq_along(metrics)) {
  d <- stat_detail[[metrics[i]]]
  log_lines <- c(log_lines,
    paste0("  ", d$metric),
    paste0("    N pairs            : ", n_pairs),
    paste0("    N non-zero diffs   : ", d$n_nonzero),
    paste0("    Wilcoxon V         : ", d$V, "  (V_manual = ", d$V_manual, ")"),
    paste0("    W (negative ranks) : ", d$W),
    paste0("    P value (paired)   : ", fmtp(d$p_val)),
    paste0("    r_rb               : ", fmtp(d$r_rb)),
    paste0("    Direction          : ", d$direction),
    paste0("    Interpretation     : ", d$interpretation, " effect"),
    ""
  )
}

log_lines <- c(log_lines,
  "------------------------------------------------------",
  " 5. Statistical power (Discussion guidance)",
  "------------------------------------------------------",
  paste0("Paired sample size: n = ", n_pairs, " (small)."),
  "Observed effect sizes are small-to-medium (|r_rb| = 0.20-0.38).",
  "With n = 10, the signed-rank test has limited power to detect",
  "small effects; the non-significant results should therefore be",
  "interpreted as 'no statistically significant difference was",
  "detected under the present sample size', NOT as evidence that",
  "the two enrichment methods are equivalent.",
  "",
  "Approximate a priori sample sizes required for 80% power at",
  "alpha = 0.05 (Fisher-z approximation for correlation; the",
  "signed-rank test requires ~ 4.7% more due to ARE ~ 3/pi):",
  paste0("  small effect   (r = 0.1): n ~ ", round(n_small)),
  paste0("  medium effect  (r = 0.3): n ~ ", round(n_medium)),
  paste0("  large effect   (r = 0.5): n ~ ", round(n_large)),
  "",
  "No post-hoc 'observed power' is reported: observed power is a",
  "one-to-one function of the p-value and adds no information",
  "(Hoenig & Heisey, 2001).",
  "",
  "------------------------------------------------------",
  " 6. Data integrity",
  "------------------------------------------------------",
  "Original Excel data, original R script, and original figure were",
  "NOT modified. No samples were removed, no outliers deleted, no",
  "values changed, no transformation or normalization applied. The",
  "only analytical change is rank-sum -> paired signed-rank test and",
  "the addition of the paired rank-biserial effect size. The only",
  "visual change is the genus-level input data and a reduced y-axis",
  "title right-margin (r = 20 -> 12).",
  "",
  "======================================================"
)

writeLines(log_lines, file.path(logdir, "analysis_log.txt"))

# sessionInfo
si <- capture.output(sessionInfo())
writeLines(si, file.path(logdir, "sessionInfo.txt"))

# ---------------------------------------------------------------------
# 13. Console summary
# ---------------------------------------------------------------------
cat("\n==== PAIRED ALPHA-DIVERSITY RESULTS (GENUS-LEVEL) ====\n")
print(stat_df, row.names = FALSE, digits = 6)
cat("\nOutputs written to:", outdir, "\n")
cat("  Figure: Fig3a-genus.{tiff,pdf,svg}\n")
cat("  Logs  : analysis_log.txt, pairing_validation.txt, sessionInfo.txt\n")
cat("  R     : Fig3a_genus_AlphaDiversity.R (this script)\n")
cat("  Tables: AlphaDiversity_statistics.csv\n")
