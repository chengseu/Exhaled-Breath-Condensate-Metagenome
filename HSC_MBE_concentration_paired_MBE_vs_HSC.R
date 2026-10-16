# ---------------------------------------------------------------------
# 1. Load packages
# ---------------------------------------------------------------------
suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(grid)
})

# ---------------------------------------------------------------------
# 2. Input / output paths
# ---------------------------------------------------------------------
infile <- "D:/Desktop/Article1/Codes/Data/HSC-MBE concentration.xlsx"
outdir <- "D:/Desktop/Article1/Figures/Revised/Fig2c"
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
raw <- read_excel(infile, sheet = 1) %>%
  rename(Sample = 1, Concentration = 2) %>%
  mutate(
    Group  = ifelse(str_detect(Sample, "^HSC"), "HSC", "MBE"),
    Group  = factor(Group, levels = c("HSC", "MBE")),
    PairID = str_remove(Sample, "^(HSC|MBE)")
  )

# ---------------------------------------------------------------------
# 4. Pairing validation
# ---------------------------------------------------------------------
pair_tab    <- table(raw$PairID, raw$Group)
hsc_count   <- sum(raw$Group == "HSC")
mbe_count   <- sum(raw$Group == "MBE")
n_pairs     <- length(unique(raw$PairID))
n_samples   <- nrow(raw)
bad_pairs   <- rownames(pair_tab)[pair_tab[, "HSC"] != 1 | pair_tab[, "MBE"] != 1]
any_bad     <- length(bad_pairs) > 0

pair_lines <- c(
  "======================================================",
  " PAIRING VALIDATION LOG  (HSC vs MBE DNA concentration)",
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
pair_lines <- c(pair_lines, capture.output(print(pair_tab)), "")

if (!any_bad) {
  pair_lines <- c(
    pair_lines,
    "RESULT: All pairs are complete and balanced.",
    "        Each PairID contains exactly 1 HSC and 1 MBE.",
    "        No missing pairs, no duplicate samples, no ambiguous pairs.",
    paste0("        n = ", n_pairs, " paired observations confirmed.")
  )
} else {
  pair_lines <- c(
    pair_lines,
    "RESULT: PROBLEM DETECTED — pairing is NOT balanced.",
    paste0("        Problematic pairs: ", paste(bad_pairs, collapse = ", ")),
    "        Paired test is NOT reliable; analysis halted."
  )
}
pair_lines <- c(pair_lines, "")

writeLines(pair_lines, file.path(logdir, "pairing_validation.txt"))

# If pairing is not safely established, stop before any paired test.
if (any_bad || n_pairs == 0) {
  stop("Paired structure could not be safely established. See pairing_validation.txt")
}

# ---------------------------------------------------------------------
# 5. Align values by PairID (do NOT rely on row order)
# ---------------------------------------------------------------------
hsc_df <- raw %>% filter(Group == "HSC") %>% arrange(PairID)
mbe_df <- raw %>% filter(Group == "MBE") %>% arrange(PairID)
stopifnot(identical(hsc_df$PairID, mbe_df$PairID))

HSC <- hsc_df$Concentration
MBE <- mbe_df$Concentration

# ---------------------------------------------------------------------
# 6. Paired rank-biserial correlation (matched-pairs, Kerby 2014)
#    r_rb = (V - W) / (V + W)
#    Called as paired_rank_biserial(MBE, HSC) so that d = MBE - HSC,
#    i.e. the effect size direction is "MBE relative to HSC",
#    matching Figure 4 and Figure 5 (MBE vs HSC).
#      V = sum of ranks of positive differences (MBE > HSC)
#      W = sum of ranks of negative differences (MBE < HSC)
#    Equivalent to effectsize::rank_biserial(..., paired = TRUE).
# ---------------------------------------------------------------------
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
  list(r_rb = r_rb, V = V, V_manual = V_manual, n_nonzero = n, W = W, total = total)
}

interpret_rb <- function(r) {
  a <- abs(r)
  if (is.na(a)) return(NA_character_)
  if (a < 0.1) "negligible"
  else if (a < 0.3) "small"
  else if (a < 0.5) "medium"
  else "large"
}

# ---------------------------------------------------------------------
# 7. Statistics
# ---------------------------------------------------------------------
# Original rank-sum (for reference only in the log):
orig_test <- wilcox.test(Concentration ~ Group, data = raw, exact = TRUE)
orig_p    <- orig_test$p.value

# Paired Wilcoxon signed-rank test.
# Input order MBE, HSC — reflects the MBE vs HSC contrast direction.
# (For the paired signed-rank test, swapping the two vectors does NOT
#  change the p-value, so P = 0.019 is preserved.)
wt  <- wilcox.test(MBE, HSC, paired = TRUE, exact = FALSE)
es  <- paired_rank_biserial(MBE, HSC)

p_val <- wt$p.value
r_rb  <- es$r_rb
# Direction label: r_rb > 0 => MBE > HSC (since d = MBE - HSC)
dir   <- ifelse(r_rb > 0, "MBE > HSC",
                ifelse(r_rb < 0, "MBE < HSC", "no direction (r_rb = 0)"))
interp <- interpret_rb(r_rb)

# Significance label (same thresholds as original code)
sig_label <- ifelse(p_val < 0.001, "***",
                    ifelse(p_val < 0.01, "**",
                           ifelse(p_val < 0.05, "*", "ns")))

# Statistics summary table
stat_df <- data.frame(
  Metric             = "DNA_concentration",
  N_pairs            = n_pairs,
  N_nonzero_pairs    = es$n_nonzero,
  V_statistic        = es$V,
  P_value            = p_val,
  Effect_size        = r_rb,
  Effect_size_method = "paired rank-biserial correlation",
  Effect_size_direction = dir,
  Effect_size_interpretation = interp,
  Test_method        = "Wilcoxon signed-rank test (paired, exact = FALSE)",
  Contrast           = "MBE vs HSC",
  stringsAsFactors   = FALSE
)
write.csv(stat_df, file.path(tabdir, "HSC_MBE_concentration_MBE_vs_HSC_statistics.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------------
# 8. Plot (visual design IDENTICAL to previous revised figure;
#    only the statistics annotation is moved DOWN)
# ---------------------------------------------------------------------
p_label <- paste0("P = ", formatC(p_val, format = "f", digits = 3))
r_label <- paste0("r_rb = ", formatC(r_rb, format = "f", digits = 3))
bottom_label <- paste0(p_label, "\n", r_label)

p <- ggplot(raw, aes(x = Group, y = Concentration, fill = Group)) +
  geom_boxplot(width = 0.5, outlier.shape = NA, linewidth = 0.4, median.linewidth = 0.4) +
  geom_jitter(width = 0.15, size = 2, shape = 16, color = "grey50", alpha = 0.85) +
  scale_fill_manual(values = c("HSC" = "#CAF0F8", "MBE" = "#FFC2D1"), guide = "none") +
  labs(x = NULL, y = "DNA concentration (ng/μL)") +
  scale_y_continuous(breaks = seq(0, 0.75, 0.1), minor_breaks = seq(0, 0.75, 0.05)) +
  coord_cartesian(ylim = c(0, 0.75), clip = "off") +
  theme_classic(base_size = 14) +
  theme(
    axis.line   = element_blank(),
    axis.ticks  = element_line(color = "black", linewidth = 0.4),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.y = element_text(margin = margin(r = 8)),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4),
    plot.margin = margin(t = 10, r = 20, b = 20, l = 10)
  ) +
  # significance bracket (identical to previous revised figure)
  geom_segment(aes(x = 1, xend = 2, y = 0.70, yend = 0.70), linewidth = 0.3, color = "black") +
  geom_segment(aes(x = 1, xend = 1, y = 0.70, yend = 0.68), linewidth = 0.3, color = "black") +
  geom_segment(aes(x = 2, xend = 2, y = 0.70, yend = 0.68), linewidth = 0.3, color = "black") +
  annotate("text", x = 1.5, y = 0.72, label = sig_label, size = 5) +
  # P value + effect size: moved DOWN (y = 0.05, was 0.12) so the two-line
  # block sits below the HSC boxplot / scatter points and above the x-axis.
  annotate("text", x = 0.55, y = 0.05, label = bottom_label,
           size = 4, hjust = 0, vjust = 0, lineheight = 1.05)

# ---------------------------------------------------------------------
# 9. Export figures (same size as previous revised figure: 2.7 x 4 in)
# ---------------------------------------------------------------------
tiff_path <- file.path(figdir, "HSC_MBE_concentration_MBE_vs_HSC.tiff")
pdf_path  <- file.path(figdir, "HSC_MBE_concentration_MBE_vs_HSC.pdf")
svg_path  <- file.path(figdir, "HSC_MBE_concentration_MBE_vs_HSC.svg")

# TIFF (300 dpi, LZW compression)
grDevices::tiff(tiff_path, width = 2.7, height = 4, units = "in",
                res = 300, type = "cairo", compression = "lzw")
print(p)
dev.off()

# PDF (vector)
grDevices::cairo_pdf(pdf_path, width = 2.7, height = 4)
print(p)
dev.off()

# SVG (true vector)
svglite::svglite(svg_path, width = 2.7, height = 4)
print(p)
dev.off()

# ---------------------------------------------------------------------
# 10. Logs
# ---------------------------------------------------------------------
fmtp <- function(x) formatC(x, format = "f", digits = 6)

# direction_change_log.txt
dir_change_lines <- c(
  "======================================================",
  " DIRECTION CHANGE LOG",
  "======================================================",
  paste0("Generated          : ", Sys.time()),
  "",
  "Scope of this revision (minimal, controlled):",
  "  1. Contrast direction unified from HSC - MBE to MBE - HSC",
  "     (to match Figure 4 and Figure 5, both defined as MBE vs HSC).",
  "  2. The lower-left 'P = ... / r_rb = ...' annotation moved DOWN",
  "     so it no longer overlaps the HSC boxplot / scatter points.",
  "",
  "NOT changed:",
  "  - raw DNA concentration data",
  "  - samples / pairing / grouping",
  "  - statistical test type (paired Wilcoxon signed-rank test)",
  "  - P value (unchanged)",
  "  - effect-size magnitude (unchanged; only sign flips)",
  "  - all Figure visual elements (colors, axes, boxplots, points,",
  "    bracket, significance symbol, size, fonts)",
  "",
  "Original direction : HSC - MBE",
  "Revised direction  : MBE - HSC",
  "",
  "Original displayed effect size : r_rb = -0.855",
  "Revised effect size            : r_rb = +0.855",
  "",
  "P value    : P = 0.019",
  "Test       : Paired Wilcoxon signed-rank test",
  "Effect size: Paired rank-biserial correlation",
  "",
  "Implementation note:",
  "  The direction change is implemented by swapping the input order",
  "  at the function-call level (paired_rank_biserial(MBE, HSC) and",
  "  wilcox.test(MBE, HSC, ...)), NOT by multiplying r_rb by -1.",
  "",
  "======================================================"
)
writeLines(dir_change_lines, file.path(logdir, "direction_change_log.txt"))

# Approximate a priori sample sizes for correlation (Fisher-z)
z_r <- function(r) 0.5 * log((1 + r) / (1 - r))
n_for_r <- function(r, power = 0.8, alpha = 0.05) {
  (qnorm(1 - alpha / 2) + qnorm(power))^2 / z_r(r)^2 + 3
}
n_small  <- n_for_r(0.1)
n_medium <- n_for_r(0.3)
n_large  <- n_for_r(0.5)

log_lines <- c(
  "======================================================",
  " HSC vs MBE DNA CONCENTRATION ANALYSIS LOG (paired, MBE vs HSC)",
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
  " 2. Statistical method",
  "------------------------------------------------------",
  "Original analysis (rank-sum, independent samples):",
  paste0("  wilcox.test(Concentration ~ Group, exact = TRUE)  -> p = ", fmtp(orig_p)),
  "",
  "Revised analysis (paired):",
  "  wilcox.test(MBE, HSC, paired = TRUE, exact = FALSE)",
  "  [paired Wilcoxon signed-rank test]",
  "  The input order MBE, HSC expresses the MBE vs HSC contrast.",
  "",
  "------------------------------------------------------",
  " 3. Effect size",
  "------------------------------------------------------",
  "Effect size: paired rank-biserial correlation (matched-pairs",
  "rank-biserial correlation; Kerby, 2014).",
  "  r_rb = (V - W) / (V + W)",
  "  d = MBE - HSC  (direction: MBE relative to HSC)",
  "  V = sum of ranks of positive differences (MBE > HSC)",
  "  W = sum of ranks of negative differences (MBE < HSC)",
  "Sign convention: r_rb > 0 => MBE > HSC; r_rb < 0 => MBE < HSC.",
  "This matches Figure 4 / Figure 5 direction convention (MBE vs HSC).",
  "Equivalent to effectsize::rank_biserial(..., paired = TRUE).",
  "V statistic additionally cross-checked by manual average-rank",
  "computation (V_manual == V).",
  "",
  "------------------------------------------------------",
  " 4. Results",
  "------------------------------------------------------",
  "  Metric             : DNA concentration",
  "  Contrast           : MBE vs HSC",
  paste0("  N pairs            : ", n_pairs),
  paste0("  N non-zero diffs   : ", es$n_nonzero),
  paste0("  Wilcoxon V         : ", es$V, "  (V_manual = ", es$V_manual, ")"),
  paste0("  W (negative ranks) : ", es$W),
  paste0("  P value (paired)   : ", fmtp(p_val)),
  paste0("  r_rb               : ", fmtp(r_rb)),
  paste0("  Direction          : ", dir),
  paste0("  Interpretation     : ", interp, " effect"),
  "",
  "------------------------------------------------------",
  " 5. Statistical power (Discussion guidance)",
  "------------------------------------------------------",
  paste0("Paired sample size: n = ", n_pairs, " (small)."),
  paste0("Observed effect size: r_rb = ", fmtp(r_rb),
         " (large in magnitude, MBE > HSC)."),
  "Even with a large observed effect, the small sample size limits",
  "precision; confidence in the exact effect magnitude is limited.",
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
  " 6. y-axis label",
  "------------------------------------------------------",
  "Revised  : DNA concentration (ng/μL)",
  "",
  "------------------------------------------------------",
  " 7. Data integrity",
  "------------------------------------------------------",
  "Original Excel data, original R script, and original figure were",
  "NOT modified. No samples were removed, no outliers deleted, no",
  "values changed, no transformation or normalization applied.",
  "",
  "======================================================"
)

writeLines(log_lines, file.path(logdir, "analysis_log.txt"))

# sessionInfo
writeLines(capture.output(sessionInfo()), file.path(logdir, "sessionInfo.txt"))

# ---------------------------------------------------------------------
# 11. Console summary
# ---------------------------------------------------------------------
cat("\n==== PAIRED HSC vs MBE DNA CONCENTRATION (MBE vs HSC) ====\n")
print(stat_df, row.names = FALSE, digits = 6)
cat("\nSignificance label (paired test):", sig_label, "\n")
cat("P value:", fmtp(p_val), "\n")
cat("r_rb   :", fmtp(r_rb), " (MBE vs HSC)\n")
cat("\nOutputs written to:", outdir, "\n")
cat("  Figure: HSC_MBE_concentration_MBE_vs_HSC.{tiff,pdf,svg}\n")
cat("  Logs  : direction_change_log.txt, analysis_log.txt, pairing_validation.txt, sessionInfo.txt\n")
cat("  R     : HSC_MBE_concentration_paired_MBE_vs_HSC.R\n")
cat("  Tables: HSC_MBE_concentration_MBE_vs_HSC_statistics.csv\n")
