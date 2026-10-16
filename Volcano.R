suppressPackageStartupMessages({
  library(readxl)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(svglite)
  library(grid)
})

# ---------- Directories (all new files are written only under Fig4d) ----------
BASE_OUT   <- "D:/Desktop/Article1/Figures/Revised/Fig4d"
FIGURE_DIR <- file.path(BASE_OUT, "Figure")
LOGS_DIR   <- file.path(BASE_OUT, "Logs")
R_DIR      <- file.path(BASE_OUT, "R")
TABLES_DIR <- file.path(BASE_OUT, "Tables")
for (d in c(BASE_OUT, FIGURE_DIR, LOGS_DIR, R_DIR, TABLES_DIR)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# ---------- Input files (read-only) ----------
PATH_PAIRED <- "D:/Desktop/Article1/Figures/Revised/Report-Volcano Chart Testing Method/Tables/Genus_Differential_Abundance_Paired_Wilcoxon_FDR.xlsx"
PATH_TOP5   <- "D:/Desktop/Article1/Figures/Revised/Full differential genus for volcano/Tables/Full_Differential_Genera_Paired_Wilcoxon_q0.10_absLog2FC.xlsx"

# ---------- Analysis / plotting parameters ----------
Q_THRESHOLD      <- 0.10   # Significance threshold (q < 0.10)
LOG2FC_THRESHOLD <- 1      # |log2FC| > 1
MAIN_OUT_NAME    <- "Fig4d_Volcano_Final"

FIG_TIFF        <- file.path(FIGURE_DIR, paste0(MAIN_OUT_NAME, ".tiff"))
FIG_PDF         <- file.path(FIGURE_DIR, paste0(MAIN_OUT_NAME, ".pdf"))
FIG_SVG         <- file.path(FIGURE_DIR, paste0(MAIN_OUT_NAME, ".svg"))
FIG_PNG_PREVIEW <- file.path(BASE_OUT, "_preview.png")   # For visual inspection only; not a deliverable

# Point display parameters (Issue 2: size 1.8->1.5, alpha 0.8->0.70)
POINT_SIZE  <- 1.5
POINT_ALPHA <- 0.70

# Colors (unchanged)
COLOR_MAP <- c("Up" = "#E63946", "Down" = "#52B788", "Not significant" = "#D3D2CF")

# Logging utilities
LOG_FILE  <- file.path(LOGS_DIR, "analysis_log.txt")
WARN_FILE <- file.path(LOGS_DIR, "warnings_and_errors.txt")
ts_now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
logf <- function(...) {
  line <- paste0("[", ts_now(), "] ", paste0(..., collapse = ""))
  cat(line, "\n", file = LOG_FILE, append = TRUE); cat(line, "\n")
}
warnf <- function(...) {
  line <- paste0("[", ts_now(), "] WARN: ", paste0(..., collapse = ""))
  cat(line, "\n", file = WARN_FILE, append = TRUE)
  cat(line, "\n", file = LOG_FILE, append = TRUE); cat(line, "\n")
}

logf("=== Fig4d_Volcano_Final.R started (Fig4d--3 -> Fig4d) ===")

# 1) Read data
df <- readxl::read_excel(PATH_PAIRED, sheet = "Differential_results")
t5 <- readxl::read_excel(PATH_TOP5, sheet = 1)

# 2) Map fields and classify significance (fully consistent with the current rules)
main_df <- df %>% mutate(
  Regulation = case_when(
    q_value < Q_THRESHOLD & log2FC >  LOG2FC_THRESHOLD  ~ "Up",
    q_value < Q_THRESHOLD & log2FC < -LOG2FC_THRESHOLD  ~ "Down",
    TRUE ~ "Not significant"
  )
)
main_df$Regulation <- factor(main_df$Regulation,
                             levels = c("Up", "Down", "Not significant"))

n_up <- sum(main_df$Regulation == "Up")
n_dn <- sum(main_df$Regulation == "Down")
n_ns <- sum(main_df$Regulation == "Not significant")
logf(sprintf("Classification counts: Up=%d, Down=%d, Not significant=%d (total %d)",
             n_up, n_dn, n_ns, nrow(main_df)))

# 3) Top 5 labels (from the specified file; select the first five after sorting by Rank)
if ("Rank" %in% names(t5)) t5 <- t5 %>% arrange(Rank)
highlight_genera <- t5$Genus[1:5]
label_df <- main_df %>% filter(Genus %in% highlight_genera)
label_df <- label_df[match(highlight_genera, label_df$Genus), ]
logf("Top 5 genera: ", paste(highlight_genera, collapse = ", "))
if (nrow(label_df) < 5) warnf("label_df fewer than five rows matched; some top-five genera were not matched")

# 4) Colors (unchanged) + dynamic legend labels
color_map <- COLOR_MAP
legend_labels <- c(
  "Up"              = paste0("Higher in MBE: ", n_up),
  "Down"            = paste0("Higher in HSC: ", n_dn),
  "Not significant" = paste0("Not significant: ", n_ns)
)

# Expected-count validation
expected <- c(Up = 57, Down = 143, `Not significant` = 2677)
counts_match <- all(c(n_up, n_dn, n_ns) == unname(expected))
if (!counts_match)
  warnf(sprintf("Dynamic counts (%d/%d/%d) do not match the expected counts (57/143/2677)", n_up, n_dn, n_ns))

# 5) Top5 Explicit label positions and leader-line endpoints (this revision: rearrange labels and shorten the Y-axis)
#    The five genera have identical y values (-log10 q = 1.065465), so their data points lie on the same horizontal line.
#    Automatic label repulsion alone would cause overlap. Following the hand-drawn template, distribute the labels with the following three above: Bordetella/
#    Penicillium/Brucella(within the reserved space from 1.5 to 1.8), and the following two below: Neisseria/
#    Salmonella(below the threshold line); the text centers are not mechanically aligned with the data-point x coordinates.
#    The leader lines are short diagonal arrows in different directions, with a clear gap between each endpoint and the text center
#    to prevent arrowheads from touching or crossing text or intersecting other labels or arrows.
#    The actual data-point coordinates remain unchanged.
XRANGE <- c(-3.6, 13.0)      # X-axis range unchanged
YRANGE <- c(0, 1.8)          # Revision: reduce the panel upper boundary from 1.9 to 1.8

label_pos <- data.frame(
  Genus    = c("Bordetella", "Neisseria", "Penicillium", "Salmonella", "Brucella"),
  label_x  = c( 5.30,  6.30,  8.95, 10.35, 11.75),   # Text-center x coordinate (spread layout)
  label_y  = c( 1.66,  0.50,  1.72,  0.32,  1.66),   # Text-center y coordinate (staggered vertically)
  end_x    = c( 5.75,  6.62,  8.95,  9.90, 11.48),   # Leader-line endpoint x coordinate (diagonal arrow)
  end_y    = c( 1.54,  0.66,  1.58,  0.50,  1.54),   # Leader-line endpoint y coordinate (gap from text)
  stringsAsFactors = FALSE
)
label_df <- left_join(label_df, label_pos, by = "Genus")
logf("Top5 Explicit label positions (x, y): ",
     paste(sprintf("%s(%.2f,%.2f)", label_df$Genus, label_df$label_x, label_df$label_y),
           collapse = ", "))

# 6) Plot
volcano_plot <- ggplot(main_df, aes(x = log2FC, y = neg_log10_q, color = Regulation)) +
  geom_point(size = POINT_SIZE, alpha = POINT_ALPHA) +
  scale_color_manual(values = color_map, labels = legend_labels) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", linewidth = 0.6) +
  geom_hline(yintercept = -log10(Q_THRESHOLD), linetype = "dashed", color = "black", linewidth = 0.6) +
  # Arrows: actual data points -> leader-line endpoints (end_x, end_y); text centers are at (label_x, label_y).
  #       Leave a gap between each leader-line endpoint and text center to prevent arrowheads from touching or crossing text.
  geom_segment(
    data = label_df,
    aes(x = log2FC, y = neg_log10_q, xend = end_x, yend = end_y),
    inherit.aes = FALSE,
    color = "black",
    linewidth = 0.4,
    arrow = arrow(length = unit(0.2, "cm"), type = "closed")
  ) +
  geom_text(
    data = label_df,
    aes(x = label_x, y = label_y, label = Genus),
    inherit.aes = FALSE,
    color = "black",
    size = 3.5
  ) +
  labs(x = "log2(Fold Change) [MBE vs HSC]", y = "-log10(q-value)") +
  # The highest major Y-axis tick is 1.5 (the panel upper boundary of 1.8 is controlled by coord_cartesian).
  scale_y_continuous(breaks = seq(0, 1.5, by = 0.5)) +
  theme_minimal(base_size = 14) +
  theme(
    legend.title   = element_blank(),
    legend.position = "right",
    axis.text      = element_text(size = 13, color = "black"),
    axis.title     = element_text(size = 15, color = "black")
  ) +
  coord_cartesian(xlim = XRANGE, ylim = YRANGE)

# 7) Export three formats (all from the same final plot object, volcano_plot)
ggsave(FIG_TIFF, volcano_plot, width = 7, height = 5, dpi = 600,
       device = "tiff", compression = "lzw")
ggsave(FIG_PDF,  volcano_plot, width = 7, height = 5, device = cairo_pdf)
ggsave(FIG_SVG,  volcano_plot, width = 7, height = 5, device = svglite::svglite)
ggsave(FIG_PNG_PREVIEW, volcano_plot, width = 7, height = 5, dpi = 150, device = "png")

logf("Saved TIFF/PDF/SVG/PNG files: ", basename(FIG_TIFF), ", ", basename(FIG_PDF), ", ",
     basename(FIG_SVG), ", ", basename(FIG_PNG_PREVIEW))

# 8) Write plot_counts_validation.txt (data and count validation)
pcf <- c(
  "plot_counts_validation.txt",
  "=========================",
  paste0("Generation time: ", ts_now()),
  "",
  "Data integrity:",
  paste0("  Input rows: ", nrow(main_df)),
  paste0("  Number of log2FC NAs: ", sum(is.na(main_df$log2FC))),
  paste0("  Number of q_value NAs: ", sum(is.na(main_df$q_value))),
  "",
  "Classification counts (dynamic; consistent with the legend):",
  paste0("  Higher in MBE (Up)   : ", n_up),
  paste0("  Higher in HSC (Down) : ", n_dn),
  paste0("  Not significant      : ", n_ns),
  "",
  "Expected-count validation:",
  paste0("  Expected: Up=57, Down=143, Not significant=2677"),
  paste0("  Match: ", if (counts_match) "Yes" else "No (see warnings)"),
  "",
  "Coordinate ranges (this revision: Y-axis upper boundary 1.9 -> 1.8; highest major tick 1.5):",
  paste0("  xlim = ", paste(XRANGE, collapse = ", ")),
  paste0("  ylim = ", paste(YRANGE, collapse = ", ")),
  paste0("  Y-axis major ticks: ", paste(seq(0, 1.5, by = 0.5), collapse = ", "))
)
writeLines(pcf, file.path(LOGS_DIR, "plot_counts_validation.txt"))
logf("Wrote plot_counts_validation.txt")

# 9) Write top5_label_validation.txt (validate label positions, leader-line gaps, and spread layout)
label_df <- label_df %>% mutate(point_x = log2FC, point_y = neg_log10_q)
in_range <- all(label_df$label_x >= XRANGE[1] & label_df$label_x <= XRANGE[2] &
                label_df$label_y >= YRANGE[1] & label_df$label_y <= YRANGE[2] &
                label_df$end_x   >= XRANGE[1] & label_df$end_x   <= XRANGE[2] &
                label_df$end_y   >= YRANGE[1] & label_df$end_y   <= YRANGE[2])
gap_ok <- all(sqrt((label_df$end_x - label_df$label_x)^2 +
                   (label_df$end_y - label_df$label_y)^2) > 0.08)
top_genera <- c("Brucella", "Penicillium", "Bordetella")
bot_genera <- c("Salmonella", "Neisseria")
alt_ok <- all(label_df$label_y[label_df$Genus %in% top_genera] >
              label_df$point_y[label_df$Genus %in% top_genera]) &&
          all(label_df$label_y[label_df$Genus %in% bot_genera] <
              label_df$point_y[label_df$Genus %in% bot_genera])
hline_clear <- all(label_df$label_y[label_df$Genus %in% bot_genera] + 0.05 <
                   -log10(Q_THRESHOLD))
dir_ok <- all(
  (label_df$end_y[label_df$Genus %in% top_genera] > label_df$point_y[label_df$Genus %in% top_genera]) &
  (label_df$end_y[label_df$Genus %in% top_genera] < label_df$label_y[label_df$Genus %in% top_genera])
) && all(
  (label_df$end_y[label_df$Genus %in% bot_genera] < label_df$point_y[label_df$Genus %in% bot_genera]) &
  (label_df$end_y[label_df$Genus %in% bot_genera] > label_df$label_y[label_df$Genus %in% bot_genera])
)
dx <- outer(label_df$label_x, label_df$label_x, "-")
dy <- outer(label_df$label_y, label_df$label_y, "-")
d <- sqrt(dx^2 + dy^2); diag(d) <- Inf
min_dist <- min(d)
min_dx_samerow <- Inf
for (i in 1:5) for (j in 1:5) {
  if (i < j && abs(label_df$label_y[i] - label_df$label_y[j]) < 0.3)
    min_dx_samerow <- min(min_dx_samerow, abs(label_df$label_x[i] - label_df$label_x[j]))
}

tv <- c(
  paste0("top5_label_validation.txt (generated at ", ts_now(), ")"),
  "",
  "Top 5 text centers / leader-line endpoints and corresponding data points:",
  paste0("  ", sprintf("%-12s point (%.3f, %.4f) -> leader-line endpoint (%.2f, %.2f) -> text (%.2f, %.2f)",
                      label_df$Genus, label_df$point_x, label_df$point_y,
                      label_df$end_x, label_df$end_y,
                      label_df$label_x, label_df$label_y)),
  "",
  paste0("Text and leader-line endpoints are within the coordinate ranges: ", if (in_range) "Yes" else "No (!!)"),
  paste0("A clear gap (Euclidean distance > 0.08) remains between leader-line endpoints and text centers: ", if (gap_ok) "Yes" else "No (!!)"),
  paste0("Leader-line endpoints lie between the data points and text: ", if (dir_ok) "Yes" else "No (!!)"),
  paste0("Spread layout (top: Brucella/Penicillium/Bordetella; bottom: Salmonella/Neisseria): ",
         if (alt_ok) "Yes" else "No (!!)"),
  paste0("Lower labels do not overlap the y=1.0 threshold line: ", if (hline_clear) "Yes" else "No (!!)"),
  paste0("Minimum pairwise distance between label centers (data units): ", format(min_dist, digits = 4)),
  paste0("Minimum horizontal distance between labels at similar heights (data units): ", format(min_dx_samerow, digits = 4)),
  "",
  "Arrows: each is a short diagonal leader line starting at the actual data point; a gap is left between its endpoint and the text center, and the arrow does not cross the text.",
  "Layout: Bordetella/Penicillium/Brucella are above (y≈1.66-1.72), and Neisseria/Salmonella are below (y≈0.50/0.32).",
  "Note: Manually review the preview image _preview.png."
)
writeLines(tv, file.path(LOGS_DIR, "top5_label_validation.txt"))
logf("Wrote top5_label_validation.txt")

logf("=== Fig4d_Volcano_Final.R completed ===")
