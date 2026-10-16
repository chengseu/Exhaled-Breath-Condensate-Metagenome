# Load required packages
library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)

# -----------------------------
# File paths (all outputs written to the Revised/Fig6 directory only)
# -----------------------------
data_file <- "D:/Desktop/Article1/Codes/Result/KEGG/LEfSe_LDA.annotated.xlsx"

out_tiff <- "D:/Desktop/Article1/Figures/Revised/Fig6/Figure/LEfSe_LDA_threshold_revised.tiff"
out_pdf  <- "D:/Desktop/Article1/Figures/Revised/Fig6/Figure/LEfSe_LDA_threshold_revised.pdf"
out_svg  <- "D:/Desktop/Article1/Figures/Revised/Fig6/Figure/LEfSe_LDA_threshold_revised.svg"

dir.create(dirname(out_tiff), recursive = TRUE, showWarnings = FALSE)

# -----------------------------
# Read data (identical to original)
# -----------------------------
df <- read_excel(data_file, col_names = TRUE) %>%
  mutate(
    Module      = as.character(Module),
    Group       = as.character(Group),
    LDA         = suppressWarnings(as.numeric(LDA)),
    Description = as.character(Description)
  ) %>%
  filter(!is.na(Module), !is.na(Group), !is.na(LDA))

# -----------------------------
# Extract the first annotation before "||" as the plot label (identical to original)
# -----------------------------
df <- df %>%
  mutate(Description_short = str_split_fixed(Description, "\\|\\|", 2)[,1])

# -----------------------------
# Define left/right display rules (identical to original)
# HSC -> negative LDA
# MBE -> keep positive
# -----------------------------
df <- df %>%
  mutate(LDA_signed = ifelse(Group == "HSC", -abs(LDA), abs(LDA)))

# -----------------------------
# Filter |LDA| >= 3 (actual threshold, confirmed from source)
# -----------------------------
plot_df <- df %>%
  filter(abs(LDA_signed) >= 3) %>%
  arrange(LDA_signed) %>%
  mutate(Module_factor = factor(Module, levels = Module[order(LDA_signed)]))

# Actual threshold (absolute value) used for the dashed threshold lines
LDA_THRESHOLD <- 3

# -----------------------------
# Color scheme (identical to original, unchanged)
# -----------------------------
bar_colors <- c("MBE" = "#87BBA2",  # green
                "HSC" = "#48CAE4")  # blue

# -----------------------------
# Draw bar plot (identical structure; only additions are the threshold lines,
# the Arial font, and a modest font-size increase for print readability)
# -----------------------------
p <- ggplot(plot_df, aes(x = LDA_signed, y = Module_factor, fill = Group)) +
  geom_col(width = 0.7, alpha = 0.95) +
  scale_fill_manual(values = bar_colors, name = "Group") +
  labs(x = "LDA SCORE (log 10)", y = NULL) +
  scale_x_continuous(breaks = seq(-4, 4, 1), limits = c(-4.5, 4.5)) +
  # --- existing thin grey coordinate reference lines (PRESERVED) ---
  geom_vline(xintercept = seq(-4, 4, 1), linetype = "dashed",
             color = "grey80", linewidth = 0.4) +
  # --- NEW: explicit LDA threshold lines (|LDA| = 3), black & thicker ---
  geom_vline(xintercept = c(-LDA_THRESHOLD, LDA_THRESHOLD),
             linetype = "dashed", color = "black", linewidth = 0.8) +
  theme_minimal(base_size = 15, base_family = "Arial") +
  theme(
    panel.border = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 16),
    legend.position = "none",
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10)
  ) +
  # show short annotation next to the bars (identical positions; modest size increase)
  geom_text(data = plot_df %>% filter(Group == "MBE"),
            aes(x = -0.1, y = Module_factor, label = Description_short),
            hjust = 1, size = 4.6, color = "black", family = "Arial") +
  geom_text(data = plot_df %>% filter(Group == "HSC"),
            aes(x =  0.1, y = Module_factor, label = Description_short),
            hjust = 0, size = 4.6, color = "black", family = "Arial")

# -----------------------------
# Top labels: two columns at top-right (no border) (identical, modest size increase)
# -----------------------------
y_top <- length(levels(plot_df$Module_factor))

p <- p +
  # HSC label (top-right, left column)
  annotate("rect", xmin = 1.5, xmax = 2.0, ymin = y_top + 1, ymax = y_top + 1.5,
           fill = "#48CAE4", color = NA) +
  annotate("text", x = 2.1, y = y_top + 1.25, label = "HSC",
           hjust = 0, size = 4.5, fontface = "bold", family = "Arial") +
  # MBE label (top-right, right column)
  annotate("rect", xmin = 3.5, xmax = 4.0, ymin = y_top + 1, ymax = y_top + 1.5,
           fill = "#87BBA2", color = NA) +
  annotate("text", x = 4.1, y = y_top + 1.25, label = "MBE",
           hjust = 0, size = 4.5, fontface = "bold", family = "Arial") +
  coord_cartesian(clip = "off")

# -----------------------------
# Save as TIFF (300 dpi), PDF, SVG (true vector)
# Use ragg / cairo_pdf / svglite so that "Arial" resolves consistently
# across raster and vector outputs.
# -----------------------------
ragg::agg_tiff(out_tiff, width = 14, height = 12, units = "in", res = 300)
print(p)
dev.off()

cairo_pdf(out_pdf, width = 14, height = 12, family = "Arial")
print(p)
dev.off()

svglite::svglite(out_svg, width = 14, height = 12)
print(p)
dev.off()

cat("Revised figure saved to:\n", out_tiff, "\n", out_pdf, "\n", out_svg, "\n")
cat("LDA threshold used:", LDA_THRESHOLD, "(|LDA| >= 3; lines at -3 and +3)\n")
cat("Number of bars plotted:", nrow(plot_df), "\n")
