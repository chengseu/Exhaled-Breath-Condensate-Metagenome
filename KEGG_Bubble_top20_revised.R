# Load required packages
library(readxl)
library(dplyr)
library(ggplot2)
library(stringr)

# -----------------------------
# Input file paths (original analysis data, read-only)
# -----------------------------
bubble_file <- "D:/Desktop/Article1/Codes/Result/KEGG/KEGG_Bubble_Data.xlsx"
score_file  <- "D:/Desktop/Article1/Codes/Result/KEGG/MBE-vs-HSC.psig.withScore.KEGG.xlsx"

# -----------------------------
# Output file paths (only written into Revised/Fig5)
# -----------------------------
out_dir <- "D:/Desktop/Article1/Figures/Revised/Fig5"
tiff_file <- file.path(out_dir, "Figure", "KEGG_Bubble_top20_revised.tiff")
pdf_file  <- file.path(out_dir, "Figure", "KEGG_Bubble_top20_revised.pdf")
svg_file  <- file.path(out_dir, "Figure", "KEGG_Bubble_top20_revised.svg")
log_file  <- file.path(out_dir, "Logs", "KEGG_Bubble_top20_revision_log.txt")

dir.create(dirname(tiff_file), recursive = TRUE, showWarnings = FALSE)
dir.create(dirname(log_file),  recursive = TRUE, showWarnings = FALSE)

log_lines <- character()
log <- function(...) { msg <- paste0(...); log_lines <<- c(log_lines, msg); cat(msg, "\n") }

# -----------------------------
# Read data (identical to original)
# -----------------------------
bubble_data <- read_excel(bubble_file) %>% as.data.frame()
score_data  <- read_excel(score_file) %>% as.data.frame()

# Ensure column names are consistent (identical to original)
colnames(bubble_data)[1] <- "Module"
colnames(score_data)[1]  <- "Module"
colnames(score_data)[6]  <- "Score"

# Select the top 20 modules by Score (identical to original)
top20 <- score_data %>%
  arrange(desc(Score)) %>%
  slice(1:20) %>%
  select(Module, Score)

# Merge bubble_data with Description (identical to original)
plot_data <- bubble_data %>%
  filter(Module %in% top20$Module) %>%
  inner_join(top20, by = "Module") %>%
  arrange(desc(Score))

# Set wrap threshold: wrap only when > 70 characters (replace || with newline) (identical)
plot_data$Description <- ifelse(
  nchar(plot_data$Description) > 70,
  str_replace_all(plot_data$Description, "\\|\\|", "\n"),
  plot_data$Description
)

# Set y-axis order: Score from high to low, highest at top (identical)
plot_data <- plot_data %>%
  mutate(Description = factor(Description, levels = rev(unique(Description))))

# Get min and max of adjPvalue (identical)
adj_min <- min(plot_data$adjPvalues, na.rm = TRUE)
adj_max <- max(plot_data$adjPvalues, na.rm = TRUE)

# -----------------------------
# Draw bubble plot (identical structure/color/size; only font settings revised)
# -----------------------------
p <- ggplot(plot_data, aes(x = logFC, y = Description,
                           size = mean_abundance,
                           color = adjPvalues)) +
  geom_point(alpha = 0.85) +
  scale_color_gradient(low = "#FF0000", high = "#0000FF",
                       limits = c(adj_min, adj_max),
                       name = "adjPvalue (FDR)") +
  scale_size(range = c(3, 12), name = "Mean Relative Abundance") +
  theme_minimal(base_size = 15, base_family = "Arial") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.6),
    panel.grid.major = element_line(color = "grey85", linewidth = 0.6),
    panel.grid.minor = element_line(color = "grey90", linewidth = 0.6),
    axis.text.y = element_text(size = 13, face = "bold", lineheight = 0.9, margin = margin(r = 8)),
    axis.text.x = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 12)),
    axis.title.x = element_text(margin = margin(t = 12)),
    legend.position = "right",
    plot.title = element_text(hjust = 0.5),
    plot.margin = margin(t = 12, r = 5, b = 12, l = 10)
  ) +
  labs(
    x = "log2 Fold Change (MBE vs HSC)",   # x-axis label (unchanged)
    y = "KEGG Module (Description)",
    size = "Mean Relative Abundance",
    color = "adjPvalue (FDR)",
    title = "KEGG Functional Module Bubble Plot (Top 20 by Score)"
  ) +
  coord_cartesian(clip = "off")

# -----------------------------
# Save as TIFF (300 dpi), PDF, SVG (true vector)
# ragg / cairo_pdf / svglite so "Arial" resolves consistently.
# -----------------------------
ragg::agg_tiff(tiff_file, width = 13, height = 12, units = "in", res = 300)
print(p)
dev.off()

cairo_pdf(pdf_file, width = 13, height = 12, family = "Arial")
print(p)
dev.off()

svglite::svglite(svg_file, width = 13, height = 12)
print(p)
dev.off()

# -----------------------------
# Write run log
# -----------------------------
log("=== KEGG Bubble Plot (Top 20) Revision Log ===")
log("Date: ", format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
log("R version: ", R.version.string)
log("ggplot2: ", as.character(packageVersion("ggplot2")))
log("readxl : ", as.character(packageVersion("readxl")))
log("dplyr  : ", as.character(packageVersion("dplyr")))
log("stringr: ", as.character(packageVersion("stringr")))
log("ragg   : ", as.character(packageVersion("ragg")))
log("svglite: ", as.character(packageVersion("svglite")))
log("")
log("Input files:")
log("  ", bubble_file)
log("  ", score_file)
log("Output files:")
log("  TIFF: ", tiff_file)
log("  PDF : ", pdf_file)
log("  SVG : ", svg_file)
log("")
log("Top 20 Module count: ", nrow(plot_data))
log("")
log("Font settings (revised):")
log("  base_family = Arial, base_size = 15")
log("  axis.text.y = 13 (bold, lineheight 0.9)")
log("  axis.text.x = 16")
log("  axis.title / legend / plot.title: inherit Arial base_size 15")
log("")
log("TIFF generated: ", file.exists(tiff_file))
log("PDF  generated: ", file.exists(pdf_file))
log("SVG  generated: ", file.exists(svg_file))
log("")
log("Note: only font family/size + output device changed. Data, Score,")
log("Top-20 selection, colors, bubble size, alpha, axis labels, title,")
log("legend, ordering and 13x12 inch layout are UNCHANGED.")

writeLines(log_lines, log_file)

cat("Revised bubble plot saved to:\n", tiff_file, "\n", pdf_file, "\n", svg_file, "\n")
cat("Log written to:\n", log_file, "\n")
cat("Number of modules plotted:", nrow(plot_data), "\n")
