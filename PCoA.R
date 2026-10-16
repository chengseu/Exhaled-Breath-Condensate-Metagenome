# ============================
# Load required R packages
# ============================
library(readxl)       # read Excel files
library(ggplot2)      # plotting
library(dplyr)        # data manipulation
library(stringr)      # string manipulation
library(ggpubr)       # figure layout and composition

# ============================
# Set input/output paths
# ============================
pcoa_path <- "D:/Desktop/Article/Codes/Data/PCoA.xlsx"
outdir    <- "D:/Desktop/Article/Figures"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ============================
# Read PCoA coordinate data
# ============================
pcoa_data <- read_excel(pcoa_path)
colnames(pcoa_data)[1:3] <- c("Sample", "PCoA1", "PCoA2")   # rename columns to PCoA1 and PCoA2

# ============================
# Extract group information
# ============================
pcoa_data <- pcoa_data %>%
  mutate(
    Group  = ifelse(str_detect(Sample, "^HSC"), "HSC", "MBE"),
    Gender = ifelse(str_detect(Sample, "M"), "Male", "Female"),
    ID     = str_extract(Sample, "\\d+")
  )

# ============================
# Set color scheme
# ============================
fill_colors  <- c("HSC" = "#ADE8F4", "MBE" = "#C9E4CA")   # ellipse fill colors (low saturation)
point_colors <- c("HSC" = "#48CAE4", "MBE" = "#87BBA2")   # sample point colors (high saturation)

# ============================
# Draw the PCoA scatter plot
# ============================
pcoa_plot <- ggplot(pcoa_data, aes(x = PCoA1, y = PCoA2, color = Group, fill = Group)) +
  geom_point(size = 3, alpha = 0.9) +  # sample points
  stat_ellipse(type = "norm", level = 0.95, geom = "polygon", alpha = 0.4, color = NA) +  # 95% confidence ellipse (filled)
  stat_ellipse(type = "norm", level = 0.95, geom = "path", linetype = "dashed", color = "grey70", size = 0.6) +  # 95% confidence ellipse (dashed boundary)
  scale_color_manual(values = point_colors) +
  scale_fill_manual(values = fill_colors) +
  labs(x = "PCoA1 (60.91%)", y = "PCoA2 (25.88%)") +  # axis variance contributions
  annotate("text",
           x = -0.65,
           y = -0.5,
           label = "Anosim:\nR = 0.345, P = 0.003",
           size = 4.2, hjust = 0, color = "black", fontface = "plain") +  # two-line ANOSIM result at lower-left
  theme_minimal(base_size = 14, base_family = "Arial") +
  theme(
    legend.position   = "right",
    axis.title.x      = element_text(size = 14),
    axis.title.y      = element_text(size = 14),
    axis.text         = element_text(size = 11, color = "black"),
    axis.ticks.length = unit(0.2, "cm"),
    axis.ticks        = element_line(color = "black"),
    panel.border      = element_rect(color = "black", fill = NA, linewidth = 0.75),  # outer border linewidth increased
    panel.grid.major  = element_line(color = "grey80", linetype = "dashed", linewidth = 0.3),
    panel.grid.minor  = element_blank(),
    plot.margin       = margin(t = 10, r = 10, b = 10, l = 10)
  )

# ============================
# Save the figure in three formats
# ============================
ggsave(file.path(outdir, "PCoA.tiff"), plot = pcoa_plot, width = 6, height = 5, dpi = 300)
ggsave(file.path(outdir, "PCoA.pdf"), plot = pcoa_plot, width = 6, height = 5, device = cairo_pdf)
ggsave(file.path(outdir, "PCoA.svg"), plot = pcoa_plot, width = 6, height = 5, device = cairo_ps)
