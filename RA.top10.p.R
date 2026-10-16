# ============================
# Load R packages (install first via install.packages() if missing)
# ============================
library(readxl)   # read Excel
library(dplyr)    # data manipulation
library(tidyr)    # data tidying
library(stringr)  # string manipulation
library(ggplot2)  # plotting

# ============================
# 1. Read data and basic preparation
# ============================
# File path (note: on Windows use / or \\)
file_path <- "D:/Desktop/Article/Codes/Data/RA.top10.p.xlsx"

# Read phylum-level relative abundance data
# Assumption: first column = sample name, the rest = relative abundance of each phylum (0-1)
df_raw <- read_excel(file_path)

# Rename the first column to Sample for downstream processing
df <- df_raw %>%
  rename(Sample = 1) %>%
  # identify groups from sample names: HSC (high-speed centrifugation) and MBE (magnetic bead enrichment)
  mutate(
    Group = ifelse(str_detect(Sample, "^HSC"), "HSC", "MBE")
  )

# Convert wide to long format: each row = one sample's abundance in one phylum
df_long <- df %>%
  pivot_longer(
    cols = -c(Sample, Group),
    names_to = "Phylum",
    values_to = "Abundance"
  )

# ============================
# 2. Unify phylum names to match those used in the color palette
#    (if column names have a p__ prefix, remove it here for matching)
# ============================
df_long <- df_long %>%
  mutate(
    Phylum_clean = str_replace(Phylum, "^p__", ""),  # remove leading p__
    Phylum_clean = ifelse(Phylum_clean == "Proteobacteria",
                          "Pseudomonadota", Phylum_clean)  # unify legacy names here if any
  )

# ============================
# 3. Set phylum order and colors (legend order and palette)
#    Note: Others is placed last in the factor -> shown at the bottom of the legend
# ============================
phylum_levels <- c(
  "Pseudomonadota",
  "Bacillota",
  "Negarnaviricota",
  "Mycoplasmatota",
  "Bacteroidota",
  "Actinomycetota",
  "Uroviricota",
  "Verrucomicrobiota",
  "Ascomycota",
  "Euryarchaeota",
  "Others"
)

phylum_colors <- c(
  "Pseudomonadota"   = "#C9E4CA",
  "Bacillota"        = "#ADE8F4",
  "Negarnaviricota"  = "#CCDBFD",
  "Mycoplasmatota"   = "#FBD1A2",
  "Bacteroidota"     = "#FFE5EC",
  "Actinomycetota"   = "#FFD6EF",
  "Uroviricota"      = "#C5DEDD",
  "Verrucomicrobiota"= "#E7CFEF",
  "Ascomycota"       = "#BCD4E6",
  "Euryarchaeota"    = "#DFE7FA",
  "Others"           = "#F5EBE0"
)

df_long <- df_long %>%
  mutate(
    Phylum_clean = factor(Phylum_clean, levels = phylum_levels)
  )

# ============================
# 4. To place Others at the top of the stacked bars:
#    use position_stack(reverse = TRUE) and put Others last in the factor
#    -> legend follows levels top-to-bottom; stacking is reversed so Others is on top
# ============================

# ============================
# 5. Compute overall group mean relative abundance (only two bars: HSC and MBE)
#    If you want one bar per sample, skip this step and map x to Sample
# ============================
df_group <- df_long %>%
  group_by(Group, Phylum_clean) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop")

# Ensure Group order: HSC on the left, MBE on the right
df_group <- df_group %>%
  mutate(Group = factor(Group, levels = c("HSC", "MBE")))

# ============================
# 6. Plot
# ============================
p <- ggplot(df_group, aes(x = Group, y = Abundance, fill = Phylum_clean)) +
  # stacked bar; reverse = TRUE keeps Others on top
  geom_bar(
    stat = "identity",
    position = position_stack(reverse = TRUE),
    width = 0.6
  ) +
  # manual colors
  scale_fill_manual(values = phylum_colors, name = NULL) +
  # set y-axis breaks and labels
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = sprintf("%.2f", c(0, 0.25, 0.50, 0.75, 1.00)),
    limits = c(0, 1.00)
  ) +
  # x-axis shows group names only, no title
  labs(
    x = NULL,                      # no x-axis title
    y = "Relative abundance"       # y-axis title
  ) +
  # theme: solid outer border + dashed inner horizontal lines
  theme_minimal(base_size = 14, base_family = "Arial") +
  theme(
    # axis text style
    axis.text.x      = element_text(size = 12, color = "black"),
    axis.text.y      = element_text(size = 12, color = "black"),
    axis.title.y     = element_text(size = 14, color = "black"),
    axis.title.x     = element_blank(),  # ensure no x-axis title again
    # legend style (right)
    legend.position  = "right",
    legend.text      = element_text(size = 11, color = "black"),
    # panel border: solid outer border
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
    # major grid: keep horizontal dashed lines only
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey80", linetype = "dashed", linewidth = 0.4),
    # remove extra margins
    plot.margin      = margin(t = 5, r = 5, b = 5, l = 5)
  )

# ============================
# 7. Save the figure as TIFF / PDF / SVG
# ============================
outdir <- "D:/Desktop/Article/Figures"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

tiff_file <- file.path(outdir, "RA.top10.p.tiff")
pdf_file  <- file.path(outdir, "RA.top10.p.pdf")
svg_file  <- file.path(outdir, "RA.top10.p.svg")

# save as high-resolution TIFF (commonly used for journals)
ggsave(tiff_file, p, width = 4.5, height = 4.5, dpi = 300)
# save as PDF vector
ggsave(pdf_file,  p, width = 4.5, height = 4.5, device = cairo_pdf)
# save as SVG vector
ggsave(svg_file,  p, width = 4.5, height = 4.5)
