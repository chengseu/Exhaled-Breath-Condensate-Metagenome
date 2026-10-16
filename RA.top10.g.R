# ============================================================
# Load required R packages
# ============================================================
library(readxl)    # read Excel files
library(dplyr)     # data manipulation
library(tidyr)     # data tidying
library(stringr)   # string manipulation
library(ggplot2)   # plotting

# ============================================================
# 1. Set file and output paths
# ============================================================
file_path <- "D:/Desktop/Article/Codes/Data/RA.top10.g.xlsx"
outdir    <- "D:/Desktop/Article/Figures"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ============================================================
# 2. Read genus-level data
#    First column is the sample name; the remaining columns are relative abundances of each genus
# ============================================================
df_raw <- read_excel(file_path)
colnames(df_raw)[1] <- "Sample"

# Extract group info (HSC = high-speed centrifugation; MBE = magnetic bead enrichment)
df <- df_raw %>%
  mutate(Group = ifelse(str_detect(Sample, "^HSC"), "HSC", "MBE"))

# Convert to long format (each row = one sample x one genus)
df_long <- df %>%
  pivot_longer(cols = -c(Sample, Group),
               names_to = "Genus",
               values_to = "Abundance")

# ============================================================
# 3. Set genus order (legend order) and colors
#    Others must be placed last -> bottom of the legend
# ============================================================
genus_levels <- c(
  "Salmonella",
  "Alphainfluenzavirus",
  "Streptococcus",
  "Mycoplasma",
  "Listeria",
  "Faecalibacterium",
  "Prevotella",
  "Bacteroides",
  "Herbaspirillum",
  "Staphylococcus",
  "Others"
)

genus_colors <- c(
  "Salmonella"         = "#C9E4CA",
  "Alphainfluenzavirus"= "#ADE8F4",
  "Streptococcus"      = "#CCDBFD",
  "Mycoplasma"         = "#FBD1A2",
  "Listeria"           = "#FFE5EC",
  "Faecalibacterium"   = "#FFD6EF",
  "Prevotella"         = "#C5DEDD",
  "Bacteroides"        = "#E7CFEF",
  "Herbaspirillum"     = "#BCD4E6",
  "Staphylococcus"     = "#DFE7FA",
  "Others"             = "#F5EBE0"
)

df_long <- df_long %>%
  mutate(Genus = factor(Genus, levels = genus_levels))

# ============================================================
# 4. Compute group mean relative abundance (HSC vs MBE)
#    If you want one bar per sample, use Sample instead
# ============================================================
df_group <- df_long %>%
  group_by(Group, Genus) %>%
  summarise(Abundance = mean(Abundance, na.rm = TRUE), .groups = "drop") %>%
  mutate(Group = factor(Group, levels = c("HSC", "MBE")))

# ============================================================
# 5. Draw the genus-level relative abundance bar plot
# ============================================================
p <- ggplot(df_group, aes(x = Group, y = Abundance, fill = Genus)) +
  
  # stacked bar; reverse=TRUE -> Others always on top
  geom_bar(stat = "identity",
           position = position_stack(reverse = TRUE),
           width = 0.6) +
  
  # manual colors (exact HEX values as provided)
  scale_fill_manual(values = genus_colors, name = NULL) +
  
  # set y-axis breaks
  scale_y_continuous(
    breaks = c(0, 0.25, 0.50, 0.75, 1.00),
    labels = sprintf("%.2f", c(0, 0.25, 0.50, 0.75, 1.00)),
    limits = c(0, 1.00)
  ) +
  
  # axis titles (no x-axis title)
  labs(x = NULL, y = "Relative abundance") +
  
  # ============================================================
# Theme settings (solid outer border + dashed inner horizontal lines)
# ============================================================
theme_minimal(base_size = 14, base_family = "Arial") +
  theme(
    axis.text.x      = element_text(size = 12, color = "black"),
    axis.text.y      = element_text(size = 12, color = "black"),
    axis.title.y     = element_text(size = 14, color = "black"),
    
    # legend settings
    legend.position  = "right",
    legend.text      = element_text(size = 11, color = "black"),
    
    # outer border (solid)
    panel.border     = element_rect(color = "black", fill = NA, linewidth = 0.6),
    
    # inner horizontal dashed lines
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_line(color = "grey80",
                                      linetype = "dashed",
                                      linewidth = 0.4),
    
    plot.margin      = margin(t = 5, r = 5, b = 5, l = 5)
  )

# ============================================================
# 6. Save the figure (TIFF / PDF / SVG)
# ============================================================
ggsave(file.path(outdir, "RA.top10.g.tiff"),
       p, width = 4.5, height = 4.5, dpi = 300)

ggsave(file.path(outdir, "RA.top10.g.pdf"),
       p, width = 4.5, height = 4.5, device = cairo_pdf)

ggsave(file.path(outdir, "RA.top10.g.svg"),
       p, width = 4.5, height = 4.5)
