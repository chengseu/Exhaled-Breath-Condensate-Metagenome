# Packages
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# Paths
file_path <- "D:/Desktop/Article/Codes/Data/CE-Concentration of PBA-MBs.xlsx"
out_base  <- "D:/Desktop/Article/Figures/CE-Concentration of PBA_MBs"

# Helper
to_num <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "%", "")
  x <- str_trim(x)
  suppressWarnings(as.numeric(x))
}

# Read (unchanged)
raw <- read_excel(file_path, sheet = 1) %>%
  rename(conc = 1, capture = 2, mean = 3, sd = 4) %>%
  mutate(conc = to_num(conc), capture = to_num(capture), mean = to_num(mean), sd = to_num(sd)) %>%
  fill(conc, .direction = "down")

# Points (unchanged)
points_df <- raw %>% filter(!is.na(capture)) %>% select(conc, capture)

# Summary (unchanged, with NA guard)
summary_df <- raw %>%
  group_by(conc) %>%
  summarise(
    mean = ifelse(all(is.na(mean)), mean(capture, na.rm = TRUE), first(na.omit(mean))),
    sd   = ifelse(all(is.na(sd)),   sd(capture,  na.rm = TRUE), first(na.omit(sd))),
    .groups = "drop"
  ) %>%
  arrange(conc) %>%
  filter(!is.na(mean) & !is.na(sd))

# Factor levels and labels (unchanged)
lvl_str <- sprintf("%.1f", c(0.5, 1.0, 1.5, 2.0))
summary_df$conc <- factor(sprintf("%.1f", summary_df$conc), levels = lvl_str)
points_df$conc  <- factor(sprintf("%.1f", points_df$conc),  levels = lvl_str)

# Palette (unchanged)
pal <- c("0.5" = "#CAF0F8", "1.0" = "#FFC2D1", "1.5" = "#CAF0F8", "2.0" = "#CAF0F8")

# Plot (only the requested tweak: match error bar thickness to axis ticks)
p <- ggplot(summary_df, aes(x = conc, y = mean, fill = conc)) +
  geom_col(width = 0.45, color = NA) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.15, size = 0.3, color = "black") +
  geom_point(data = points_df,
             aes(x = conc, y = capture),
             position = position_jitter(width = 0.08, height = 0),
             size = 2, shape = 16, color = "grey50", alpha = 0.85) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(
    breaks = seq(80, 100, 5),
    minor_breaks = seq(80, 100, 2.5)
  ) +
  labs(x = "Concentration of PBA-MBs  (mg/mL)", y = "Capture efficiency (%)") +
  coord_cartesian(ylim = c(80, 100)) +
  theme_classic(base_size = 14) +
  theme(
    axis.line   = element_blank(),                         # replace default left/bottom axes with a unified outer border
    axis.ticks  = element_line(color = "black", linewidth = 0.6),  # tick linewidth 0.6
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)  # thicker outer border
  )

# Save (unchanged)
ggsave(paste0(out_base, ".tiff"), plot = p, width = 4.5, height = 4, dpi = 300, compression = "lzw")
ggsave(paste0(out_base, ".pdf"),  plot = p, width = 4.5, height = 4)
ggsave(paste0(out_base, ".svg"),  plot = p, width = 4.5, height = 4)

# Show
p
