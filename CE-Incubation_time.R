# Packages
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# Paths
file_path <- "D:/Desktop/Article/Codes/Data/CE-Incubation time.xlsx"
out_base  <- "D:/Desktop/Article/Figures/CE-Incubation_time"

# Helper
to_num <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "%", "")
  x <- str_trim(x)
  suppressWarnings(as.numeric(x))
}

# Read
raw <- read_excel(file_path, sheet = 1) %>%
  rename(time = 1, capture = 2, mean = 3, sd = 4) %>%
  mutate(time = to_num(time), capture = to_num(capture), mean = to_num(mean), sd = to_num(sd)) %>%
  fill(time, .direction = "down")

# Points
points_df <- raw %>% filter(!is.na(capture)) %>% select(time, capture)

# Summary
summary_df <- raw %>%
  group_by(time) %>%
  summarise(
    mean = ifelse(all(is.na(mean)), mean(capture, na.rm = TRUE), first(na.omit(mean))),
    sd   = ifelse(all(is.na(sd)),   sd(capture,  na.rm = TRUE), first(na.omit(sd))),
    .groups = "drop"
  ) %>%
  arrange(time) %>%
  filter(!is.na(mean) & !is.na(sd))

# Factor levels and labels
lvl_str <- as.character(c(10, 20, 30, 40, 50, 60))
summary_df$time <- factor(summary_df$time, levels = c(10, 20, 30, 40, 50, 60))
points_df$time  <- factor(points_df$time,  levels = c(10, 20, 30, 40, 50, 60))

# Palette (consistent style)
pal <- c("10" = "#CAF0F8", "20" = "#CAF0F8", "30" = "#FFC2D1",
         "40" = "#FFC2D1", "50" = "#CAF0F8", "60" = "#CAF0F8")

# Plot
p <- ggplot(summary_df, aes(x = time, y = mean, fill = time)) +
  geom_col(width = 0.45, color = NA) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.15, size = 0.3, color = "black") +
  geom_point(data = points_df,
             aes(x = time, y = capture),
             position = position_jitter(width = 0.08, height = 0),
             size = 2, shape = 16, color = "grey50", alpha = 0.85) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(
    breaks = seq(80, 100, 5),
    minor_breaks = seq(80, 100, 2.5)
  ) +
  labs(x = "Incubation time (min)", y = "Capture efficiency (%)") +
  coord_cartesian(ylim = c(80, 100)) +
  theme_classic(base_size = 14) +
  theme(
    axis.line   = element_blank(),                         
    axis.ticks  = element_line(color = "black", linewidth = 0.6),  # tick linewidth 0.6
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(color = "black"),
    axis.title.x = element_text(margin = margin(t = 8)),
    axis.title.y = element_text(margin = margin(r = 8)),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8)  # thicker outer border
  )

# Save
ggsave(paste0(out_base, ".tiff"), plot = p, width = 4.5, height = 4, dpi = 300, compression = "lzw")
ggsave(paste0(out_base, ".pdf"),  plot = p, width = 4.5, height = 4)
ggsave(paste0(out_base, ".svg"),  plot = p, width = 4.5, height = 4)

# Show
p
