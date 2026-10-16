# Packages
library(readxl)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)

# Paths
file_path <- "D:/Desktop/Article/Codes/Data/CE-Storage time.xlsx"
out_base  <- "D:/Desktop/Article/Figures/CE-Storage_time"

# Helper
to_num <- function(x) {
  x <- as.character(x)
  x <- str_replace_all(x, "%", "")
  x <- str_trim(x)
  suppressWarnings(as.numeric(x))
}

# Read
raw <- read_excel(file_path, sheet = 1) %>%
  rename(storage = 1, capture = 2, mean = 3, sd = 4) %>%
  mutate(storage = to_num(storage), capture = to_num(capture), mean = to_num(mean), sd = to_num(sd)) %>%
  fill(storage, .direction = "down")

# Points
points_df <- raw %>% filter(!is.na(capture)) %>% select(storage, capture)

# Summary
summary_df <- raw %>%
  group_by(storage) %>%
  summarise(
    mean = ifelse(all(is.na(mean)), mean(capture, na.rm = TRUE), first(na.omit(mean))),
    sd   = ifelse(all(is.na(sd)),   sd(capture,  na.rm = TRUE), first(na.omit(sd))),
    .groups = "drop"
  ) %>%
  arrange(storage) %>%
  filter(!is.na(mean) & !is.na(sd))

# Factor levels and labels
lvl_str <- as.character(c(1, 2, 3, 4, 5))
summary_df$storage <- factor(summary_df$storage, levels = c(1, 2, 3, 4, 5))
points_df$storage  <- factor(points_df$storage,  levels = c(1, 2, 3, 4, 5))

# Palette (consistent style)
pal <- c("1" = "#CAF0F8", "2" = "#CAF0F8", "3" = "#CAF0F8",
         "4" = "#CAF0F8", "5" = "#CAF0F8")

# Plot
p <- ggplot(summary_df, aes(x = storage, y = mean, fill = storage)) +
  geom_col(width = 0.45, color = NA) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.15, size = 0.3, color = "black") +
  geom_point(data = points_df,
             aes(x = storage, y = capture),
             position = position_jitter(width = 0.08, height = 0),
             size = 2, shape = 16, color = "grey50", alpha = 0.85) +
  scale_fill_manual(values = pal, guide = "none") +
  scale_y_continuous(
    breaks = seq(50, 100, 10),          # major ticks: 50 to 100, step 10
    minor_breaks = seq(50, 100, 5)   # minor ticks: 50 to 100, step 5
  ) +
  labs(x = "Storage time (month)", y = "Capture efficiency (%)") +
  coord_cartesian(ylim = c(50, 100)) +  # display range changed to 50-100
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
