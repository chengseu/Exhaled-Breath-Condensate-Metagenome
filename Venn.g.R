# ============================================================
# 0. Load required R packages
# ============================================================
pkgs <- c("ggplot2", "cowplot", "readxl", "dplyr", "stringr")
for (pkg in pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  library(pkg, character.only = TRUE)
}

# ============================================================
# 1. Read genus-level data and count (corrected)
# ============================================================
file_path <- "D:/Desktop/Article/Codes/Data/Venn.g.xlsx"
ext <- tools::file_ext(file_path)

if (ext %in% c("xlsx", "xls")) {
  df_raw <- readxl::read_excel(file_path)
} else if (ext == "csv") {
  df_raw <- read.csv(file_path, stringsAsFactors = FALSE)
} else if (ext == "txt") {
  df_raw <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
} else {
  stop("不支持的文件格式：", ext)
}

# ---- 1.1 Normalize column names ----
colnames(df_raw)[1:3] <- c("Genus", "HSC", "MBE")

# ---- 1.2 Extract genus from taxonomy ----
df_raw$Genus <- stringr::str_extract(df_raw$Genus, "g__[^;]+")
df_raw$Genus <- stringr::str_replace(df_raw$Genus, "g__", "")
df_raw$Genus <- trimws(df_raw$Genus)

# remove empty genus
df_raw <- df_raw[!is.na(df_raw$Genus) & df_raw$Genus != "", ]

# ---- 1.3 Ensure HSC / MBE are numeric ----
df_raw$HSC <- suppressWarnings(as.numeric(df_raw$HSC))
df_raw$MBE <- suppressWarnings(as.numeric(df_raw$MBE))

# ---- 1.4 Remove genera that are 0 or NA in both groups ----
df_use <- df_raw %>%
  filter(
    (!is.na(HSC) & HSC > 0) |
      (!is.na(MBE) & MBE > 0)
  )

# ---- 1.5 Build unique genus sets ----
genera_HSC <- df_use %>%
  filter(!is.na(HSC) & HSC > 0) %>%
  pull(Genus) %>%
  unique()

genera_MBE <- df_use %>%
  filter(!is.na(MBE) & MBE > 0) %>%
  pull(Genus) %>%
  unique()

# ---- 1.6 Count the three categories ----
only_HSC <- setdiff(genera_HSC, genera_MBE)
only_MBE <- setdiff(genera_MBE, genera_HSC)
shared    <- intersect(genera_HSC, genera_MBE)

only_HSC_count <- length(only_HSC)
only_MBE_count <- length(only_MBE)
shared_count   <- length(shared)

cat("HSC 特有属数量：", only_HSC_count, "\n")
cat("MBE 特有属数量：", only_MBE_count, "\n")
cat("共有属数量：", shared_count, "\n")

# ============================================================
# 2. Build ellipse data (perfect proportions)
# ============================================================
ellipse_points <- function(cx, cy, rx = 1.3, ry = 0.9, n = 400) {
  theta <- seq(0, 2 * pi, length.out = n)
  data.frame(x = cx + rx * cos(theta), y = cy + ry * sin(theta))
}

cx1 <- -0.55
cx2 <-  0.55
cy  <-  0

ellipse_HSC <- ellipse_points(cx1, cy); ellipse_HSC$group <- "HSC"
ellipse_MBE <- ellipse_points(cx2, cy); ellipse_MBE$group <- "MBE"
ellipses_df <- rbind(ellipse_HSC, ellipse_MBE)

# ============================================================
# 3. Colors (soft blend)
# ============================================================
col_HSC <- "#CAF0F8"  # light blue
col_MBE <- "#D8F3DC"  # light green

x_HSC_text    <- cx1 - 0.65
x_shared_text <- 0
x_MBE_text    <- cx2 + 0.65
y_text        <- 0

# ============================================================
# 4. Main plot: two ellipses + numeric labels
# ============================================================
p_venn <- ggplot() +
  geom_polygon(data = subset(ellipses_df, group == "HSC"),
               aes(x, y), fill = col_HSC, alpha = 0.7) +
  geom_polygon(data = subset(ellipses_df, group == "MBE"),
               aes(x, y), fill = col_MBE, alpha = 0.7) +
  annotate("text", x = x_HSC_text,    y = y_text, label = only_HSC_count, size = 5) +
  annotate("text", x = x_shared_text, y = y_text, label = shared_count,   size = 5) +
  annotate("text", x = x_MBE_text,    y = y_text, label = only_MBE_count, size = 5) +
  coord_fixed(xlim = c(-2.5, 2.5), ylim = c(-1.6, 1.6), expand = FALSE) +
  theme_void() +
  theme(plot.margin = margin(t = 10, r = 10, b = 5, l = 10))

# ============================================================
# 5. Legend
# ============================================================
legend_plot <- ggplot() +
  annotate("point", x = 0.46, y = 0.8, shape = 21, size = 6,
           fill = col_HSC, color = "#CAF0F8" ) +
  annotate("text", x = 0.49, y = 0.8, label = "HSC", size = 5, hjust = 0) +
  annotate("point", x = 0.46, y = 0.5, shape = 21, size = 6,
           fill = col_MBE, color = "#D8F3DC" ) +
  annotate("text", x = 0.49, y = 0.5, label = "MBE", size = 5, hjust = 0) +
  xlim(0, 1) + ylim(0, 1) +
  theme_void() +
  theme(plot.margin = margin(t = -50, r = 0, b = 0, l = 0))

# ============================================================
# 6. Combine main plot + legend
# ============================================================
combined_plot <- cowplot::plot_grid(
  p_venn,
  legend_plot,
  ncol = 1,
  rel_heights = c(0.9, 0.1)
)

# ============================================================
# 7. Save figure
# ============================================================
outdir <- "D:/Desktop/Article/Figures"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

ggsave(file.path(outdir, "Venn.g.tiff"),
       combined_plot, width = 6, height = 4, dpi = 300, device = "tiff")

ggsave(file.path(outdir, "Venn.g.pdf"),
       combined_plot, width = 6, height = 4, device = cairo_pdf)

ggsave(file.path(outdir, "Venn.g.svg"),
       combined_plot, width = 6, height =4, device = "svg")

cat("✅ VThe Venn diagram has been saved toD:/Desktop/Article/Figures\n")
