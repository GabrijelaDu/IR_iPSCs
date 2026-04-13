# ───────────────────────────────────────────────────────────────────────────────
# SRSF7 boxplot from 2-column Excel file
# Columns in file: NT, INT3_DYRK
# Style matched to previous boxplot script
# Stats: two-sided Welch's t-test
# ───────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(readxl)

analysis_dir <- here::here("R", "Suppl_Fig_10B")
setwd(analysis_dir)

input_path <- file.path(analysis_dir, "2026_03_02_SRSF7_FISH_quantification_FINAL.xlsx")
if (!file.exists(input_path)) {
  stop("Missing input file: ", input_path, call. = FALSE)
}

# ── Load file ──────────────────────────────────────────────────────────────────
df <- read_excel(input_path)

# If needed, check column names
print(names(df))

# Keep only the two columns and make sure numeric
df <- df %>%
  select(NT, INT3_DYRK) %>%
  mutate(
    NT = suppressWarnings(as.numeric(NT)),
    INT3_DYRK = suppressWarnings(as.numeric(INT3_DYRK))
  )

# ── Long format ────────────────────────────────────────────────────────────────
df_long <- df %>%
  pivot_longer(
    cols = everything(),
    names_to = "Treatment",
    values_to = "value"
  ) %>%
  filter(is.finite(value)) %>%
  mutate(
    Treatment = factor(Treatment, levels = c("NT", "INT3_DYRK"))
  )

# ── Plot ───────────────────────────────────────────────────────────────────────
p <- ggplot(df_long, aes(x = Treatment, y = value)) +
  geom_boxplot(
    outlier.shape = NA,
    fill = "grey90",
    color = "black",
    width = 0.55,
    linewidth = 0.9,
    coef = 1.5
  ) +
  geom_jitter(
    width = 0.08,
    height = 0,
    alpha = 0.7,
    color = "#990066",
    size = 1.3
  ) +
  labs(
    title = "SRSF7",
    x = "Treatment",
    y = "Value"
  ) +
  scale_y_continuous(limits = c(0, NA)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks.y     = element_line(colour = "black"),
    axis.ticks.x     = element_line(colour = "black")
  )

print(p)

ggsave(
  "SRSF7_boxplot_NT_vs_INT3_DYRK.pdf",
  p,
  width = 4,
  height = 4
)

# ───────────────────────────────────────────────────────────────────────────────
# Summary table
# ───────────────────────────────────────────────────────────────────────────────
summary_table <- df_long %>%
  group_by(Treatment) %>%
  summarise(
    n = n(),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    .groups = "drop"
  )

print("SUMMARY TABLE")
print(summary_table)

write.csv(
  summary_table,
  "SRSF7_summary_table.csv",
  row.names = FALSE
)

# ───────────────────────────────────────────────────────────────────────────────
# Compact stats table: two-sided Welch's t-test
# ───────────────────────────────────────────────────────────────────────────────
x <- df_long %>% filter(Treatment == "NT") %>% pull(value)
y <- df_long %>% filter(Treatment == "INT3_DYRK") %>% pull(value)

tt <- t.test(
  x, y,
  paired = FALSE,
  var.equal = FALSE,
  alternative = "two.sided"
)

stats_compact <- tibble(
  comparison = "NT vs INT3_DYRK",
  value = "SRSF7",
  group = "all",
  test = "Welch's t-test",
  statistic = unname(tt$statistic),
  `p-value` = tt$p.value,
  p.adj = p.adjust(tt$p.value, method = "BH")
)

print("COMPACT STATS TABLE")
print(stats_compact)

write.csv(
  stats_compact,
  "SRSF7_stats_compact_Welch_ttest.csv",
  row.names = FALSE
)

# ───────────────────────────────────────────────────────────────────────────────
# Optional: paper-style reporting line
# ───────────────────────────────────────────────────────────────────────────────
cat(
  "\nWelch's two-sided t-test:\n",
  "t = ", round(unname(tt$statistic), 3),
  ", p = ", signif(tt$p.value, 3),
  "\n",
  sep = ""
)
