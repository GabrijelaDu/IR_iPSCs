# =========================================================
# Figure 1: intron length, GC content, and nuclear enrichment
#
# Expected repository structure:
# Figure_1/
#   ├── data/
#   │   └── merged_data_GC_filtered.csv
# ├── results/
#   └── scripts/
#   └── Figure_1_length_gc_nuclear_enrichment.R
# This script assumes it is run from the scripts/ folder.

# Input file:
# merged_data_GC_filtered.csv is a filtered table in which introns were
# retained only if they had:
# - at least 5 non-missing nuclear PIR values across the time course
# - Nuc_UT_TPM_mean >= 5
#
# In this script, introns are classified as:
# - IR:   Nuc_UT_mean >= 30
# - noIR: Nuc_UT_mean < 30
#
# For nuclear enrichment analysis, one intron per gene is retained by selecting
# the intron with the highest Nuc_UT_mean.
# =========================================================

# -----------------------------
# Libraries
# -----------------------------
library(dplyr)
library(ggplot2)
library(readr)
library(tibble)
library(grid)

# -----------------------------
# Paths
# -----------------------------
input_file <- "../data/merged_data_GC_filtered.csv"
output_dir <- "../results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# -----------------------------
# Load data
# -----------------------------
merged_data_GC_filtered <- read_csv(
  input_file,
  show_col_types = FALSE
) %>%
  mutate(
    Cyt_UT_mean = as.numeric(Cyt_UT_mean),
    Cyt_30_mean = as.numeric(Cyt_30_mean),
    Cyt_2h_mean = as.numeric(Cyt_2h_mean),
    Cyt_4h_mean = as.numeric(Cyt_4h_mean),
    Nuc_UT_mean = as.numeric(Nuc_UT_mean),
    Nuc_30_mean = as.numeric(Nuc_30_mean),
    Nuc_2h_mean = as.numeric(Nuc_2h_mean),
    Nuc_4h_mean = as.numeric(Nuc_4h_mean),
    Nuc_UT_TPM_mean = as.numeric(Nuc_UT_TPM_mean),
    Cyt_UT_TPM_mean = as.numeric(Cyt_UT_TPM_mean),
    LENGTH = as.numeric(LENGTH),
    GC_Content = as.numeric(GC_Content)
  )

cat("Loaded rows:", nrow(merged_data_GC_filtered), "\n")
cat("Duplicated EVENT entries:", sum(duplicated(merged_data_GC_filtered$EVENT)), "\n")

# -----------------------------
# Shared plotting settings
# -----------------------------
group_cols <- c("noIR" = "darkgray", "IR" = "darkred")
group_cols_ne <- c("noIR (<30)" = "darkgray", "IR (>=30)" = "darkred")

theme_axes_out <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    axis.line = element_line(linewidth = 0.6, color = "black"),
    axis.ticks = element_line(color = "black"),
    axis.ticks.length = grid::unit(0.2, "cm"),
    axis.text = element_text(color = "black"),
    axis.title = element_text(color = "black"),
    legend.position = "none"
  )

# =========================================================
# PART 1: Intron-level length and GC-content analysis
# =========================================================

# IR is defined as Nuc_UT_mean >= 30
merged_data_GC_IR <- merged_data_GC_filtered %>%
  filter(Nuc_UT_mean >= 30) %>%
  mutate(Group = "IR")

merged_data_GC_noIR <- merged_data_GC_filtered %>%
  filter(Nuc_UT_mean < 30) %>%
  mutate(Group = "noIR")

combined_df <- bind_rows(merged_data_GC_IR, merged_data_GC_noIR) %>%
  mutate(Group = factor(Group, levels = c("noIR", "IR")))

# Summary tables for labels
length_means <- combined_df %>%
  group_by(Group) %>%
  summarise(
    mean_value = mean(LENGTH, na.rm = TRUE),
    .groups = "drop"
  )

gc_means <- combined_df %>%
  group_by(Group) %>%
  summarise(
    mean_value = mean(GC_Content, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# Intron length plot
# -----------------------------
p_len <- ggplot(combined_df, aes(x = Group, y = LENGTH, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  geom_text(
    data = length_means,
    aes(x = Group, y = mean_value, label = round(mean_value, 1)),
    inherit.aes = FALSE,
    vjust = -1.2,
    color = "black",
    size = 3.5
  ) +
  scale_fill_manual(values = group_cols) +
  scale_y_log10(labels = function(x) format(x, scientific = FALSE)) +
  labs(
    title = "Comparison of intron length by IR and noIR groups",
    x = "Group",
    y = "Intron length (log10)"
  ) +
  theme_axes_out

print(p_len)

# -----------------------------
# GC-content plot
# -----------------------------
p_gc <- ggplot(combined_df, aes(x = Group, y = GC_Content, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  geom_text(
    data = gc_means,
    aes(x = Group, y = mean_value, label = round(mean_value, 1)),
    inherit.aes = FALSE,
    vjust = -1.2,
    color = "black",
    size = 3.5
  ) +
  scale_fill_manual(values = group_cols) +
  labs(
    title = "Comparison of intron GC content by IR and noIR groups",
    x = "Group",
    y = "Intron GC (%)"
  ) +
  theme_axes_out

print(p_gc)

# -----------------------------
# Summary tables
# -----------------------------
mean_length_table <- combined_df %>%
  group_by(Group) %>%
  summarise(
    Mean_LENGTH = round(mean(LENGTH, na.rm = TRUE), 1),
    .groups = "drop"
  )

mean_GC_table <- combined_df %>%
  group_by(Group) %>%
  summarise(
    Mean_GC = round(mean(GC_Content, na.rm = TRUE), 1),
    .groups = "drop"
  )

print(mean_length_table)
print(mean_GC_table)

# -----------------------------
# Statistical tests
# -----------------------------
length_ttest <- t.test(log10(LENGTH) ~ Group, data = combined_df)
gc_ttest <- t.test(GC_Content ~ Group, data = combined_df)

length_p <- length_ttest$p.value
gc_p <- gc_ttest$p.value

cat(
  ifelse(
    length_p < .Machine$double.eps,
    "Welch's t-test p-value (log10 length): < 2.2e-16\n",
    sprintf("Welch's t-test p-value (log10 length): %.3e\n", length_p)
  )
)

cat(
  ifelse(
    gc_p < .Machine$double.eps,
    "Welch's t-test p-value (GC%%): < 2.2e-16\n",
    sprintf("Welch's t-test p-value (GC%%): %.3e\n", gc_p)
  )
)

# =========================================================
# PART 2: Gene-level nuclear enrichment analysis
# =========================================================

# Keep one row per gene by selecting the intron with the highest Nuc_UT_mean
unique_genes <- merged_data_GC_filtered %>%
  group_by(GENE) %>%
  slice_max(order_by = Nuc_UT_mean, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  filter(Cyt_UT_TPM_mean > 0) %>%
  mutate(
    log2FC = log2(Nuc_UT_TPM_mean / Cyt_UT_TPM_mean),
    Group = if_else(Nuc_UT_mean >= 30, "IR (>=30)", "noIR (<30)")
  )

cat("\nGene-level nuclear enrichment groups:\n")
print(table(unique_genes$Group, useNA = "ifany"))

unique_genes_plot <- unique_genes %>%
  filter(!is.na(Group), is.finite(log2FC)) %>%
  mutate(Group = factor(Group, levels = c("noIR (<30)", "IR (>=30)")))

nuc_means <- unique_genes_plot %>%
  group_by(Group) %>%
  summarise(
    mean_value = mean(log2FC, na.rm = TRUE),
    .groups = "drop"
  )

# -----------------------------
# Nuclear enrichment plot
# -----------------------------
p_nuc <- ggplot(unique_genes_plot, aes(x = Group, y = log2FC, fill = Group)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.shape = NA, fill = "white") +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 3, color = "black") +
  geom_text(
    data = nuc_means,
    aes(x = Group, y = mean_value, label = round(mean_value, 2)),
    inherit.aes = FALSE,
    vjust = -1.2,
    color = "black",
    size = 3.5
  ) +
  scale_fill_manual(values = group_cols_ne, drop = TRUE) +
  scale_x_discrete(drop = TRUE) +
  labs(
    title = "Nuclear enrichment",
    x = "Group",
    y = expression(log[2]~"(Nuc_UT_TPM_mean / Cyt_UT_TPM_mean)")
  ) +
  theme_axes_out

print(p_nuc)

# -----------------------------
# Nuclear enrichment test
# -----------------------------
unique_genes_plot_test <- unique_genes_plot %>%
  filter(!is.na(log2FC), is.finite(log2FC))

nuc_ttest <- t.test(log2FC ~ Group, data = unique_genes_plot_test)

cat(sprintf(
  "Welch's t-test (nuclear enrichment): t = %.8f, p = %.3e\n",
  unname(nuc_ttest$statistic), nuc_ttest$p.value
))
