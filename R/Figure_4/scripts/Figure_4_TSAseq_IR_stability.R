# =========================================================
# Figure 4: TSA-seq overlap by IR stability class
#
# Expected repository structure:
# Figure_4/
# ├── data/
# │   ├── merged_data_GC_filtered.csv
# │   └── GSE81553_SON_TSA-Seq_Decile_Color_Condition2.bed
# ├── results/
# └── scripts/
#     └── Figure_4_TSAseq_IR_stability.R
#
# This script assumes it is run from the scripts/ folder.
#
# Input file notes:
# merged_data_GC_filtered.csv is a prefiltered table in which introns were
# retained only if they had:
# - at least 5 non-missing nuclear PIR values across the time course
# - Nuc_UT_TPM_mean >= 5
# =========================================================

# -----------------------------
# Libraries
# -----------------------------
library(GenomicRanges)
library(dplyr)
library(ggplot2)
library(tidyr)
library(RColorBrewer)
library(readr)

# -----------------------------
# Paths
# -----------------------------
bed_file <- "../data/GSE81553_SON_TSA-Seq_Decile_Color_Condition2.bed"
input_file <- "../data/merged_data_GC_filtered.csv"
output_dir <- "../results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

if (!file.exists(bed_file)) {
  stop("BED file not found: ", bed_file)
}

if (!file.exists(input_file)) {
  stop("Input file not found: ", input_file)
}

# -----------------------------
# Load data
# -----------------------------
bed_data <- read.table(
  bed_file,
  header = FALSE
)

merged_data_GC_filtered <- read_csv(
  input_file,
  show_col_types = FALSE
)

# -----------------------------
# Convert relevant columns to numeric
# -----------------------------
merged_data_GC_filtered <- merged_data_GC_filtered %>%
  mutate(
    Nuc_UT_mean = as.numeric(as.character(Nuc_UT_mean)),
    Nuc_30_mean = as.numeric(as.character(Nuc_30_mean)),
    Nuc_2h_mean = as.numeric(as.character(Nuc_2h_mean)),
    Nuc_4h_mean = as.numeric(as.character(Nuc_4h_mean)),
    Cyt_UT_mean = as.numeric(as.character(Cyt_UT_mean)),
    Cyt_30_mean = as.numeric(as.character(Cyt_30_mean)),
    Cyt_2h_mean = as.numeric(as.character(Cyt_2h_mean)),
    Cyt_4h_mean = as.numeric(as.character(Cyt_4h_mean)),
    Nuc_UT_TPM_mean = as.numeric(as.character(Nuc_UT_TPM_mean)),
    Nuc_30_TPM_mean = as.numeric(as.character(Nuc_30_TPM_mean)),
    Nuc_2h_TPM_mean = as.numeric(as.character(Nuc_2h_TPM_mean)),
    Nuc_4h_TPM_mean = as.numeric(as.character(Nuc_4h_TPM_mean)),
    Cyt_UT_TPM_mean = as.numeric(as.character(Cyt_UT_TPM_mean)),
    Cyt_30_TPM_mean = as.numeric(as.character(Cyt_30_TPM_mean)),
    Cyt_2h_TPM_mean = as.numeric(as.character(Cyt_2h_TPM_mean)),
    Cyt_4h_TPM_mean = as.numeric(as.character(Cyt_4h_TPM_mean))
  )

cat("Loaded rows:", nrow(merged_data_GC_filtered), "\n")
cat("Duplicated EVENT entries:", sum(duplicated(merged_data_GC_filtered$EVENT)), "\n")

# -----------------------------
# Recreate stability columns
# -----------------------------
merged_data_GC_filtered <- merged_data_GC_filtered %>%
  mutate(
    Nuc_A_UT_IR_stability = Nuc_A_UT_1.x / Nuc_UT_mean,
    Nuc_B_UT_IR_stability = Nuc_B_UT_1.x / Nuc_UT_mean,
    Nuc_A_30_IR_stability = Nuc_A_30_1.x / Nuc_UT_mean,
    Nuc_B_30_IR_stability = Nuc_B_30_1.x / Nuc_UT_mean,
    Nuc_A_2h_IR_stability = Nuc_A_2h_1.x / Nuc_UT_mean,
    Nuc_B_2h_IR_stability = Nuc_B_2h_1.x / Nuc_UT_mean,
    Nuc_A_4h_IR_stability = Nuc_A_4h_1.x / Nuc_UT_mean,
    Nuc_B_4h_IR_stability = Nuc_B_4h_1.x / Nuc_UT_mean
  ) %>%
  mutate(
    Nuc_UT_TPM_stability = Nuc_UT_TPM_mean / Nuc_UT_TPM_mean,
    Nuc_30_TPM_stability = Nuc_30_TPM_mean / Nuc_UT_TPM_mean,
    Nuc_2h_TPM_stability = Nuc_2h_TPM_mean / Nuc_UT_TPM_mean,
    Nuc_4h_TPM_stability = Nuc_4h_TPM_mean / Nuc_UT_TPM_mean
  ) %>%
  mutate(
    Cyt_UT_TPM_stability = Cyt_UT_TPM_mean / Cyt_UT_TPM_mean,
    Cyt_30_TPM_stability = Cyt_30_TPM_mean / Cyt_UT_TPM_mean,
    Cyt_2h_TPM_stability = Cyt_2h_TPM_mean / Cyt_UT_TPM_mean,
    Cyt_4h_TPM_stability = Cyt_4h_TPM_mean / Cyt_UT_TPM_mean
  ) %>%
  mutate(
    Nuc_UT_IR_stability = Nuc_UT_mean / Nuc_UT_mean,
    Nuc_30_IR_stability = Nuc_30_mean / Nuc_UT_mean,
    Nuc_2h_IR_stability = Nuc_2h_mean / Nuc_UT_mean,
    Nuc_4h_IR_stability = Nuc_4h_mean / Nuc_UT_mean
  ) %>%
  mutate(
    Cyt_UT_IR_stability = Cyt_UT_mean / Cyt_UT_mean,
    Cyt_30_IR_stability = Cyt_30_mean / Cyt_UT_mean,
    Cyt_2h_IR_stability = Cyt_2h_mean / Cyt_UT_mean,
    Cyt_4h_IR_stability = Cyt_4h_mean / Cyt_UT_mean
  )

# -----------------------------
# One row per gene when needed
# -----------------------------
merged_data_GC_filtered_GENE <- merged_data_GC_filtered %>%
  group_by(GENE) %>%
  top_n(1, Nuc_UT_mean) %>%
  ungroup()

# -----------------------------
# Assign stability classes
# -----------------------------
IR_UT <- merged_data_GC_filtered %>%
  filter(
    Nuc_UT_mean >= 30 &
      (
        (Nuc_A_UT_1.x >= 25 & Nuc_B_UT_1.x >= 25) |
          (is.na(Nuc_A_UT_1.x) & Nuc_B_UT_1.x >= 25) |
          (Nuc_A_UT_1.x >= 25 & is.na(Nuc_B_UT_1.x))
      )
  ) %>%
  mutate(stability = "IR")

stableR4h <- merged_data_GC_filtered %>%
  filter(
    Nuc_UT_mean >= 30 &
      Nuc_4h_IR_stability >= 0.5 &
      Nuc_30_IR_stability >= 0.5 &
      Nuc_2h_IR_stability >= 0.5 &
      (
        (Nuc_A_UT_1.x >= 25 & Nuc_B_UT_1.x >= 25) |
          (is.na(Nuc_A_UT_1.x) & Nuc_B_UT_1.x >= 25) |
          (Nuc_A_UT_1.x >= 25 & is.na(Nuc_B_UT_1.x))
      ) &
      Nuc_30_mean >= 25 &
      Nuc_2h_mean >= 25 &
      Nuc_4h_mean >= 25
  ) %>%
  mutate(stability = "stable 4h")

unstableR30min <- merged_data_GC_filtered %>%
  filter(
    Nuc_UT_mean >= 30 &
      Nuc_30_IR_stability < 0.5 &
      Nuc_2h_IR_stability < 0.5 &
      Nuc_4h_IR_stability < 0.5 &
      (
        (Nuc_A_UT_1.x >= 25 & Nuc_B_UT_1.x >= 25) |
          (is.na(Nuc_A_UT_1.x) & Nuc_B_UT_1.x >= 25) |
          (Nuc_A_UT_1.x >= 25 & is.na(Nuc_B_UT_1.x))
      )
  ) %>%
  mutate(stability = "unstable30min")

unstableR2h <- merged_data_GC_filtered %>%
  filter(
    Nuc_UT_mean >= 30 &
      Nuc_30_IR_stability >= 0.5 &
      Nuc_2h_IR_stability < 0.5 &
      Nuc_4h_IR_stability < 0.5 &
      Nuc_30_mean >= 25 &
      (
        (Nuc_A_UT_1.x >= 25 & Nuc_B_UT_1.x >= 25) |
          (is.na(Nuc_A_UT_1.x) & Nuc_B_UT_1.x >= 25) |
          (Nuc_A_UT_1.x >= 25 & is.na(Nuc_B_UT_1.x))
      )
  ) %>%
  mutate(stability = "unstable2h")

unstableR4h <- merged_data_GC_filtered %>%
  filter(
    Nuc_UT_mean >= 30 &
      Nuc_4h_IR_stability < 0.5 &
      Nuc_30_IR_stability >= 0.5 &
      Nuc_2h_IR_stability >= 0.5 &
      Nuc_30_mean >= 25 &
      Nuc_2h_mean >= 25 &
      (
        (Nuc_A_UT_1.x >= 25 & Nuc_B_UT_1.x >= 25) |
          (is.na(Nuc_A_UT_1.x) & Nuc_B_UT_1.x >= 25) |
          (Nuc_A_UT_1.x >= 25 & is.na(Nuc_B_UT_1.x))
      )
  ) %>%
  mutate(stability = "unstable 4h")

noIR <- merged_data_GC_filtered %>%
  filter(Nuc_UT_mean < 30) %>%
  mutate(stability = "no IR")

noIR_GENE <- merged_data_GC_filtered_GENE %>%
  filter(Nuc_UT_mean < 30) %>%
  mutate(stability = "no_IR_GENE") %>%
  distinct(GENE, .keep_all = TRUE)

cat("Number of noIR_GENE rows:", nrow(noIR_GENE), "\n")

# -----------------------------
# Build filtered combined table
# -----------------------------
combined_data <- bind_rows(
  stableR4h,
  unstableR4h,
  unstableR30min,
  unstableR2h,
  noIR
)

combined_dataTPMno100 <- combined_data %>%
  filter(
    is.na(Cyt_A_UT_1.x) | is.na(Cyt_B_UT_1.x) |
      (Cyt_A_UT_1.x < 99 & Cyt_B_UT_1.x < 99)
  )

cat("Rows in combined_dataTPMno100:", nrow(combined_dataTPMno100), "\n")
cat("Duplicated EVENTs in combined_dataTPMno100:", sum(duplicated(combined_dataTPMno100$EVENT)), "\n")

# -----------------------------
# Remove overlapping genes between IR classes
# -----------------------------
sets_list <- list(
  "Stable 4h" = unique(stableR4h$GENE),
  "Unstable 4h" = unique(unstableR4h$GENE),
  "Unstable 30min" = unique(unstableR30min$GENE),
  "Unstable 2h" = unique(unstableR2h$GENE)
)

exclusive_unstable30min <- setdiff(
  sets_list[["Unstable 30min"]],
  union(union(sets_list[["Stable 4h"]], sets_list[["Unstable 4h"]]), sets_list[["Unstable 2h"]])
)

exclusive_unstable2h <- setdiff(
  sets_list[["Unstable 2h"]],
  union(sets_list[["Stable 4h"]], sets_list[["Unstable 4h"]])
)

exclusive_unstable4h <- setdiff(
  unstableR4h$GENE,
  stableR4h$GENE
)

df_stable4h <- stableR4h %>%
  dplyr::select(GENE) %>%
  mutate(GENE_stability = "stable_4h")

df_unstable4h <- data.frame(
  GENE = exclusive_unstable4h,
  GENE_stability = "unstable_4h"
)

df_unstable30min <- data.frame(
  GENE = exclusive_unstable30min,
  GENE_stability = "unstable_30min"
)

df_unstable2h <- data.frame(
  GENE = exclusive_unstable2h,
  GENE_stability = "unstable_2h"
)

gene_stability_df <- bind_rows(
  df_stable4h,
  df_unstable4h,
  df_unstable30min,
  df_unstable2h
)

gene_stability_df_unique <- gene_stability_df %>%
  distinct(GENE, .keep_all = TRUE)

combined_exclusive <- combined_dataTPMno100 %>%
  filter(GENE %in% gene_stability_df_unique$GENE) %>%
  left_join(gene_stability_df_unique, by = "GENE") %>%
  filter(stability != "no IR")

combined_exclusive_unique <- combined_exclusive %>%
  distinct(GENE, .keep_all = TRUE)

# -----------------------------
# Build TSA-Seq GRanges
# -----------------------------
bed_data <- bed_data[grepl("Group_\\d+$", bed_data$V4), ]

bed_granges <- GRanges(
  seqnames = bed_data$V1,
  ranges = IRanges(start = bed_data$V2, end = bed_data$V3)
)
mcols(bed_granges)$group <- bed_data$V4

# -----------------------------
# Split by stability class
# -----------------------------
unstable_30IR <- combined_exclusive_unique[combined_exclusive_unique$GENE_stability == "unstable_30min", ]
unstable_2hIR <- combined_exclusive_unique[combined_exclusive_unique$GENE_stability == "unstable_2h", ]
unstable_4hIR <- combined_exclusive_unique[combined_exclusive_unique$GENE_stability == "unstable_4h", ]
stable_4hr    <- combined_exclusive_unique[combined_exclusive_unique$GENE_stability == "stable_4h", ]
noIR_df_input <- noIR_GENE[noIR_GENE$stability == "no_IR_GENE", ]

# -----------------------------
# Keep one representative intron per gene
# -----------------------------
extract_representative_intron <- function(df) {
  df <- df[!is.na(df$COORD), ]
  df <- df[!duplicated(df$GENE), ]
  df
}

unstable_30_df <- extract_representative_intron(unstable_30IR)
unstable_2h_df <- extract_representative_intron(unstable_2hIR)
unstable_4h_df <- extract_representative_intron(unstable_4hIR)
stable_df      <- extract_representative_intron(stable_4hr)
noIR_df        <- extract_representative_intron(noIR_df_input)

# -----------------------------
# Convert coordinates to GRanges
# -----------------------------
coord_to_granges <- function(df) {
  coord_parts <- do.call(rbind, strsplit(df$COORD, ":|-"))
  GRanges(
    seqnames = coord_parts[, 1],
    ranges = IRanges(
      start = as.integer(coord_parts[, 2]),
      end = as.integer(coord_parts[, 3])
    ),
    gene = df$GENE
  )
}

results_list_granges <- list(
  "unstable 30min" = coord_to_granges(unstable_30_df),
  "unstable 2h"    = coord_to_granges(unstable_2h_df),
  "unstable 4h"    = coord_to_granges(unstable_4h_df),
  "stable"         = coord_to_granges(stable_df),
  "noIR"           = coord_to_granges(noIR_df)
)

# -----------------------------
# Overlap percentages by TSA group
# -----------------------------
get_overlap_percentages_overlap_only <- function(gr, label) {
  hits <- findOverlaps(gr, bed_granges)
  
  df <- data.frame(
    intron_id = queryHits(hits),
    group = bed_granges$group[subjectHits(hits)],
    stringsAsFactors = FALSE
  )
  
  n_overlapping_introns <- n_distinct(df$intron_id)
  
  summary <- df %>%
    group_by(group) %>%
    summarise(count = n_distinct(intron_id), .groups = "drop") %>%
    mutate(
      percent = 100 * count / n_overlapping_introns,
      stability = label,
      overlap_positive_n = n_overlapping_introns
    )
  
  return(summary)
}

tsa_overlap_plot <- bind_rows(
  get_overlap_percentages_overlap_only(results_list_granges[["unstable 30min"]], "unstable 30min"),
  get_overlap_percentages_overlap_only(results_list_granges[["unstable 2h"]], "unstable 2h"),
  get_overlap_percentages_overlap_only(results_list_granges[["unstable 4h"]], "unstable 4h"),
  get_overlap_percentages_overlap_only(results_list_granges[["stable"]], "stable"),
  get_overlap_percentages_overlap_only(results_list_granges[["noIR"]], "noIR")
)

tsa_overlap_plot$group <- factor(tsa_overlap_plot$group, levels = paste0("Group_", 1:10))
tsa_overlap_plot$stability <- factor(
  tsa_overlap_plot$stability,
  levels = c("stable", "unstable 30min", "unstable 2h", "unstable 4h", "noIR")
)

tsa_overlap_plot %>%
  group_by(stability) %>%
  summarise(total_percent = sum(percent), .groups = "drop") %>%
  print()

# -----------------------------
# Plot
# -----------------------------
t <- ggplot(tsa_overlap_plot, aes(x = stability, y = percent, fill = group)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(9, "Blues"))(10)) +
  labs(
    title = "TSA-Seq Decile Distribution by IR Stability Group",
    x = "IR Stability Group",
    y = "Overlap Percentage (%)",
    fill = "TSA-Seq Decile"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(t)
output_dir <- "../results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(
  filename = file.path(output_dir, "TSASeq_IR_stability_stacked_barplot.pdf"),
  plot = t,
  width = 7,
  height = 5
)

# -----------------------------
# Export percentage table
# -----------------------------
output_dir <- "../results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

percent_table <- tsa_overlap_plot %>%
  dplyr::select(group, stability, percent) %>%
  pivot_wider(names_from = stability, values_from = percent)

print(percent_table, n = Inf)

write.csv(
  percent_table,
  file.path(output_dir, "TSAseq_stability_group_percentages_overlap_positive_corrected.csv"),
  row.names = FALSE
)

# -----------------------------
# Proximal overlap summary
# -----------------------------
get_overlap_positive_proximal <- function(gr, label) {
  hits <- findOverlaps(gr, bed_granges)
  
  df <- data.frame(
    intron_id = queryHits(hits),
    group = bed_granges$group[subjectHits(hits)],
    stringsAsFactors = FALSE
  ) %>%
    group_by(intron_id) %>%
    summarise(
      proximal = any(group %in% c("Group_9", "Group_10")),
      .groups = "drop"
    ) %>%
    mutate(stability = label)
  
  return(df)
}

tsa_proximal_overlap_only <- bind_rows(
  get_overlap_positive_proximal(results_list_granges[["stable"]], "stable"),
  get_overlap_positive_proximal(results_list_granges[["unstable 30min"]], "unstable 30min"),
  get_overlap_positive_proximal(results_list_granges[["unstable 2h"]], "unstable 2h"),
  get_overlap_positive_proximal(results_list_granges[["unstable 4h"]], "unstable 4h"),
  get_overlap_positive_proximal(results_list_granges[["noIR"]], "noIR")
)

proximal_overlap_positive_summary <- tsa_proximal_overlap_only %>%
  group_by(stability) %>%
  summarise(
    overlap_positive_n = n(),
    proximal_n = sum(proximal),
    non_proximal_n = sum(!proximal),
    percent_proximal = 100 * mean(proximal),
    .groups = "drop"
  )

print(proximal_overlap_positive_summary)

# -----------------------------
# Fisher tests
# -----------------------------
test_proximal_enrichment_overlap_only <- function(group1, group2, df) {
  df_sub <- df %>% filter(stability %in% c(group1, group2))
  
  tbl <- table(df_sub$stability, df_sub$proximal)
  test <- fisher.test(tbl)
  
  pct1 <- mean(df_sub$proximal[df_sub$stability == group1]) * 100
  pct2 <- mean(df_sub$proximal[df_sub$stability == group2]) * 100
  
  data.frame(
    comparison = paste(group1, "vs", group2),
    value = sprintf("Groups 9-10 overlap: %.1f%% vs %.1f%%", pct1, pct2),
    group = paste(group1, "|", group2),
    test = "Fisher exact",
    statistic = unname(test$estimate),
    `p-value` = test$p.value,
    check.names = FALSE
  )
}

fisher_results <- bind_rows(
  test_proximal_enrichment_overlap_only("stable", "unstable 30min", tsa_proximal_overlap_only),
  test_proximal_enrichment_overlap_only("stable", "unstable 2h", tsa_proximal_overlap_only),
  test_proximal_enrichment_overlap_only("stable", "unstable 4h", tsa_proximal_overlap_only),
  test_proximal_enrichment_overlap_only("stable", "noIR", tsa_proximal_overlap_only)
)

fisher_results$p.adj <- p.adjust(fisher_results$`p-value`, method = "BH")

print(fisher_results)

output_dir <- "../results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

write.csv(
  fisher_results,
  file.path(output_dir, "TSAseq_Group9_10_enrichment_vs_stable_corrected.csv"),
  row.names = FALSE
)
