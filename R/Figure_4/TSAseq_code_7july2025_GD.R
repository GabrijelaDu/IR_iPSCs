
library(GenomicRanges)
library(dplyr)
library(ggplot2)

analysis_dir <- here::here("R", "Figure_4")
setwd(analysis_dir)

bed_path <- file.path(analysis_dir, "GSE81553_SON_TSA-Seq_Decile_Color_Condition2.bed")
if (!file.exists(bed_path)) {
  stop("Missing input file: ", bed_path, call. = FALSE)
}
# 4. Load TSA-Seq bed file
bed_data <- read.table(bed_path, header = FALSE)

bed_data <- bed_data[grepl("Group_\\d+$", bed_data$V4), ]
bed_granges <- GRanges(
  seqnames = bed_data$V1,
  ranges = IRanges(start = bed_data$V2, end = bed_data$V3)
)
mcols(bed_granges)$group <- bed_data$V4
unique(bed_data$V4)


# 1. Split merged_data by stability class
# I used combined_exclusive_unique data generated in stability groups script
# filtering out only 5 TPM nuc, not cytoplasmic, and no 100% in cyto
unstable_30IR <- combined_exclusive_unique[combined_exclusive_unique$GENE_stability == "unstable_30min", ]
unstable_2hIR <- combined_exclusive_unique[combined_exclusive_unique$GENE_stability == "unstable_2h", ]
unstable_4hIR <- combined_exclusive_unique[combined_exclusive_unique$GENE_stability == "unstable_4h", ]
stable_4hr    <- combined_exclusive_unique[combined_exclusive_unique$GENE_stability == "stable_4h", ]
noIR <- noIR_GENE[noIR_GENE$stability == "no_IR_GENE", ]

# 2. Function to extract one intron per gene
extract_representative_intron <- function(df) {
  df <- df[!is.na(df$COORD), ]
  df <- df[!duplicated(df$GENE), ]
  return(df)
}

unstable_30_df <- extract_representative_intron(unstable_30IR)
unstable_2h_df <- extract_representative_intron(unstable_2hIR)
unstable_4h_df <- extract_representative_intron(unstable_4hIR)
stable_df      <- extract_representative_intron(stable_4hr)
noIR_df      <- extract_representative_intron(noIR)

# 3. Convert to GRanges
coord_to_granges <- function(df) {
  coord_parts <- do.call(rbind, strsplit(df$COORD, ":|-"))
  GRanges(
    seqnames = coord_parts[, 1],
    ranges = IRanges(start = as.integer(coord_parts[, 2]), end = as.integer(coord_parts[, 3])),
    gene = df$GENE
  )
}

results_list_granges <- list(
  "unstable 30min" = coord_to_granges(unstable_30_df),
  "unstable 2h"    = coord_to_granges(unstable_2h_df),
  "unstable 4h"    = coord_to_granges(unstable_4h_df),
  "stable"         = coord_to_granges(stable_df),
  "noIR"         = coord_to_granges(noIR_df)
)


# 5. Compute percent overlaps
get_overlap_percentages <- function(gr, label) {
  hits <- findOverlaps(gr, bed_granges)
  groups <- bed_granges[subjectHits(hits)]$group
  query_hits <- queryHits(hits)
  
  # Only include introns that overlapped any Group_1–10
  n_overlapping_introns <- length(unique(query_hits))
  
  summary <- as.data.frame(table(groups))
  summary$percent <- (summary$Freq / n_overlapping_introns) * 100
  summary$stability <- label
  colnames(summary) <- c("group", "count", "percent", "stability")
  return(summary)
}

unstable_30_tsa_df <- get_overlap_percentages(results_list_granges[["unstable 30min"]], "unstable 30min")
unstable_2h_tsa_df <- get_overlap_percentages(results_list_granges[["unstable 2h"]], "unstable 2h")
unstable_4h_tsa_df <- get_overlap_percentages(results_list_granges[["unstable 4h"]], "unstable 4h")
stable_tsa_df      <- get_overlap_percentages(results_list_granges[["stable"]], "stable")
noIR_tsa_df      <- get_overlap_percentages(results_list_granges[["noIR"]], "noIR")

# 6. Combine for plotting
tsa_overlap_all <- bind_rows(
  unstable_30_tsa_df,
  unstable_2h_tsa_df,
  unstable_4h_tsa_df,
  stable_tsa_df, noIR_tsa_df
)

tsa_overlap_all$group <- factor(tsa_overlap_all$group, levels = paste0("Group_", 1:10))
tsa_overlap_all$stability <- factor(tsa_overlap_all$stability,
                                    levels = c("stable", "unstable 30min", "unstable 2h", "unstable 4h", "noIR")
)

# 7. Plot
library(RColorBrewer)

t <-ggplot(tsa_overlap_all, aes(x = stability, y = percent, fill = group)) +
  geom_bar(stat = "identity", position = "stack", color = "black") +
  scale_fill_manual(values = colorRampPalette(brewer.pal(9, "Blues"))(10)) +
  labs(
    title = "TSA-Seq Decile Distribution by IR Stability Group (Gene Level)",
    x = "IR Stability Group",
    y = "Overlap Percentage (%)",
    fill = "TSA-Seq Decile"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("TSASeq_IR_stability_plot.pdf", plot = t, width = 7, height =5)


#### check if more proximal####
# Make a binary label: proximal (Group_9 or Group_10) vs. all others
# Function to test enrichment in Group_9/10 for any given group vs stable
test_proximal_enrichment <- function(group1, group2, df) {
  df_sub <- df[df$stability %in% c(group1, group2), ]
  tbl <- table(df_sub$stability == group1, df_sub$decile %in% c(9, 10))
  test <- fisher.test(tbl)
  
  data.frame(
    group1 = group1,
    group2 = group2,
    odds_ratio = test$estimate,
    conf_low = test$conf.int[1],
    conf_high = test$conf.int[2],
    p_value = test$p.value
  )
}

# Run comparisons: stable vs each other group
fisher_results <- bind_rows(
  test_proximal_enrichment("stable", "unstable 30min", tsa_gene_level),
  test_proximal_enrichment("stable", "unstable 2h", tsa_gene_level),
  test_proximal_enrichment("stable", "unstable 4h", tsa_gene_level),
  test_proximal_enrichment("stable", "noIR", tsa_gene_level)
)

# Save to CSV
write.csv(fisher_results, "TSAseq_Group9_10_enrichment_vs_stable.csv", row.names = FALSE)


#### Check percentages ####
############################
get_overlap_percentages_all <- function(gr, label) {
  total_introns <- length(gr)
  hits <- findOverlaps(gr, bed_granges)
  groups <- bed_granges[subjectHits(hits)]$group
  query_hits <- queryHits(hits)
  
  df <- data.frame(
    intron_id = query_hits,
    group = groups
  )
  
  summary <- df %>%
    group_by(group) %>%
    summarise(count = n_distinct(intron_id)) %>%
    ungroup()
  
  all_groups <- paste0("Group_", 1:10)
  summary <- summary %>%
    right_join(data.frame(group = all_groups), by = "group") %>%
    mutate(count = ifelse(is.na(count), 0, count)) %>%
    mutate(
      percent = (count / total_introns) * 100,
      stability = label
    )
  
  return(summary)
}

# Apply to all groups
tsa_overlap_all_full <- bind_rows(
  get_overlap_percentages_all(results_list_granges[["unstable 30min"]], "unstable 30min"),
  get_overlap_percentages_all(results_list_granges[["unstable 2h"]], "unstable 2h"),
  get_overlap_percentages_all(results_list_granges[["unstable 4h"]], "unstable 4h"),
  get_overlap_percentages_all(results_list_granges[["stable"]], "stable"),
  get_overlap_percentages_all(results_list_granges[["noIR"]], "noIR")
)

# Factor levels
tsa_overlap_all_full$group <- factor(tsa_overlap_all_full$group, levels = paste0("Group_", 1:10))
tsa_overlap_all_full$stability <- factor(tsa_overlap_all_full$stability,
                                         levels = c("stable", "unstable 30min", "unstable 2h", "unstable 4h", "noIR"))

# Optional: Check that % adds up to ~100% per group
tsa_overlap_all_full %>%
  group_by(stability) %>%
  summarise(total_percent = sum(percent))

### data frame groups and distances
get_intron_tsa_deciles <- function(gr, label) {
  hits <- findOverlaps(gr, bed_granges)
  df <- data.frame(
    gene = mcols(gr[queryHits(hits)])$gene,
    stability = label,
    group = bed_granges[subjectHits(hits)]$group
  )
  df$decile <- as.numeric(gsub("Group_", "", df$group))
  return(df)
}

tsa_gene_level <- bind_rows(
  get_intron_tsa_deciles(results_list_granges[["stable"]], "stable"),
  get_intron_tsa_deciles(results_list_granges[["unstable 30min"]], "unstable 30min"),
  get_intron_tsa_deciles(results_list_granges[["unstable 2h"]], "unstable 2h"),
  get_intron_tsa_deciles(results_list_granges[["unstable 4h"]], "unstable 4h"),
  get_intron_tsa_deciles(results_list_granges[["noIR"]], "noIR")
  
)
pairwise.wilcox.test(tsa_gene_level$decile, tsa_gene_level$stability, p.adjust.method = "BH")
# Create a wide-format table: rows = deciles, columns = stability groups, values = percentages
percent_table <- tsa_overlap_all_full %>%
  dplyr::select(stability, group, percent) %>%
  pivot_wider(names_from = stability, values_from = percent)

# View the table
print(percent_table, n = Inf)

# Save as CSV
write.csv(percent_table, "TSAseq_stability_group_percentages.csv", row.names = FALSE)

#######Check if an intron overlaps multple TSA seq regions:
# Run findOverlaps
hits <- findOverlaps(results_list_granges[["unstable 4h"]], bed_granges)

# Create a data frame of hits
df_hits <- data.frame(
  intron_id = queryHits(hits),
  group = mcols(bed_granges)[subjectHits(hits), "group"]
)

# Count how many times each intron appears (i.e. how many groups it overlaps)
overlap_counts <- df_hits %>%
  group_by(intron_id) %>%
  summarise(n_groups = n()) %>%
  arrange(desc(n_groups))

# See how many overlap more than one TSA group
table(overlap_counts$n_groups)

# Count how many unique decile groups each intron overlaps
query_group_map <- data.frame(
  query = queryHits(hits),
  group = mcols(bed_granges)$group[subjectHits(hits)]
)

# Count unique groups per intron
group_counts <- query_group_map %>%
  group_by(query) %>%
  summarise(n_groups = n_distinct(group)) %>%
  filter(n_groups > 1)

# Now retrieve the intron info
multi_overlap_introns <- unstable_30_df[group_counts$query, ]

# Show them
multi_overlap_introns


###################

# Revised function to assign TSA decile per intron
get_intron_tsa_deciles <- function(gr, label) {
  hits <- findOverlaps(gr, bed_granges)
  df <- data.frame(
    gene = mcols(gr[queryHits(hits)])$gene,
    stability = label,
    group = bed_granges[subjectHits(hits)]$group
  )
  df$decile <- as.numeric(gsub("Group_", "", df$group))
  return(df)
}

# Run findOverlaps
hits <- findOverlaps(results_list_granges[["stable"]], bed_granges)

# Create a data frame of hits
df_hits <- data.frame(
  intron_id = queryHits(hits),
  group = mcols(bed_granges)[subjectHits(hits), "group"]
)

# Count how many times each intron appears (i.e. how many groups it overlaps)
overlap_counts <- df_hits %>%
  group_by(intron_id) %>%
  summarise(n_groups = n()) %>%
  arrange(desc(n_groups))

# See how many overlap more than one TSA group
table(overlap_counts$n_groups)
# Apply to each IR group
unstable_30_gene_level <- get_intron_tsa_deciles(results_list_granges[["unstable 30min"]], "unstable 30min")
unstable_2h_gene_level <- get_intron_tsa_deciles(results_list_granges[["unstable 2h"]], "unstable 2h")
unstable_4h_gene_level <- get_intron_tsa_deciles(results_list_granges[["unstable 4h"]], "unstable 4h")
stable_gene_level      <- get_intron_tsa_deciles(results_list_granges[["stable"]], "stable")
noIR_gene_level      <- get_intron_tsa_deciles(results_list_granges[["noIR"]], "noIR")

# Combine into one long data frame
tsa_deciles_gene_level <- bind_rows(
  unstable_30_gene_level,
  unstable_2h_gene_level,
  unstable_4h_gene_level,
  stable_gene_level, noIR_gene_level
)


# Boxplot
ggplot(tsa_deciles_gene_level, aes(x = stability, y = decile, fill = stability)) +
  geom_boxplot(color = "black") +
  scale_y_continuous(breaks = 1:10) +
  labs(
    title = "TSA-Seq Deciles by IR Stability Group (Gene-level)",
    y = "TSA-Seq Decile (1 = far, 10 = near speckles)",
    x = "IR Stability"
  ) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_fill_brewer(palette = "Set2")

# Violin plot
ggplot(tsa_deciles_gene_level, aes(x = stability, y = decile, fill = stability)) +
  geom_violin(trim = FALSE) +
  geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
  scale_y_continuous(breaks = 1:10) +
  labs(
    title = "Distribution of TSA-Seq Deciles (Gene-level)",
    y = "TSA-Seq Decile",
    x = "IR Stability"
  ) +
  theme_minimal(base_size = 14) +
  scale_fill_brewer(palette = "Set2")

# CDF plot
ggplot(tsa_deciles_gene_level, aes(x = decile, color = stability)) +
  stat_ecdf(geom = "step", size = 1.2) +
  labs(
    title = "Cumulative Distribution of TSA-Seq Deciles",
    x = "TSA-Seq Decile (1 = far, 10 = near speckles)",
    y = "Cumulative Fraction"
  ) +
  theme_minimal(base_size = 14) +
  scale_color_brewer(palette = "Set2")

