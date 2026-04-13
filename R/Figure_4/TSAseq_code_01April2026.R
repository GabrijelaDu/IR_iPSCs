library(GenomicRanges);library(dplyr);library(ggplot2)

analysis_dir <- here::here("R", "Figure_4")
setwd(analysis_dir)

bed_path <- file.path(analysis_dir, "GSE81553_SON_TSA-Seq_Decile_Color_Condition2.bed")
merged_data_path <- here::here("data", "merged_data_GC2.tsv")

if (!file.exists(bed_path)) {
  stop("Missing input file: ", bed_path, call. = FALSE)
}
if (!file.exists(merged_data_path)) {
  stop("Missing input file: ", merged_data_path, call. = FALSE)
}

# 4. Load TSA-Seq bed file
bed_data <- read.table(bed_path, header = FALSE)

merged_data_GC <- read.delim(merged_data_path, stringsAsFactors = FALSE)

merged_data_GC <- merged_data_GC %>% 
  mutate(
    Nuc_UT_mean = as.numeric(as.character(Nuc_UT_mean)),
    Nuc_30_mean = as.numeric(as.character(Nuc_30_mean)),
    Nuc_2h_mean = as.numeric(as.character(Nuc_2h_mean)),
    Nuc_4h_mean = as.numeric(as.character(Nuc_4h_mean)),
    Cyt_UT_mean = as.numeric(as.character(Cyt_UT_mean)),
    Cyt_30_mean = as.numeric(as.character(Cyt_30_mean)),
    Cyt_2h_mean = as.numeric(as.character(Cyt_2h_mean)),
    Cyt_4h_mean = as.numeric(as.character(Cyt_4h_mean)))

merged_data_GC <- merged_data_GC %>%
  mutate(
    Nuc_A_UT_IR_stability = Nuc_A_UT_1.x / Nuc_UT_mean,
    Nuc_B_UT_IR_stability = Nuc_B_UT_1.x / Nuc_UT_mean,
    
    Nuc_A_30_IR_stability = Nuc_A_30_1.x / Nuc_UT_mean,
    Nuc_B_30_IR_stability = Nuc_B_30_1.x / Nuc_UT_mean,
    
    Nuc_A_2h_IR_stability = Nuc_A_2h_1.x / Nuc_UT_mean,
    Nuc_B_2h_IR_stability = Nuc_B_2h_1.x / Nuc_UT_mean,
    
    Nuc_A_4h_IR_stability = Nuc_A_4h_1.x / Nuc_UT_mean,
    Nuc_B_4h_IR_stability = Nuc_B_4h_1.x / Nuc_UT_mean
  )

# Create new columns comparing the stability of RNA in the nucleus
merged_data_GC_filtered <- merged_data_GC %>%
  mutate(Nuc_UT_TPM_stability = Nuc_UT_TPM_mean / Nuc_UT_TPM_mean,  # Reference of 1 for UT
         Nuc_30_TPM_stability = Nuc_30_TPM_mean / Nuc_UT_TPM_mean,
         Nuc_2h_TPM_stability = Nuc_2h_TPM_mean / Nuc_UT_TPM_mean,
         Nuc_4h_TPM_stability = Nuc_4h_TPM_mean / Nuc_UT_TPM_mean)

# Create new columns comparing the stability of RNA in the cytoplasm
merged_data_GC_filtered <- merged_data_GC_filtered %>%
  mutate(Cyt_UT_mean = as.numeric(as.character(Cyt_UT_mean)),
         Cyt_30_mean = as.numeric(as.character(Cyt_30_mean)),
         Cyt_2h_mean = as.numeric(as.character(Cyt_2h_mean)),
         Cyt_4h_mean = as.numeric(as.character(Cyt_4h_mean)))

merged_data_GC_filtered <- merged_data_GC_filtered %>%
  mutate(Cyt_UT_TPM_stability = Cyt_UT_TPM_mean / Cyt_UT_TPM_mean,  # Reference of 1 for UT
         Cyt_30_TPM_stability = Cyt_30_TPM_mean / Cyt_UT_TPM_mean,
         Cyt_2h_TPM_stability = Cyt_2h_TPM_mean / Cyt_UT_TPM_mean,
         Cyt_4h_TPM_stability = Cyt_4h_TPM_mean / Cyt_UT_TPM_mean)

# Create new columns comparing the stability of IR in the nucleus
merged_data_GC_filtered <- merged_data_GC_filtered %>%
  mutate(Nuc_UT_IR_stability = Nuc_UT_mean / Nuc_UT_mean,  # Reference of 1 for UT
         Nuc_30_IR_stability = Nuc_30_mean / Nuc_UT_mean,
         Nuc_2h_IR_stability = Nuc_2h_mean / Nuc_UT_mean,
         Nuc_4h_IR_stability = Nuc_4h_mean / Nuc_UT_mean)

# Create new columns comparing the stability of IR in the cytoplasm
merged_data_GC_filtered <- merged_data_GC_filtered %>%
  mutate(Cyt_UT_IR_stability = Cyt_UT_mean / Cyt_UT_mean,  # Reference of 1 for UT
         Cyt_30_IR_stability = Cyt_30_mean / Cyt_UT_mean,
         Cyt_2h_IR_stability = Cyt_2h_mean / Cyt_UT_mean,
         Cyt_4h_IR_stability = Cyt_4h_mean / Cyt_UT_mean)

#Filter out rows with 3 or more non-NAs in Nuc IR fractions
merged_data_GC_filtered <- merged_data_GC_filtered %>%
  mutate(non_na_count = rowSums(!is.na(dplyr::select(., Nuc_A_UT_1.x, Nuc_B_UT_1.x, Nuc_A_30_1.x, Nuc_B_30_1.x, Nuc_A_2h_1.x, Nuc_B_2h_1.x, Nuc_A_4h_1.x, Nuc_B_4h_1.x)))) %>%
  filter(non_na_count >= 5)

# For gene_orientated studies sometimes you need to filter the data to keep only the intron with the highest Nuc_UT_mean for each gene
merged_data_GC_filtered_GENE <- merged_data_GC_filtered %>%
  group_by(GENE) %>%
  top_n(1, Nuc_UT_mean) %>%  # Select the top 1 row with the highest Nuc_UT_mean per GENE
  ungroup()

# Assign the stability to introns
#--------------------------------------------#
#### Assign stability, here I care about introns so I dont need to filter overlaping genes
#--------------------------------------------#
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
  ) %>% mutate(stability = "stable 4h")


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
      )) %>%
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
          (Nuc_A_UT_1.x >= 25 & is.na(Nuc_B_UT_1.x)))
  ) %>% mutate(stability = "unstable2h")

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
  ) %>% mutate(stability = "unstable 4h")

noIR <- merged_data_GC_filtered %>%
  filter(Nuc_UT_mean < 30) %>%
  mutate(stability = "no IR")

##no IR GENE
noIR_GENE <- merged_data_GC_filtered_GENE %>%
  filter((Nuc_UT_mean) < 30) %>%
  mutate(stability = "no_IR_GENE")

noIR_GENE <- noIR_GENE %>%
  distinct(GENE, .keep_all = TRUE)
sum(is.na(noIR_GENE$Nuc_UT_mean))

combined_data <- bind_rows(stableR4h, unstableR4h, unstableR30min, unstableR2h, noIR)
combined_dataTPMno100 <- combined_data %>%
  filter(Nuc_UT_TPM_mean >= 5) %>%
  filter(is.na(Cyt_A_UT_1.x) | is.na(Cyt_B_UT_1.x) | #keep is.na because otherwise they get moved
           (Cyt_A_UT_1.x < 99 & Cyt_B_UT_1.x < 99))

#############
#generating stable and unstable GENES, 30% cutoff
#############
#Check overlapping introns between groups
sets_list <- list(
  "Stable 4h" = unique(stableR4h$COORD),
  "Unstable 4h" = unique(unstableR4h$COORD),
  "Unstable 30min" = unique(unstableR30min$COORD),
  "Unstable 2h" = unique(unstableR2h$COORD),
  "No IR" = unique(noIR$COORD))

### Check overlap between GENES in different groups ####
# First, make a list of genes for each stability category of introns
sets_list <- list(
  "Stable 4h" = unique(stableR4h$GENE),
  "Unstable 4h" = unique(unstableR4h$GENE),
  "Unstable 30min" = unique(unstableR30min$GENE),
  "Unstable 2h" = unique(unstableR2h$GENE))

###------------------------------------------------------------------------------####
# ---- Now eliminate overlapping GENES favoring the longer stability group -----######
# ---- Option 1: Extract exclusive genes from each set individually ----
exclusive_by_set <- lapply(names(sets_list), function(set_name) {
  current_set <- sets_list[[set_name]]
  # Combine all the genes in the other sets into one vector:
  other_sets <- unlist(sets_list[names(sets_list) != set_name])
  setdiff(current_set, other_sets)})
names(exclusive_by_set) <- names(sets_list)

# View the exclusive (non-overlapping) genes in each category:
exclusive_by_set  # This is a list where each element contains the non-overlapping genes for that category
exclusive_by_set[["Unstable 30min"]]

### CLEANING over lapping genes ###
exclusive_unstable30min <- setdiff(
  sets_list[["Unstable 30min"]],
  union(union(sets_list[["Stable 4h"]], sets_list[["Unstable 4h"]]), sets_list[["Unstable 2h"]]))

exclusive_unstable2h <- setdiff(
  sets_list[["Unstable 2h"]],
  union(sets_list[["Stable 4h"]], sets_list[["Unstable 4h"]]))

exclusive_unstable4h <- setdiff(unstableR4h$GENE, stableR4h$GENE)

### STEP 2. Check if GENES still overlapp
sets_list_updated <- list(
  "Stable 4h"      = unique(stableR4h$GENE),
  "Unstable 4h"    = unique(exclusive_unstable4h),
  "Unstable 30min" = unique(exclusive_unstable30min),
  "Unstable 2h"    = unique(exclusive_unstable2h))

### STEP 4. Compare the Union of Genes Before and After Filtering
initial_genes <- unique(c(
  stableR4h$GENE,
  unstableR4h$GENE,
  unstableR30min$GENE,
  unstableR2h$GENE))

# For the filtered union, we use the updated sets (note that unstable4h_exclusive is part of the union)
filtered_genes <- unique(c(
  stableR4h$GENE,
  exclusive_unstable4h,              # Already removed genes overlapping stableR4h
  exclusive_unstable30min,
  exclusive_unstable2h))

# Print out the total counts for comparison
cat("Total unique genes before filtering:", length(initial_genes), "\n")
cat("Total unique genes after filtering: ", length(filtered_genes), "\n")

# You can also check that the sets are identical (if the overall union remains unchanged)
if (setequal(initial_genes, filtered_genes)) {
  cat("The union of genes is identical before and after filtering.\n")
} else {
  cat("The union of genes differs between the initial and filtered sets.\n")}

# Create labeled data frame from stableR4h (which is a data frame)
df_stable4h <- stableR4h %>%
  dplyr::select(GENE) %>%
  mutate(GENE_stability = "stable_4h")

# For the others, turn them into data frames directly
df_unstable4h <- data.frame(GENE = exclusive_unstable4h,
                            GENE_stability = "unstable_4h")

df_unstable30min <- data.frame(GENE = exclusive_unstable30min,
                               GENE_stability = "unstable_30min")

df_unstable2h <- data.frame(GENE = exclusive_unstable2h,
                            GENE_stability = "unstable_2h")

# Combine all into one data frame
gene_stability_df <- bind_rows(
  df_stable4h,
  df_unstable4h,
  df_unstable30min,
  df_unstable2h)

gene_stability_df_unique <- gene_stability_df %>%
  distinct(GENE, .keep_all = TRUE)  # Keep only one row per gene

combined_exclusive <- combined_dataTPMno100 %>%
  filter(GENE %in% gene_stability_df_unique$GENE) %>%
  left_join(gene_stability_df_unique, by = "GENE") %>%
  filter(stability != "no IR")

#I dont want to have repeated rows (intros) for the same genes:
combined_exclusive_unique <- combined_exclusive %>%
  distinct(GENE, .keep_all = TRUE)

#### TSA 
sum(duplicated(combined_exclusive_unique$GENE))
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

t


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

#### check if more proximal####
# Make a binary label: proximal (Group_9 or Group_10) vs. all others
# Function to test enrichment in Group_9/10 for any given group vs stable
test_proximal_enrichment <- function(group1, group2, df) {
  df_sub <- df[df$stability %in% c(group1, group2), ]
  tbl <- table(df_sub$stability == group1, df_sub$decile %in% c(9, 10))
  test <- fisher.test(tbl)
  
  pct1 <- mean(df_sub$decile[df_sub$stability == group1] %in% c(9, 10)) * 100
  pct2 <- mean(df_sub$decile[df_sub$stability == group2] %in% c(9, 10)) * 100
  
  data.frame(
    comparison = paste(group1, "vs", group2),
    value = sprintf("Groups 9-10 overlap: %.1f%% vs %.1f%%", pct1, pct2),
    group = paste(group1, "|", group2),
    test = "Fisher exact",
    statistic = unname(test$estimate),   # odds ratio
    `p-value` = test$p.value,
    check.names = FALSE
  )
}

# Run comparisons
fisher_results <- bind_rows(
  test_proximal_enrichment("stable", "unstable 30min", tsa_gene_level),
  test_proximal_enrichment("stable", "unstable 2h", tsa_gene_level),
  test_proximal_enrichment("stable", "unstable 4h", tsa_gene_level)
)

# Add multiple-testing correction
fisher_results$p.adj <- p.adjust(fisher_results$`p-value`, method = "BH")

# Save
write.csv(
  fisher_results,
  "TSAseq_Group9_10_enrichment_vs_stable_corrected.csv",
  row.names = FALSE
)

print(fisher_results)

# Save to CSV
#write.csv(fisher_results, "TSAseq_Group9_10_enrichment_vs_stable.csv", row.names = FALSE)

pairwise.wilcox.test(tsa_gene_level$decile, tsa_gene_level$stability, p.adjust.method = "BH")
# Create a wide-format table: rows = deciles, columns = stability groups, values = percentages
percent_table <- tsa_overlap_all_full %>%
  dplyr::select(stability, group, percent) %>%
  pivot_wider(names_from = stability, values_from = percent)

# View the table
print(percent_table, n = Inf)

# Save as CSV
#write.csv(percent_table, "TSAseq_stability_group_percentages.csv", row.names = FALSE)

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
