suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

analysis_dir <- here::here("R", "Suppl_Fig_7")
dir.create(analysis_dir, recursive = TRUE, showWarnings = FALSE)
setwd(analysis_dir)

input_path <- file.path(analysis_dir, "Correlation_Suppl_Fig. 7B.tsv")
if (!file.exists(input_path)) {
  stop("Missing input file: ", input_path, call. = FALSE)
}

# Load the data
df <- read.delim(input_path,
                 header = TRUE,
                 sep = "\t",
                 check.names = FALSE)

# Inspect column names
colnames(df)

# Split into 3 clean tables
df_intron <- df[, c(1, 2, 3)]
colnames(df_intron) <- c("Feature", "Distance_sc35", "logFC_SON")

df_exon <- df[, c(5, 6, 7)]
colnames(df_exon) <- c("Feature", "Distance_sc35", "logFC_SON")

df_total <- df[, c(9, 10, 11)]
colnames(df_total) <- c("Feature", "Distance_sc35", "logFC_SON")

# Function to run Pearson correlation
run_corr <- function(dat, label) {
  dat$Distance_sc35 <- as.numeric(gsub(",", ".", dat$Distance_sc35))
  dat$logFC_SON <- as.numeric(gsub(",", ".", dat$logFC_SON))
  
  tmp <- dat[complete.cases(dat[, c("Distance_sc35", "logFC_SON")]), ]
  
  ct <- cor.test(tmp$Distance_sc35, tmp$logFC_SON, method = "pearson")
  
  data.frame(
    comparison = label,
    N = nrow(tmp),
    Pearson_R = unname(ct$estimate),
    R2 = unname(ct$estimate)^2,
    p_value = ct$p.value,
    stringsAsFactors = FALSE
  )
}

# Run correlations
res <- bind_rows(
  run_corr(df_intron, "SC35 distance vs SON intron logFC"),
  run_corr(df_exon,   "SC35 distance vs SON exon logFC"),
  run_corr(df_total,  "SC35 distance vs SON total logFC")
)

# Adjust across the 3 tests
res$p_adj_BH <- p.adjust(res$p_value, method = "BH")

# Print results
print(res)

# Save results
write.table(
  res,
  "SC35_SON_correlation_results_Suppl_Fig_7B.tsv",
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Recreate clean data for plotting
df_intron <- df[, c(1, 2, 3)]
colnames(df_intron) <- c("Feature", "Distance_sc35", "logFC_SON")

df_exon <- df[, c(5, 6, 7)]
colnames(df_exon) <- c("Feature", "Distance_sc35", "logFC_SON")

df_total <- df[, c(9, 10, 11)]
colnames(df_total) <- c("Feature", "Distance_sc35", "logFC_SON")

# Convert to numeric and keep complete cases
clean_corr_data <- function(dat) {
  dat$Distance_sc35 <- as.numeric(gsub(",", ".", dat$Distance_sc35))
  dat$logFC_SON <- as.numeric(gsub(",", ".", dat$logFC_SON))
  dat <- dat[complete.cases(dat[, c("Distance_sc35", "logFC_SON")]), ]
  return(dat)
}

df_intron_clean <- clean_corr_data(df_intron)
df_exon_clean   <- clean_corr_data(df_exon)
df_total_clean  <- clean_corr_data(df_total)

# Plotting function
make_corr_plot <- function(dat, xlab, ylab, title_prefix) {
  ct <- cor.test(dat$Distance_sc35, dat$logFC_SON, method = "pearson")
  r_val <- signif(unname(ct$estimate), 3)
  p_val <- format(ct$p.value, scientific = TRUE, digits = 3)
  
  ggplot(dat, aes(x = Distance_sc35, y = logFC_SON)) +
    geom_point(size = 2, alpha = 0.8) +
    geom_smooth(method = "lm", se = TRUE) +
    labs(
      title = paste0(title_prefix, "\nPearson r = ", r_val, ", p = ", p_val),
      x = xlab,
      y = ylab
    ) +
    theme_minimal(base_size = 14)
}

# Create plots
p_intron <- make_corr_plot(
  df_intron_clean,
  xlab = "Median distance to SC35",
  ylab = "log2FC SON intron",
  title_prefix = "Distance to SC35 vs SON intron signal"
)

p_exon <- make_corr_plot(
  df_exon_clean,
  xlab = "Median distance to SC35",
  ylab = "log2FC SON exon",
  title_prefix = "Distance to SC35 vs SON exon signal"
)

p_total <- make_corr_plot(
  df_total_clean,
  xlab = "Median distance to SC35",
  ylab = "log2FC SON total",
  title_prefix = "Distance to SC35 vs SON total signal"
)

# Show plots
print(p_intron)
print(p_exon)
print(p_total)

# Save plots
ggsave("SC35_vs_SON_intron_correlation.pdf", plot = p_intron, width = 6, height = 5)
ggsave("SC35_vs_SON_exon_correlation.pdf", plot = p_exon, width = 6, height = 5)
ggsave("SC35_vs_SON_total_correlation.pdf", plot = p_total, width = 6, height = 5)
