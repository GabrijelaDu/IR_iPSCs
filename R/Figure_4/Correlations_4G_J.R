analysis_dir <- here::here("R", "Figure_4")
setwd(analysis_dir)

fig4_corr_path <- file.path(analysis_dir, "Fig_4_correlation_HUVEC_iPS.tsv")
tss_tsa_path <- file.path(analysis_dir, "TSS_TSA_Correlation.tsv")

if (!file.exists(fig4_corr_path)) {
  stop("Missing input file: ", fig4_corr_path, call. = FALSE)
}
if (!file.exists(tss_tsa_path)) {
  stop("Missing input file: ", tss_tsa_path, call. = FALSE)
}

# Load the data
df <- read.delim(fig4_corr_path,
                 header = TRUE,
                 sep = "\t",
                 check.names = FALSE)

# Inspect column names
colnames(df)

# Define the two variables to correlate
x_col <- "IPS_Median_D.SC35"
y_col <- "HUVEC_Median_D.SC35"

# Keep only complete cases
tmp <- df[complete.cases(df[, c(x_col, y_col)]), c(x_col, y_col)]

# Pearson correlation
ct <- cor.test(tmp[[x_col]], tmp[[y_col]], method = "pearson")

# Extract r and calculate R2
r <- unname(ct$estimate)
r2 <- r^2

# Build result table
res <- data.frame(
  X = x_col,
  Y = y_col,
  N = nrow(tmp),
  Pearson_R = r,
  R2 = r2,
  p_value = ct$p.value,
  stringsAsFactors = FALSE
)

# Adjust p-value
# If you only test this one correlation, p_adj will be identical to p_value
res$p_adj <- p.adjust(res$p_value, method = "BH")

# Print result
print(res)

# Optional: save result
write.table(res, "pearson_correlation_results.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)


# Load the file (no header in this file)
df <- read.delim(tss_tsa_path, header = FALSE, sep = "\t", check.names = FALSE)

# Name the columns
colnames(df) <- c("Feature", "Distance", "TSA_decile")

# Convert decimal commas to decimal points, then to numeric
df$Distance <- as.numeric(gsub(",", ".", df$Distance))
df$TSA_decile <- as.numeric(gsub(",", ".", df$TSA_decile))

# Inspect
print(df)

# Keep complete cases only
tmp <- df[complete.cases(df[, c("Distance", "TSA_decile")]), ]

# Pearson correlation
ct <- cor.test(tmp$Distance, tmp$TSA_decile, method = "pearson")

# Build result table
res <- data.frame(
  X = "Distance",
  Y = "TSA_decile",
  N = nrow(tmp),
  Pearson_R = unname(ct$estimate),
  R2 = unname(ct$estimate)^2,
  p_value = ct$p.value,
  stringsAsFactors = FALSE
)

# Adjust p-value
res$p_adj <- p.adjust(res$p_value, method = "BH")

# Print result
print(res)

# Save result
write.table(res,
            "TSS_TSA_correlation_results.tsv",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE)
