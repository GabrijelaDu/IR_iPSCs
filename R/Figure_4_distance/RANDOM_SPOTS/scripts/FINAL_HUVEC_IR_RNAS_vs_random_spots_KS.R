suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(patchwork)
})

# =========================================================
# HUVEC IR-RNAs vs RANDOM spots (SC35 & U2)
#
# Primary test:
#   - Build equal-weight mixture null across 10 random nuclei
#   - Convert IR distances to percentiles p = F_mix(d)
#   - Under randomness: p ~ Uniform(0,1)
#   - One-sided 1-sample KS test (alternative = "greater")
#     to test enrichment near speckles
#   - BH-FDR correction on p_ks_closer within each target
#     (SC35, U2)
#
# Output includes:
#   - delta_median and delta_mean vs random mixture expectation
#   - closer / not closer calls based on median and mean shifts
#   - overlay density plots for U2 and SC35
#
# Expected repository structure:
# Figure_4_distance/
#   └── RANDOM_SPOTS/
#       ├── HUVEC/
#       ├── iPS/
#       ├── results/
#       ├── scripts/
#       │   ├── FINAL_HUVEC_IR_RNAs_vs_random_spots_KS.R
#       │   └── FINAL_iPS_IR_RNAs_vs_random_spots_KS.R
#       ├── Distance_to_speckles_HUVECs.xlsx
#       └── Distance_to_speckles_iPSCs.xlsx
#
# This script assumes it is run from:
# RANDOM_SPOTS/scripts/
# =========================================================

# -------------------------
# INPUTS
# -------------------------
file_path  <- "../Distance_to_speckles_HUVECs.xlsx"
random_dir <- "../HUVEC"
output_dir <- "../results"

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

IR_RNAs_to_plot <- c(
  "CENPT", "COG4", "DDX17", "FANCA", "LAMA5", "MEG3", "METTL3",
  "PTBP1", "RAD52", "SFPQ", "SON", "TELO2", "TUG1", "LINC00106"
)

if (!file.exists(file_path)) {
  stop("Input file not found: ", file_path)
}
if (!dir.exists(random_dir)) {
  stop("Random directory not found: ", random_dir)
}

available_sheets <- excel_sheets(file_path)
missing_sheets <- setdiff(IR_RNAs_to_plot, available_sheets)

if (length(missing_sheets) > 0) {
  warning(
    "These sheets were not found and will be skipped: ",
    paste(missing_sheets, collapse = ", ")
  )
}

IR_RNAs_to_plot <- intersect(IR_RNAs_to_plot, available_sheets)

# -------------------------
# RANDOM FILE SETTINGS
# -------------------------
file_stub <- "HUVEC"
K <- 10
dist_col_random <- "Shortest Distance to Surfaces"

# -------------------------
# STATS SETTINGS
# -------------------------
set.seed(1)

alpha <- 0.05
effect_um <- 0.10
distance_threshold_um <- 0.10

# =========================================================
# FUNCTIONS
# =========================================================

read_random_one <- function(fp, target = c("SC35", "U2")) {
  target <- match.arg(target)
  dat <- read_excel(fp)
  
  required_cols <- c("Category", "Surfaces", dist_col_random)
  
  if (!all(required_cols %in% names(dat))) {
    stop(
      "Random file is missing expected columns. File: ", basename(fp),
      "\nExpected columns: ", paste(required_cols, collapse = ", "),
      "\nFound columns: ", paste(names(dat), collapse = ", ")
    )
  }
  
  x <- dat %>%
    filter(Category == "Spot", Surfaces == target) %>%
    pull(.data[[dist_col_random]])
  
  x <- x[is.finite(x)]
  
  if (length(x) < 20) {
    stop(
      "Too few random distances after filtering in: ",
      basename(fp), " (n = ", length(x), ")"
    )
  }
  
  x
}

build_mixture_null <- function(target = c("SC35", "U2")) {
  target <- match.arg(target)
  
  files <- file.path(random_dir, paste0(file_stub, "_", 1:K, "_", target, ".xls"))
  missing <- files[!file.exists(files)]
  
  if (length(missing) > 0) {
    stop(
      "Missing random files for ", target, ":\n- ",
      paste(missing, collapse = "\n- ")
    )
  }
  
  rand_by_nuc <- lapply(files, function(fp) read_random_one(fp, target = target))
  ecdfs <- lapply(rand_by_nuc, ecdf)
  
  mix_cdf <- function(x) {
    x <- as.numeric(x)
    rowMeans(sapply(ecdfs, function(Fi) Fi(x)))
  }
  
  n_per <- min(5000, min(sapply(rand_by_nuc, length)))
  rand_mix_sample <- unlist(
    lapply(rand_by_nuc, function(v) sample(v, size = n_per, replace = TRUE))
  )
  
  list(
    target = target,
    files = files,
    rand_by_nuc = rand_by_nuc,
    mix_cdf = mix_cdf,
    rand_mix_sample = rand_mix_sample
  )
}

analyze_rna_vs_null <- function(rna_vec, null_obj, RNA, target,
                                distance_threshold_um = 0.10) {
  rna_vec <- rna_vec[is.finite(rna_vec)]
  n <- length(rna_vec)
  
  if (n < 3) {
    return(NULL)
  }
  
  p <- null_obj$mix_cdf(rna_vec)
  
  ksg <- suppressWarnings(ks.test(p, "punif", alternative = "greater"))
  ks2 <- suppressWarnings(ks.test(p, "punif", alternative = "two.sided"))
  ksl <- suppressWarnings(ks.test(p, "punif", alternative = "less"))
  
  med_dist <- median(rna_vec)
  mean_dist <- mean(rna_vec)
  
  rand_samp <- null_obj$rand_mix_sample
  rand_med <- median(rand_samp)
  rand_mean <- mean(rand_samp)
  
  ir_frac <- mean(rna_vec <= distance_threshold_um)
  rand_frac <- mean(rand_samp <= distance_threshold_um)
  enrich_within_X <- ifelse(rand_frac > 0, ir_frac / rand_frac, NA_real_)
  
  tibble(
    RNA = RNA,
    target = target,
    n_rna = n,
    
    rna_median = med_dist,
    rand_median_exp = rand_med,
    delta_median = med_dist - rand_med,
    
    rna_mean = mean_dist,
    rand_mean_exp = rand_mean,
    delta_mean = mean_dist - rand_mean,
    
    ks_stat = as.numeric(ks2$statistic),
    p_ks_two_sided = as.numeric(ks2$p.value),
    p_ks_closer = as.numeric(ksg$p.value),
    p_ks_farther = as.numeric(ksl$p.value),
    
    frac_within_X = ir_frac,
    rand_frac_within_X = rand_frac,
    enrich_within_X = enrich_within_X,
    X_um = distance_threshold_um
  )
}

make_density_df <- function(target_col, null_obj) {
  out <- list()
  
  for (sheet in IR_RNAs_to_plot) {
    df <- read_excel(file_path, sheet = sheet)
    
    if (!target_col %in% names(df)) next
    
    v <- df[[target_col]]
    v <- v[is.finite(v)]
    
    if (length(v) < 3) next
    
    out[[sheet]] <- bind_rows(
      tibble(RNA = sheet, group = "IR-RNA", distance = v),
      tibble(RNA = sheet, group = "Random (mixture)", distance = null_obj$rand_mix_sample)
    )
  }
  
  bind_rows(out)
}

order_by_median <- function(plot_df) {
  plot_df %>%
    filter(group == "IR-RNA") %>%
    group_by(RNA) %>%
    summarise(rna_median = median(distance, na.rm = TRUE), .groups = "drop") %>%
    arrange(rna_median) %>%
    pull(RNA)
}

# =========================================================
# BUILD NULLS
# =========================================================
null_SC35 <- build_mixture_null("SC35")
null_U2   <- build_mixture_null("U2")

# =========================================================
# RUN ANALYSIS FOR ALL RNAs
# =========================================================
out <- list()

for (sheet in IR_RNAs_to_plot) {
  df <- read_excel(file_path, sheet = sheet)
  
  if ("SC35" %in% names(df)) {
    res <- analyze_rna_vs_null(
      rna_vec = df$SC35,
      null_obj = null_SC35,
      RNA = sheet,
      target = "SC35",
      distance_threshold_um = distance_threshold_um
    )
    if (!is.null(res)) {
      out[[paste0(sheet, "_SC35")]] <- res
    }
  }
  
  if ("U2" %in% names(df)) {
    res <- analyze_rna_vs_null(
      rna_vec = df$U2,
      null_obj = null_U2,
      RNA = sheet,
      target = "U2",
      distance_threshold_um = distance_threshold_um
    )
    if (!is.null(res)) {
      out[[paste0(sheet, "_U2")]] <- res
    }
  }
}

res_df <- bind_rows(out) %>%
  arrange(target, p_ks_closer)

res_df <- res_df %>%
  group_by(target) %>%
  mutate(
    fdr_ks_closer = p.adjust(p_ks_closer, method = "BH")
  ) %>%
  ungroup()

res_df <- res_df %>%
  mutate(
    p_ks_closer_txt   = format.pval(p_ks_closer, digits = 3, eps = 1e-300),
    fdr_ks_closer_txt = format.pval(fdr_ks_closer, digits = 3, eps = 1e-300),
    
    closer_median_report = ifelse(
      fdr_ks_closer < alpha & delta_median <= -effect_um,
      "closer", "not closer"
    ),
    closer_mean_report = ifelse(
      fdr_ks_closer < alpha & delta_mean <= -effect_um,
      "closer", "not closer"
    ),
    
    class_median = case_when(
      fdr_ks_closer >= alpha ~ "indistinguishable from random",
      fdr_ks_closer < alpha & delta_median <= -effect_um ~ "non-random: closer",
      fdr_ks_closer < alpha & abs(delta_median) < effect_um ~ paste0("sig, tiny median shift (<", effect_um, " µm)"),
      TRUE ~ "other"
    ),
    class_mean = case_when(
      fdr_ks_closer >= alpha ~ "indistinguishable from random",
      fdr_ks_closer < alpha & delta_mean <= -effect_um ~ "non-random: closer",
      fdr_ks_closer < alpha & abs(delta_mean) < effect_um ~ paste0("sig, tiny mean shift (<", effect_um, " µm)"),
      TRUE ~ "other"
    )
  ) %>%
  arrange(target, fdr_ks_closer, RNA)

print(res_df)

write.csv(
  res_df,
  file.path(output_dir, "HUVEC_IRRNA_equalweight_mixture_results_KS.csv"),
  row.names = FALSE
)

# =========================================================
# OVERLAY DENSITY PLOTS
# =========================================================
dens_u2 <- make_density_df("U2", null_U2)
u2_order <- order_by_median(dens_u2)

dens_u2 <- dens_u2 %>%
  mutate(RNA = factor(RNA, levels = u2_order))

dens_sc35 <- make_density_df("SC35", null_SC35) %>%
  mutate(RNA = factor(RNA, levels = u2_order))

ir_u2 <- dens_u2 %>% filter(group == "IR-RNA")
ir_sc35 <- dens_sc35 %>% filter(group == "IR-RNA")

rand_u2 <- tibble(distance = null_U2$rand_mix_sample)
rand_sc35 <- tibble(distance = null_SC35$rand_mix_sample)

green_pal <- colorRampPalette(c("#0B3D2E", "#1B5E20", "#2E7D32", "#66BB6A", "#C8E6C9"))(length(u2_order))
names(green_pal) <- u2_order

red_pal <- colorRampPalette(c("#4A0D0D", "#7F0000", "#B22222", "#E57373", "#FAD4D4"))(length(u2_order))
names(red_pal) <- u2_order

p_u2_overlay <- ggplot() +
  geom_density(
    data = ir_u2,
    aes(x = distance, colour = RNA, fill = RNA, group = RNA),
    linewidth = 0.8,
    alpha = 0.3
  ) +
  geom_density(
    data = rand_u2,
    aes(x = distance),
    colour = "black",
    fill = "black",
    alpha = 0.15,
    linewidth = 1
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6, colour = "black") +
  scale_colour_manual(values = red_pal) +
  scale_fill_manual(values = red_pal) +
  coord_cartesian(xlim = c(-0.5, 1.5)) +
  theme_classic(base_size = 12) +
  labs(
    title = "U2: all IR-RNAs overlaid with random spots",
    x = "Distance to speckle surface (µm)",
    y = "Density",
    colour = "RNA",
    fill = "RNA"
  )

p_sc35_overlay <- ggplot() +
  geom_density(
    data = ir_sc35,
    aes(x = distance, colour = RNA, fill = RNA, group = RNA),
    linewidth = 0.8,
    alpha = 0.3
  ) +
  geom_density(
    data = rand_sc35,
    aes(x = distance),
    colour = "black",
    fill = "black",
    alpha = 0.15,
    linewidth = 1
  ) +
  geom_vline(xintercept = 0, linetype = "dashed", linewidth = 0.6, colour = "black") +
  scale_colour_manual(values = green_pal) +
  scale_fill_manual(values = green_pal) +
  coord_cartesian(xlim = c(-0.5, 1.5)) +
  theme_classic(base_size = 12) +
  labs(
    title = "SC35: all IR-RNAs overlaid with random spots",
    x = "Distance to speckle surface (µm)",
    y = "Density",
    colour = "RNA",
    fill = "RNA"
  )

p_overlay_side_by_side <- p_u2_overlay + p_sc35_overlay +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

ggsave(
  filename = file.path(output_dir, "p_u2_dens_overlay_red_HUVEC.pdf"),
  plot = p_u2_overlay,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(output_dir, "p_sc35_dens_overlay_green_HUVEC.pdf"),
  plot = p_sc35_overlay,
  width = 8,
  height = 6
)

ggsave(
  filename = file.path(output_dir, "p_u2_sc35_side_by_side_HUVEC.pdf"),
  plot = p_overlay_side_by_side,
  width = 16,
  height = 6
)

message("Saved results to: ", output_dir)