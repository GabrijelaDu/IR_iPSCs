suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tibble)
  library(ggplot2)
})

analysis_dir <- here::here("R", "RANDOM_SPOTS")
setwd(analysis_dir)

# =========================================================
# HUVEC IR-RNAs vs RANDOM spots (SC35 & U2)
#
# Primary test:
#   - Build equal-weight mixture null across 10 random nuclei
#   - Convert IR distances to percentiles p = F_mix(d)
#   - Under randomness: p ~ Uniform(0,1)
#   - ONE-SIDED 1-sample KS test (alternative="greater") = enrichment near speckles
#   - BH-FDR on p_ks_closer within each target (SC35, U2)
#
# Output includes:
#   - delta_median and delta_mean vs random mixture expectations
#   - closer/not closer calls based on median and mean shifts (same cutoff)
# =========================================================

# -------------------------
# INPUTS (HUVEC)
# -------------------------
file_path  <- file.path(analysis_dir, "10.07.2025_FINAL_Calculations_Dist_to_Speckles_HUVECs_AC.xlsx")
random_dir <- file.path(analysis_dir, "HUVEC")

if (!file.exists(file_path)) {
  stop("Missing input file: ", file_path, call. = FALSE)
}
if (!dir.exists(random_dir)) {
  stop("Missing random-spot directory: ", random_dir, call. = FALSE)
}

IR_RNAs_to_plot <- c(
  "CENPT", "COG4", "DDX17", "FANCA", "LAMA5", "MEG3", "METTL3",
  "PTBP1", "RAD52_I9", "SFPQ", "SON", "TELO2", "TUG1", "LINC00106"
)

# -------------------------
# RANDOM FILE SETTINGS (HUVEC)
# -------------------------
file_stub <- "HUVEC"
K <- 10
dist_col_random <- "Shortest Distance to Surfaces"

# -------------------------
# STATS SETTINGS
# -------------------------
set.seed(1)

alpha <- 0.05
effect_um <- 0.10              # effect-size cutoff for calling "closer"
distance_threshold_um <- 0.10  # enrichment-within-X metric

# =========================================================
# FUNCTIONS
# =========================================================

read_random_one <- function(fp, target = c("SC35","U2")) {
  target <- match.arg(target)
  dat <- read_excel(fp)
  
  if (!all(c("Category", "Surfaces", dist_col_random) %in% names(dat))) {
    stop("Random file is missing expected columns (Category, Surfaces, distance). File: ",
         basename(fp), "\nColumns are:\n", paste(names(dat), collapse = ", "))
  }
  
  x <- dat %>%
    filter(Category == "Spot", Surfaces == target) %>%
    pull(.data[[dist_col_random]])
  
  x <- x[is.finite(x)]
  if (length(x) < 20) {
    stop("Too few random distances after filtering in: ", basename(fp),
         " (n=", length(x), "). Check filters/column names.")
  }
  x
}

build_mixture_null <- function(target = c("SC35","U2")) {
  target <- match.arg(target)
  
  files <- file.path(random_dir, paste0(file_stub, "_", 1:K, "_", target, ".xls"))
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) {
    stop("Missing random files for ", target, ":\n- ", paste(missing, collapse="\n- "))
  }
  
  rand_by_nuc <- lapply(files, function(fp) read_random_one(fp, target = target))
  ecdfs <- lapply(rand_by_nuc, ecdf)
  
  mix_cdf <- function(x) {
    x <- as.numeric(x)
    rowMeans(sapply(ecdfs, function(Fi) Fi(x)))
  }
  
  # equal-weight sample for plotting + summary expectations
  n_per <- min(5000, min(sapply(rand_by_nuc, length)))
  rand_mix_sample <- unlist(lapply(rand_by_nuc, function(v) sample(v, size = n_per, replace = TRUE)))
  
  list(
    target = target,
    files = files,
    rand_by_nuc = rand_by_nuc,
    mix_cdf = mix_cdf,
    rand_mix_sample = rand_mix_sample
  )
}

plot_random_ecdfs <- function(null_obj) {
  target <- null_obj$target
  rand_by_nuc <- null_obj$rand_by_nuc
  ecdfs <- lapply(rand_by_nuc, ecdf)
  
  all_rand <- unlist(rand_by_nuc)
  xs <- seq(min(all_rand), max(all_rand), length.out = 400)
  
  df <- do.call(rbind, lapply(seq_along(ecdfs), function(i) {
    data.frame(nucleus = paste0("N", i), x = xs, F = ecdfs[[i]](xs))
  }))
  
  mixF <- null_obj$mix_cdf(xs)
  df_mix <- data.frame(nucleus = "Mixture (equal-weight)", x = xs, F = mixF)
  
  ggplot(df, aes(x = x, y = F, group = nucleus, linetype = nucleus)) +
    geom_line(color = "black", linewidth = 0.35) +
    geom_line(data = df_mix, aes(x = x, y = F),
              color = "black", linewidth = 1.0, linetype = "solid", inherit.aes = FALSE) +
    theme_classic(base_size = 12) +
    labs(
      title = paste0("HUVEC random ECDFs per nucleus: ", target),
      subtitle = "Thin: each nucleus; Thick: equal-weight mixture null",
      x = "Distance to speckle surface (µm)",
      y = "ECDF",
      linetype = NULL
    ) +
    guides(linetype = "none")
}

analyze_rna_vs_null <- function(rna_vec, null_obj, RNA, target,
                                distance_threshold_um = 0.10) {
  rna_vec <- rna_vec[is.finite(rna_vec)]
  n <- length(rna_vec)
  if (n < 3) return(NULL)
  
  # percentile transform vs random mixture null
  p <- null_obj$mix_cdf(rna_vec)
  
  # ONE-SIDED KS test for enrichment near speckles
  ksg <- suppressWarnings(ks.test(p, "punif", alternative = "greater"))
  
  # optional diagnostics (kept)
  ks2 <- suppressWarnings(ks.test(p, "punif", alternative = "two.sided"))
  ksl <- suppressWarnings(ks.test(p, "punif", alternative = "less"))
  
  # effect sizes
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

print(plot_random_ecdfs(null_SC35))
print(plot_random_ecdfs(null_U2))

# =========================================================
# RUN ANALYSIS FOR ALL RNAs
# =========================================================
out <- list()

for (sheet in IR_RNAs_to_plot) {
  df <- read_excel(file_path, sheet = sheet)
  
  if ("SC35" %in% names(df)) {
    v <- df$SC35
    res <- analyze_rna_vs_null(v, null_SC35, RNA = sheet, target = "SC35",
                               distance_threshold_um = distance_threshold_um)
    if (!is.null(res)) out[[paste0(sheet, "_SC35")]] <- res
  }
  
  if ("U2" %in% names(df)) {
    v <- df$U2
    res <- analyze_rna_vs_null(v, null_U2, RNA = sheet, target = "U2",
                               distance_threshold_um = distance_threshold_um)
    if (!is.null(res)) out[[paste0(sheet, "_U2")]] <- res
  }
}

res_df <- bind_rows(out) %>% arrange(target, p_ks_closer)

# FDR correction within target (PRIMARY: one-sided closer p-values)
res_df <- res_df %>%
  group_by(target) %>%
  mutate(
    fdr_ks_closer = p.adjust(p_ks_closer, method = "BH")
  ) %>%
  ungroup()

# Pretty printing for tiny p-values (avoids showing 0)
res_df <- res_df %>%
  mutate(
    p_ks_closer_txt   = format.pval(p_ks_closer, digits = 3, eps = 1e-300),
    fdr_ks_closer_txt = format.pval(fdr_ks_closer, digits = 3, eps = 1e-300)
  )

# Add BOTH median-based and mean-based "closer / not closer" reports
res_df <- res_df %>%
  mutate(
    closer_median_report = ifelse(fdr_ks_closer < alpha & delta_median <= -effect_um,
                                  "closer", "not closer"),
    closer_mean_report   = ifelse(fdr_ks_closer < alpha & delta_mean <= -effect_um,
                                  "closer", "not closer"),
    
    class_median = case_when(
      fdr_ks_closer >= alpha ~ "indistinguishable from random",
      fdr_ks_closer <  alpha & delta_median <= -effect_um ~ "non-random: closer",
      fdr_ks_closer <  alpha & abs(delta_median) < effect_um ~ paste0("sig, tiny median shift (<", effect_um, " µm)"),
      TRUE ~ "other"
    ),
    class_mean = case_when(
      fdr_ks_closer >= alpha ~ "indistinguishable from random",
      fdr_ks_closer <  alpha & delta_mean <= -effect_um ~ "non-random: closer",
      fdr_ks_closer <  alpha & abs(delta_mean) < effect_um ~ paste0("sig, tiny mean shift (<", effect_um, " µm)"),
      TRUE ~ "other"
    )
  ) %>%
  arrange(target, fdr_ks_closer, RNA)

print(res_df)

# =========================================================
# OPTIONAL: density overlays IR vs Random (mixture)
# =========================================================
dens_u2 <- make_density_df("U2", null_U2)

# define RNA order from U2
u2_order <- order_by_median(dens_u2)

dens_u2 <- dens_u2 %>%
  mutate(RNA = factor(RNA, levels = u2_order))

p_u2_dens <- ggplot(dens_u2, aes(x = distance, linetype = group)) +
  geom_density(color = "black", linewidth = 0.8, adjust = 1) +
  theme_classic(base_size = 12) +
  facet_wrap(~RNA, ncol = 4) +
  labs(
    title = "U2: IR-RNA vs Random (mixture null), ordered by IR median",
    subtitle = paste0("Primary test: one-sided KS on percentiles; B=", B),
    x = "Distance to speckle surface (µm)",
    y = "Density",
    linetype = NULL
  )

dens_sc35 <- make_density_df("SC35", null_SC35) %>%
  mutate(RNA = factor(RNA, levels = u2_order))

p_sc35_dens <- ggplot(dens_sc35, aes(x = distance, linetype = group)) +
  geom_density(color = "black", linewidth = 0.8, adjust = 1) +
  theme_classic(base_size = 12) +
  facet_wrap(~RNA, ncol = 4) +
  labs(
    title = "SC35: IR-RNA vs Random (mixture null), ordered as in U2",
    subtitle = paste0("Primary test: one-sided KS on percentiles; B=", B),
    x = "Distance to speckle surface (µm)",
    y = "Density",
    linetype = NULL
  )

print(p_u2_dens)
print(p_sc35_dens)

# Save SC35 density plot
ggsave(
  filename = "p_sc35_dens_huvec.pdf",
  plot = p_sc35_dens,
  width = 6,
  height = 6
)

ggsave(
  filename = "p_u2_dens_huvec.pdf",
  plot = p_u2_dens,
  width = 6,
  height = 6
)

###### OVERLAY plot #####

dens_u2 <- make_density_df("U2", null_U2)

# define RNA order from U2
u2_order <- order_by_median(dens_u2)

dens_u2 <- dens_u2 %>%
  mutate(RNA = factor(RNA, levels = u2_order))

dens_sc35 <- make_density_df("SC35", null_SC35) %>%
  mutate(RNA = factor(RNA, levels = u2_order))

# Split IR-RNA and random
ir_u2 <- dens_u2 %>% filter(group == "IR-RNA")
ir_sc35 <- dens_sc35 %>% filter(group == "IR-RNA")

rand_u2 <- tibble(distance = null_U2$rand_mix_sample)
rand_sc35 <- tibble(distance = null_SC35$rand_mix_sample)

# palettes
green_pal <- colorRampPalette(c("#0B3D2E", "#1B5E20", "#2E7D32", "#66BB6A", "#C8E6C9"))(length(u2_order))
names(green_pal) <- u2_order

red_pal <- colorRampPalette(c("#4A0D0D", "#7F0000", "#B22222", "#E57373", "#FAD4D4"))(length(u2_order))
names(red_pal) <- u2_order

# U2 = red
p_u2_dens <- ggplot() +
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
  xlim(-0.5, 1.5) +
  theme_classic(base_size = 12) +
  labs(
    title = "U2: all IR-RNAs overlaid with random spots",
    subtitle = paste0("Primary test: one-sided KS on percentiles; B=", B),
    x = "Distance to speckle surface (µm)",
    y = "Density",
    colour = "RNA",
    fill = "RNA"
  )

# SC35 = green
p_sc35_dens <- ggplot() +
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
  xlim(-0.5, 1.5) +
  theme_classic(base_size = 12) +
  labs(
    title = "SC35: all IR-RNAs overlaid with random spots",
    subtitle = paste0("Primary test: one-sided KS on percentiles; B=", B),
    x = "Distance to speckle surface (µm)",
    y = "Density",
    colour = "RNA",
    fill = "RNA"
  )

print(p_u2_dens)
print(p_sc35_dens)

ggsave("p_u2_dens_overlay_red_HUVEC.pdf", plot = p_u2_dens, width = 8, height = 6)
ggsave("p_sc35_dens_overlay_green_HUVEC.pdf", plot = p_sc35_dens, width = 8, height = 6)

library(patchwork)

p_side_by_side <- p_u2_dens + p_sc35_dens +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "right")

print(p_side_by_side)

ggsave(
  "p_u2_sc35_side_by_side_HUVEC_for_legend.pdf",
  plot = p_side_by_side,
  width = 16,
  height = 16
)

# =========================================================
# SAVE RESULTS
# =========================================================
output_path <- file.path(analysis_dir, "HUVEC_IRRNA_equalweight_mixture_results_KS.csv")
write.csv(res_df, output_path, row.names = FALSE)
