# ───────────────────────────────────────────────────────────────────────────────
# BOX PLOTS ONLY — one complete script (NO patchwork)
# X order: for each gene -> NT, siSON, siSRRM2, siSRRM2+siSON, then next gene
# Genes shown: 4 genes (TERT_I11, CENPT, RAD52_I10, COG4)
# Style: KEEP your theme (no grey grid lines, black axes, x-axis ticks kept)
# Jitter: reduced
# Boxplot: pure standard Tukey boxplot; whisker end ticks more visible via linewidth
# ADD: Exon nuc/cyt abundance in the SAME one-line style
# ───────────────────────────────────────────────────────────────────────────────

library(tidyverse)
library(tidyr)

analysis_dir <- here::here("R", "Figure_5")
setwd(analysis_dir)

input_path <- file.path(analysis_dir, "son_srrm2_quantif_exp_unified_R.csv")
if (!file.exists(input_path)) {
  stop("Missing input file: ", input_path, call. = FALSE)
}

# ── Load cleaned CSV file ─────────────────────────────────────────────────────
df_si <- read.csv(
  input_path,
  check.names = FALSE
)

# Convert all columns to numeric (safe coercion)
df_si[] <- lapply(df_si, function(x) suppressWarnings(as.numeric(x)))

# ── Long format + metadata ────────────────────────────────────────────────────
df_si_long <- df_si %>%
  pivot_longer(
    cols      = everything(),
    names_to  = "variable",
    values_to = "value"
  ) %>%
  mutate(
    Feature = case_when(
      str_detect(variable, "_PIR$") ~ "PIR",
      str_detect(variable, "_Ex")   ~ "Exon",
      str_detect(variable, "_I")    ~ "Intron",
      TRUE ~ NA_character_
    ),
    Gene = str_extract(variable, "^[A-Z0-9]+(_I\\d+)?"),
    Treatment = case_when(
      str_detect(variable, "_SON_SRRM2_") ~ "siSRRM2+siSON",
      str_detect(variable, "_SRRM2_")     ~ "siSRRM2",
      str_detect(variable, "_SON_")       ~ "siSON",
      TRUE ~ "NT"
    ),
    Compartment = case_when(
      str_detect(variable, "_PIR$")  ~ "total",
      str_detect(variable, "_Nuc_")  ~ "nuc",
      str_detect(variable, "_Cyto_") ~ "cyt",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(Feature %in% c("Exon", "Intron", "PIR")) %>%
  filter(!is.na(Gene), !is.na(Treatment), !is.na(value))

# ── Global factor levels (treatment order fixed) ──────────────────────────────
treat_order <- c("NT", "siSON", "siSRRM2", "siSRRM2+siSON")
df_si_long <- df_si_long %>%
  mutate(Treatment = factor(Treatment, levels = treat_order))

# ───────────────────────────────────────────────────────────────────────────────
# BOXplot function (genes in ONE ROW, shared y, separators between genes)
# IMPORTANT: x-axis order is forced: Gene-major then Treatment within Gene
# Pure boxplot only, thicker whiskers/caps via linewidth
# ───────────────────────────────────────────────────────────────────────────────
plot_box_feature <- function(feature_label,
                             y_title,
                             compartment = NULL,
                             ylim_override = NULL,
                             genes = NULL,
                             gene_order = c("TERT_I11", "CENPT", "RAD52_I10", "COG4"),
                             jitter_width = 0.08,
                             box_linewidth = 0.9) {
  
  df_plot <- df_si_long %>%
    filter(Feature == feature_label) %>%
    { if (!is.null(compartment)) filter(., Compartment == compartment) else . } %>%
    { if (!is.null(genes))       filter(., as.character(Gene) %in% genes) else . } %>%
    mutate(
      Gene      = factor(as.character(Gene), levels = gene_order),
      Treatment = factor(as.character(Treatment), levels = treat_order)
    ) %>%
    filter(!is.na(Gene)) %>%
    droplevels()
  
  if (nrow(df_plot) == 0) stop("No data after filtering (feature/compartment/genes).")
  
  # Explicit x levels: for each gene, all treatments
  x_levels <- unlist(lapply(gene_order, function(g) paste(g, treat_order)))
  df_plot <- df_plot %>%
    mutate(
      x_group_chr = paste(as.character(Gene), as.character(Treatment)),
      x_group     = factor(x_group_chr, levels = x_levels)
    )
  
  # separators between genes (thin black)
  n_treat <- length(treat_order)
  boundaries <- if (length(gene_order) > 1) {
    seq(from = n_treat + 0.5, by = n_treat, length.out = length(gene_order) - 1)
  } else numeric(0)
  
  p <- ggplot(df_plot, aes(x = x_group, y = value)) +
    geom_boxplot(
      outlier.shape = NA,
      fill      = "grey90",
      color     = "black",
      width     = 0.5,
      linewidth = box_linewidth,  # ⬅ thicker whiskers + whisker ticks
      coef      = 1.5
    ) +
    geom_jitter(
      width = jitter_width, height = 0,
      alpha = 0.7, color = "#990066", size = 1.3
    ) +
    { if (length(boundaries) > 0) geom_vline(xintercept = boundaries, linewidth = 0.25) } +
    labs(
      title = paste0(feature_label, if (!is.null(compartment)) paste0(" [", compartment, "]")),
      y = y_title,
      x = "Treatment"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(colour = "black"),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      axis.ticks.y     = element_line(colour = "black"),
      axis.ticks.x     = element_line(colour = "black")
    )
  
  if (!is.null(ylim_override)) {
    p <- p + scale_y_continuous(limits = ylim_override)
  } else {
    p <- p + scale_y_continuous(limits = c(0, NA))
  }
  
  p
}

# ───────────────────────────────────────────────────────────────────────────────
# FEATURE BOX PLOTS (4 genes on one row; each gene has NT, siSON, siSRRM2, combo)
# ───────────────────────────────────────────────────────────────────────────────
gene_order_4 <- c("TERT_I11", "CENPT", "RAD52_I10", "COG4")

plot_exon_nuc_4  <- plot_box_feature("Exon",   "Exon (nuc)",   compartment = "nuc",
                                     genes = gene_order_4, gene_order = gene_order_4)

plot_exon_cyt_4  <- plot_box_feature("Exon",   "Exon (cyt)",   compartment = "cyt",
                                     genes = gene_order_4, gene_order = gene_order_4)

plot_intron_nuc_4 <- plot_box_feature("Intron", "Intron (nuc)", compartment = "nuc",
                                      genes = gene_order_4, gene_order = gene_order_4)

plot_pir_4 <- plot_box_feature("PIR", "PIR (%)", compartment = "total",
                               genes = gene_order_4, gene_order = gene_order_4,
                               ylim_override = c(0, 100))

ggsave("BOX_4genes_Exon_nuc_revised.pdf",   plot_exon_nuc_4,   width = 11, height = 4)
ggsave("BOX_4genes_Exon_cyt_revised.pdf",   plot_exon_cyt_4,   width = 11, height = 4)
ggsave("BOX_4genes_Intron_nuc_revised.pdf", plot_intron_nuc_4, width = 11, height = 4)
ggsave("BOX_4genes_PIR_revised.pdf",        plot_pir_4,        width = 11, height = 4)

# ───────────────────────────────────────────────────────────────────────────────
# EXON NUC/CYT ABUNDANCE (boxplots) — SAME one-line style as above
# shared y 0..1, same x order, same separators, same theme
# ───────────────────────────────────────────────────────────────────────────────

# Keep only Exon rows
df_exon <- df_si_long %>% filter(Feature == "Exon")

# Replicate id per Gene/Treatment/Compartment
df_exon <- df_exon %>%
  group_by(Gene, Treatment, Compartment) %>%
  mutate(rep_id = row_number()) %>%
  ungroup()

# Wide nuc/cyt per replicate
exon_abund_wide <- df_exon %>%
  select(Gene, Treatment, rep_id, Compartment, value) %>%
  pivot_wider(
    names_from = Compartment,
    values_from = value,
    values_fill = 0
  ) %>%
  mutate(
    Total          = nuc + cyt,
    Nuc_abundance  = if_else(Total > 0, nuc / Total, NA_real_),
    Cyto_abundance = if_else(Total > 0, cyt / Total, NA_real_)
  ) %>%
  mutate(
    Gene      = factor(as.character(Gene), levels = gene_order_4),
    Treatment = factor(as.character(Treatment), levels = treat_order)
  ) %>%
  filter(!is.na(Gene)) %>%
  droplevels()

# Explicit x levels (same as feature plots)
x_levels_abund <- unlist(lapply(gene_order_4, function(g) paste(g, treat_order)))
exon_abund_wide <- exon_abund_wide %>%
  mutate(
    x_group_chr = paste(as.character(Gene), as.character(Treatment)),
    x_group     = factor(x_group_chr, levels = x_levels_abund)
  )

abund_boxplot_one_line <- function(dat, yvar, ylab, title_txt,
                                   jitter_width = 0.08,
                                   box_linewidth = 0.9) {
  
  # separators between genes
  n_treat <- length(treat_order)
  boundaries <- if (length(gene_order_4) > 1) {
    seq(from = n_treat + 0.5, by = n_treat, length.out = length(gene_order_4) - 1)
  } else numeric(0)
  
  ggplot(dat, aes(x = x_group, y = .data[[yvar]])) +
    geom_boxplot(
      outlier.shape = NA,
      width     = 0.55,
      fill      = "grey90",
      color     = "black",
      linewidth = box_linewidth,  # ⬅ thicker whiskers + ticks
      coef      = 1.5
    ) +
    geom_jitter(width = jitter_width, alpha = 0.7, color = "#990066", size = 1.2) +
    { if (length(boundaries) > 0) geom_vline(xintercept = boundaries, linewidth = 0.25) } +
    labs(title = title_txt, y = ylab, x = NULL) +
    scale_y_continuous(limits = c(0, 1)) +
    theme_minimal(base_size = 14) +
    theme(
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line        = element_line(colour = "black"),
      axis.text.x      = element_text(angle = 45, hjust = 1),
      axis.ticks.y     = element_line(colour = "black"),
      axis.ticks.x     = element_line(colour = "black")
    ) +
    theme(legend.position = "none")
}

plot_abund_nuc_4 <- abund_boxplot_one_line(
  exon_abund_wide, "Nuc_abundance", "Nuclear fraction",
  "Nuclear exon abundance (box): 4 genes"
)

plot_abund_cyt_4 <- abund_boxplot_one_line(
  exon_abund_wide, "Cyto_abundance", "Cytoplasmic fraction",
  "Cytoplasmic exon abundance (box): 4 genes"
)

ggsave("BOX_4genes_Nuc_abundance_oneLine_revised.pdf", plot_abund_nuc_4, width = 11, height = 4)
ggsave("BOX_4genes_Cyt_abundance_oneLine_revised.pdf", plot_abund_cyt_4, width = 11, height = 4)


library(patchwork)

library(patchwork)

combined_5_vertical <-
  plot_exon_nuc_4 /
  plot_exon_cyt_4 /
  plot_intron_nuc_4 /
  plot_pir_4 /
  plot_abund_nuc_4 +
  plot_layout(ncol = 1)

ggsave("BOX_4genes_ALL_vertical_withNucAbundance.pdf",
       combined_5_vertical,
       width = 11,
       height = 20)

# ───────────────────────────────────────────────────────────────────────────────
# TERT_I8 ONLY — same one-line style, separate outputs
# (Assumes df_si_long, treat_order, and plot_box_feature() already exist)
# ───────────────────────────────────────────────────────────────────────────────

gene_i8 <- c("TERT_I8")

plot_i8_exon_nuc  <- plot_box_feature(
  feature_label = "Exon",
  y_title       = "Nuclear exon spots (count)",
  compartment   = "nuc",
  genes         = gene_i8,
  gene_order    = gene_i8,
  jitter_width  = 0.08,
  box_linewidth = 0.9
)

plot_i8_exon_cyt  <- plot_box_feature(
  feature_label = "Exon",
  y_title       = "Cytoplasmic exon spots (count)",
  compartment   = "cyt",
  genes         = gene_i8,
  gene_order    = gene_i8,
  jitter_width  = 0.08,
  box_linewidth = 0.9
)

plot_i8_intron_nuc <- plot_box_feature(
  feature_label = "Intron",
  y_title       = "Nuclear retained intron spots (count)",
  compartment   = "nuc",
  genes         = gene_i8,
  gene_order    = gene_i8,
  jitter_width  = 0.08,
  box_linewidth = 0.9
)

plot_i8_pir <- plot_box_feature(
  feature_label = "PIR",
  y_title       = "PIR (%)",
  compartment   = "total",
  genes         = gene_i8,
  gene_order    = gene_i8,
  ylim_override = c(0, 100),
  jitter_width  = 0.08,
  box_linewidth = 0.9
)
combined_I8_vertical <-
  plot_i8_exon_nuc /
  plot_i8_exon_cyt /
  plot_i8_intron_nuc /
  plot_i8_pir +
  plot_layout(ncol = 1)

ggsave("BOX_I8_vertical_withNucAbundance.pdf",
       combined_I8_vertical,
       width = 11,
       height = 20)
#ggsave("BOX_TERT_I8_Exon_nuc.pdf",   plot_i8_exon_nuc,   width = 4, height = 4)
#ggsave("BOX_TERT_I8_Exon_cyt.pdf",   plot_i8_exon_cyt,   width = 4, height = 4)
#ggsave("BOX_TERT_I8_Intron_nuc.pdf", plot_i8_intron_nuc, width = 4, height = 4)
#ggsave("BOX_TERT_I8_PIR.pdf",        plot_i8_pir,        width = 4, height = 4)

# ───────────────────────────────────────────────────────────────────────────────
# TERT_I8 exon nuc/cyt fraction — same style (0..1)
# ───────────────────────────────────────────────────────────────────────────────

df_exon_i8 <- df_si_long %>%
  filter(Feature == "Exon", Gene == "TERT_I8", Compartment %in% c("nuc", "cyt")) %>%
  group_by(Gene, Treatment, Compartment) %>%
  mutate(rep_id = row_number()) %>%
  ungroup()

exon_i8_wide <- df_exon_i8 %>%
  select(Gene, Treatment, rep_id, Compartment, value) %>%
  pivot_wider(names_from = Compartment, values_from = value, values_fill = 0) %>%
  mutate(
    Total          = nuc + cyt,
    Nuc_fraction   = if_else(Total > 0, nuc / Total, NA_real_),
    Cyto_fraction  = if_else(Total > 0, cyt / Total, NA_real_)
  ) %>%
  mutate(
    Treatment = factor(as.character(Treatment), levels = treat_order),
    x_group   = factor(as.character(Treatment), levels = treat_order)
  )

plot_i8_nuc_frac <- ggplot(exon_i8_wide, aes(x = x_group, y = Nuc_fraction)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", color = "black",
               width = 0.55, linewidth = 0.9, coef = 1.5) +
  geom_jitter(width = 0.08, height = 0, alpha = 0.7, color = "#990066", size = 1.2) +
  labs(title = "TERT_I8", y = "Nuclear transcript fraction", x = "Treatment") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks.y     = element_line(colour = "black"),
    axis.ticks.x     = element_line(colour = "black")
  )

plot_i8_cyt_frac <- ggplot(exon_i8_wide, aes(x = x_group, y = Cyto_fraction)) +
  geom_boxplot(outlier.shape = NA, fill = "grey90", color = "black",
               width = 0.55, linewidth = 0.9, coef = 1.5) +
  geom_jitter(width = 0.08, height = 0, alpha = 0.7, color = "#990066", size = 1.2) +
  labs(title = "TERT_I8", y = "Cytoplasmic transcript fraction", x = "Treatment") +
  scale_y_continuous(limits = c(0, 1)) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line        = element_line(colour = "black"),
    axis.text.x      = element_text(angle = 45, hjust = 1),
    axis.ticks.y     = element_line(colour = "black"),
    axis.ticks.x     = element_line(colour = "black")
  )
