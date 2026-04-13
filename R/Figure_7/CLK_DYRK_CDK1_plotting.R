### This script plots Clk DYRK CDK1 treatments for our IR-RNAs in interphase and mitosis

# Load required libraries
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)

analysis_dir <- here::here("R", "Figure_7")
setwd(analysis_dir)

# Set the Excel file path
file_path <- file.path(analysis_dir, "Quantification_CDK_CLK_08092025_forR_FINAL6_JB2.xlsx")
if (!file.exists(file_path)) {
  stop("Missing input file: ", file_path, call. = FALSE)
}

# Read the data (first sheet)
df_CLK_CDK_DYRK <- read_excel(file_path)
head(df_CLK_CDK_DYRK)

# Step 1: Select metaphase columns for NT, CDK1, CLK_DYRK, exclude 3XD
metaphase_cols <- grep("_M_(NT|CDK1|Clk_DYRK)_.*total_PIR$", names(df_CLK_CDK_DYRK), value = TRUE)
metaphase_cols <- metaphase_cols[!grepl("3XD", metaphase_cols)]

# Step 2: Subset and reshape
df_m <- df_CLK_CDK_DYRK %>% select(all_of(metaphase_cols))

df_m_long <- df_m %>%
  pivot_longer(cols = everything(), names_to = "Condition", values_to = "total_PIR") %>%
  mutate(
    Gene = sub("_M_.*", "", Condition),
    Treatment = case_when(
      grepl("_M_NT_", Condition) ~ "NT",
      grepl("_M_CDK1_", Condition) ~ "CDK1",
      grepl("_M_Clk_DYRK_", Condition) ~ "CLK_DYRK"
    )
  )

# Step 3: Set treatment order manually
df_m_long$Treatment <- factor(df_m_long$Treatment, levels = c("NT", "CDK1", "CLK_DYRK"))

library(dplyr)
library(ggplot2)
library(patchwork)

# Palette
pal_teal <- c(NT = "#005f73", CDK1 = "#0a9396", CLK_DYRK = "#6FB3A8")
pal_purple <- c(NT = "#5e3c99", CDK1 = "#8c6bb1", CLK_DYRK = "#d0c1e1")
pal_slate <- c(NT = "#37474f", CDK1 = "#78909c", CLK_DYRK = "#cfd8dc")
pal_indigo <- c(NT = "#2c3e99", CDK1 = "#5b74c6", CLK_DYRK = "#b0bdf0")

# (Optional) enforce Treatment order
df_m_long <- df_m_long %>%
  mutate(Treatment = factor(Treatment, levels = c("NT","CDK1","CLK_DYRK")))

# Function to make a single gene plot (uses palette)
make_gene_plot <- function(gene_name) {
  ggplot(df_m_long %>% filter(Gene == gene_name),
         aes(x = Treatment, y = total_PIR, color = Treatment, fill = Treatment)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2, show.legend = FALSE) +
    geom_boxplot(alpha = 0.25, width = 0.5, outlier.shape = NA, show.legend = FALSE) +
    scale_color_manual(values = pal_teal) +
    scale_fill_manual(values = pal_teal) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    theme_minimal() +
    labs(title = gene_name, x = "Treatment", y = "total_PIR") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      text = element_text(size = 12),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.6)
    )
}

# List of unique genes and assemble grid
genes <- unique(df_m_long$Gene)
plots <- lapply(genes, make_gene_plot)
final_plot <- wrap_plots(plots, nrow = 3, ncol = 2)

final_plot


ggplot2::ggsave("CDK1_CLK_DYRK_mitosis.pdf", final_plot, width = 8, height = 9)  # letter size


#### interphase #####

# Step 1: Select interphase total_PIR columns (NT, CDK1, CLK_DYRK), excluding 3XD
interphase_cols <- grep("_Inter_(NT|CDK1|Clk_DYRK)_.*total_PIR$", names(df_CLK_CDK_DYRK), value = TRUE)
interphase_cols <- interphase_cols[!grepl("3XD", interphase_cols)]

# Step 2: Subset and reshape
df_inter <- df_CLK_CDK_DYRK %>% select(all_of(interphase_cols))

df_inter_long <- df_inter %>%
  pivot_longer(cols = everything(), names_to = "Condition", values_to = "total_PIR") %>%
  mutate(
    Gene = sub("_Inter_.*", "", Condition),
    Treatment = case_when(
      grepl("_Inter_NT_", Condition) ~ "NT",
      grepl("_Inter_CDK1_", Condition) ~ "CDK1",
      grepl("_Inter_Clk_DYRK_", Condition) ~ "CLK_DYRK"
    )
  )

# Step 3: Set treatment order
df_inter_long$Treatment <- factor(df_inter_long$Treatment, levels = c("NT", "CDK1", "CLK_DYRK"))

# Step 4: Plot with facet wrap (3 rows × 2 columns)
library(dplyr)
library(ggplot2)
library(patchwork)

# optional: consistent treatment colors (light-green forward)
pal <- c("NT"="#2e7d32","CDK1"="#66bb6a","CLK_DYRK"="#a5d6a7")

make_inter_plot <- function(gene_name) {
  ggplot(df_inter_long %>% filter(Gene == gene_name),
         aes(x = Treatment, y = total_PIR, color = Treatment)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2) +
    geom_boxplot(alpha = 0.2, width = 0.5, outlier.shape = NA) +
    scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
    # comment out next line if you don't want custom colors
    scale_color_manual(values = pal) +
    theme_minimal() +
    labs(title = gene_name, x = "Treatment", y = "total_PIR") +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      axis.text.x = element_text(angle = 45, hjust = 1),
      text = element_text(size = 12),
      legend.position = "none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.6)
    )
}

# Pick the six genes you want (order matters). For all genes, use unique(df_inter_long$Gene)
genes_to_plot <- df_inter_long %>%
  distinct(Gene) %>%
  arrange(Gene) %>%      # or set your preferred order
  pull(Gene) %>%
  head(6)

plots <- lapply(genes_to_plot, make_inter_plot)

p_inter <- wrap_plots(plots, nrow = 3, ncol = 2) +
  plot_annotation(title = "Interphase total_PIR across Treatments") &
  theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"))

p_inter


ggplot2::ggsave("CDK1_CLK_DYRK_interphase.pdf", p_inter, width = 8, height = 9)  # letter size


########## AVERAGE VALUES ############

### Mitosis: calculate ΔPIR = Treatment - NT

# --- your existing code up to df_mitosis_long ---

# --- compute per-gene p-values: Treatment vs NT ---
# Welch t-test per gene for CDK1 vs NT and CLK_DYRK vs NT
# --- compute per-gene RAW p-values: Treatment vs NT (Welch unpaired t-test) ---

# helper: format number as mantissa (1 decimal) x10-EXP
format_sci1 <- function(x) {
  if (is.na(x)) return("n/a")
  if (x == 0) return("0.0")
  e <- floor(log10(abs(x)))
  m <- round(x / (10^e), 1)
  # handle rounding to 10.0
  if (m >= 10) { m <- 1.0; e <- e + 1 }
  paste0(format(m, nsmall = 1, trim = TRUE), "x10", e)  # e prints like -6
}

# --- Select mitosis PIR columns (exclude 3XD) ---
mitosis_cols <- grep("_M_(NT|CDK1|Clk_DYRK)_.*total_PIR$", names(df_CLK_CDK_DYRK), value = TRUE)
mitosis_cols <- mitosis_cols[!grepl("3XD", mitosis_cols)]

df_mitosis <- df_CLK_CDK_DYRK %>%
  select(all_of(mitosis_cols))

# --- Long format with Gene and Treatment ---
df_mitosis_long <- df_mitosis %>%
  pivot_longer(cols = everything(), names_to = "Condition", values_to = "PIR") %>%
  mutate(
    Gene = sub("_M_.*", "", Condition),
    Treatment = case_when(
      grepl("_M_NT_", Condition) ~ "NT",
      grepl("_M_CDK1_", Condition) ~ "CDK1",
      grepl("_M_Clk_DYRK_", Condition) ~ "CLK_DYRK",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Treatment))

# --- RAW p-values (Welch unpaired t-test) vs NT, per Gene ---
df_pvals <- df_mitosis_long %>%
  filter(Treatment %in% c("NT", "CDK1", "CLK_DYRK")) %>%
  group_by(Gene) %>%
  group_modify(~{
    dat <- .x
    cmp <- function(tr) {
      x <- dat$PIR[dat$Treatment == tr]
      y <- dat$PIR[dat$Treatment == "NT"]
      p <- if (length(x) >= 2 && length(y) >= 2)
        tryCatch(t.test(x, y, var.equal = FALSE)$p.value, error = function(e) NA_real_)
      else NA_real_
      tibble(Treatment = tr, p_value = p)  # RAW p
    }
    bind_rows(cmp("CDK1"), cmp("CLK_DYRK"))
  }) %>%
  ungroup()

# --- Mean PIRs, deltas vs NT, join raw p-values ---
df_mitosis_delta <- df_mitosis_long %>%
  group_by(Gene, Treatment) %>%
  summarise(mean_PIR = mean(PIR, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Treatment, values_from = mean_PIR) %>%
  mutate(
    delta_CDK1 = CDK1 - NT,
    delta_CLK_DYRK = `CLK_DYRK` - NT
  ) %>%
  pivot_longer(cols = starts_with("delta_"), names_to = "Treatment", values_to = "delta_PIR") %>%
  mutate(Treatment = gsub("delta_", "", Treatment)) %>%
  left_join(df_pvals, by = c("Gene", "Treatment")) %>%
  mutate(
    p_label = ifelse(is.na(p_value), "n/a", paste0("p=", format_sci1(p_value)))
  )

# --- Label positions ---
label_offset <- max(abs(df_mitosis_delta$delta_PIR), na.rm = TRUE) * 0.05
df_mitosis_delta <- df_mitosis_delta %>%
  mutate(y_pos = delta_PIR + ifelse(delta_PIR >= 0, label_offset, -label_offset))

# --- Optional: order genes by mean absolute delta for nicer plotting ---
gene_order <- df_mitosis_delta %>%
  group_by(Gene) %>%
  summarise(mag = mean(abs(delta_PIR), na.rm = TRUE), .groups = "drop") %>%
  arrange(desc(mag)) %>%
  pull(Gene)

df_mitosis_delta <- df_mitosis_delta %>%
  mutate(Gene = factor(Gene, levels = gene_order))

# --- Plot ΔPIR with RAW p-value labels (1-dec mantissa, x10-EXP) ---
ggplot(df_mitosis_delta, aes(x = Gene, y = delta_PIR, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(y = y_pos, label = p_label),
            position = position_dodge(width = 0.7),
            vjust = ifelse(df_mitosis_delta$delta_PIR >= 0, -0.2, 1.2),
            size = 3.8) +
  labs(title = "ΔPIR in Mitosis (Treatment vs NT)",
       y = "ΔPIR (mean)", x = "") +
  scale_fill_manual(values = c("CDK1" = "#2e7d32", "CLK_DYRK" = "#81c784")) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5, size = 16),
    text = element_text(size = 14)
  ) +
  coord_cartesian(ylim = range(df_mitosis_delta$delta_PIR, na.rm = TRUE) + c(-label_offset*3, label_offset*3))


####CONFIDENCE INTERVALS MITOSIS ####
######################################
####

# Build per-gene Welch tests returning p and 95% CI for (Treatment - NT)
df_tests <- df_mitosis_long %>%
  filter(Treatment %in% c("NT","CDK1","CLK_DYRK")) %>%
  group_by(Gene) %>%
  group_modify(~{
    dat <- .x
    test_one <- function(tr) {
      x <- dat$PIR[dat$Treatment == tr]
      y <- dat$PIR[dat$Treatment == "NT"]
      if (length(x) >= 2 && length(y) >= 2) {
        res <- tryCatch(t.test(x, y, var.equal = FALSE), error = function(e) NULL)
        if (!is.null(res)) {
          return(tibble(
            Treatment = tr,
            p_value   = res$p.value,
            ci_low    = unname(res$conf.int[1]),
            ci_high   = unname(res$conf.int[2]),
            n_trt     = length(x),
            n_nt      = length(y)
          ))
        }
      }
      tibble(Treatment = tr, p_value = NA_real_, ci_low = NA_real_, ci_high = NA_real_,
             n_trt = length(x), n_nt = length(y))
    }
    dplyr::bind_rows(test_one("CDK1"), test_one("CLK_DYRK"))
  }) %>% ungroup()

# Join CIs to your delta table (df_mitosis_delta already has delta_PIR and p_label)
df_mitosis_delta2 <- df_mitosis_delta %>%
  select(Gene, Treatment, delta_PIR, p_label) %>%
  left_join(df_tests, by = c("Gene","Treatment"))

# Plot: mean Δ with 95% CI error bars + p labels
ggplot(df_mitosis_delta2, aes(x = Gene, y = delta_PIR, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_errorbar(aes(ymin = ci_low, ymax = ci_high),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.5) +
  geom_text(aes(y = delta_PIR, label = p_label),
            position = position_dodge(width = 0.7), vjust = -1.0, size = 3.6) +
  labs(title = "ΔPIR in Mitosis (Treatment vs NT)",
       y = "ΔPIR (mean ± 95% CI)", x = "") +
  scale_fill_manual(values = c("CDK1" = "#0a9396", "CLK_DYRK" = "#6FB3A8")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title  = element_text(hjust = 0.5, size = 16),
        text        = element_text(size = 14)) +
  coord_cartesian(ylim = range(c(df_mitosis_delta2$ci_low, df_mitosis_delta2$ci_high), na.rm = TRUE) * c(1,1.05))

ggplot2::ggsave("mitosis_delta_ci.pdf", plot = ggplot2::last_plot(),
                width = 8, height = 6, units = "in", device = cairo_pdf)
###### SAME CALCULUS FOR INTERPHASE ######

##############################################
# -------- Interphase ΔPIR (Treatment vs NT) with raw Welch p-values --------
# helper for labels
format_sci1 <- function(x) {
  if (is.na(x)) return("n/a")
  if (x == 0) return("0.0")
  e <- floor(log10(abs(x)))
  m <- round(x / (10^e), 1)
  if (m >= 10) { m <- 1.0; e <- e + 1 }
  paste0(format(m, nsmall = 1, trim = TRUE), "x10", e)  # e like -6
}

# 1) Select interphase columns (NT, CDK1, CLK_DYRK), exclude 3XD
inter_cols <- grep("_Inter_(NT|CDK1|Clk_DYRK)_.*total_PIR$", names(df_CLK_CDK_DYRK), value = TRUE)
inter_cols <- inter_cols[!grepl("3XD", inter_cols)]
stopifnot(length(inter_cols) > 0)

# 2) Subset & long
df_inter <- df_CLK_CDK_DYRK %>% dplyr::select(all_of(inter_cols))

df_inter_long <- df_inter %>%
  tidyr::pivot_longer(cols = dplyr::everything(), names_to = "Condition", values_to = "PIR") %>%
  dplyr::mutate(
    Gene = sub("_Inter_.*", "", Condition),
    Treatment = dplyr::case_when(
      grepl("_Inter_NT_", Condition) ~ "NT",
      grepl("_Inter_CDK1_", Condition) ~ "CDK1",
      grepl("_Inter_Clk_DYRK_", Condition) ~ "CLK_DYRK",
      TRUE ~ NA_character_
    )
  ) %>%
  dplyr::filter(!is.na(Treatment)) %>%
  dplyr::mutate(Treatment = factor(Treatment, levels = c("NT","CDK1","CLK_DYRK")))

# 3) Per-gene raw p-values (Welch unpaired t-test) vs NT
df_pvals_inter <- df_inter_long %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(
    p_CDK1 = {
      x <- PIR[Treatment == "CDK1"]; y <- PIR[Treatment == "NT"]
      if (length(x) >= 2 && length(y) >= 2) tryCatch(t.test(x,y, var.equal = FALSE)$p.value, error = function(e) NA_real_) else NA_real_
    },
    p_CLK_DYRK = {
      x <- PIR[Treatment == "CLK_DYRK"]; y <- PIR[Treatment == "NT"]
      if (length(x) >= 2 && length(y) >= 2) tryCatch(t.test(x,y, var.equal = FALSE)$p.value, error = function(e) NA_real_) else NA_real_
    },
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(cols = dplyr::starts_with("p_"), names_to = "Treatment", values_to = "p_value") %>%
  dplyr::mutate(Treatment = dplyr::recode(Treatment, p_CDK1 = "CDK1", p_CLK_DYRK = "CLK_DYRK"))

# 4) Means, deltas, join p-values
df_inter_delta <- df_inter_long %>%
  dplyr::group_by(Gene, Treatment) %>%
  dplyr::summarise(mean_PIR = mean(PIR, na.rm = TRUE), .groups = "drop") %>%
  tidyr::pivot_wider(names_from = Treatment, values_from = mean_PIR) %>%
  dplyr::mutate(
    delta_CDK1 = CDK1 - NT,
    delta_CLK_DYRK = `CLK_DYRK` - NT
  ) %>%
  tidyr::pivot_longer(dplyr::starts_with("delta_"), names_to = "Treatment", values_to = "delta_PIR") %>%
  dplyr::mutate(Treatment = sub("^delta_", "", Treatment)) %>%
  dplyr::left_join(df_pvals_inter, by = c("Gene", "Treatment")) %>%
  dplyr::mutate(p_label = ifelse(is.na(p_value), "n/a", paste0("p=", format_sci1(p_value))))

# 5) Label positions & order
label_offset_inter <- max(abs(df_inter_delta$delta_PIR), na.rm = TRUE) * 0.05
df_inter_delta <- df_inter_delta %>%
  dplyr::mutate(y_pos = delta_PIR + ifelse(delta_PIR >= 0, label_offset_inter, -label_offset_inter))

gene_order_inter <- df_inter_delta %>%
  dplyr::group_by(Gene) %>%
  dplyr::summarise(mag = mean(abs(delta_PIR), na.rm = TRUE), .groups = "drop") %>%
  dplyr::arrange(dplyr::desc(mag)) %>% dplyr::pull(Gene)

df_inter_delta <- df_inter_delta %>%
  dplyr::mutate(Gene = factor(Gene, levels = gene_order_inter))

# 6) Plot (Interphase)
p_inter_delta <- ggplot(df_inter_delta, aes(x = Gene, y = delta_PIR, fill = Treatment)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(y = y_pos, label = p_label),
            position = position_dodge(width = 0.7),
            vjust = ifelse(df_inter_delta$delta_PIR >= 0, -0.2, 1.2),
            size = 3.8) +
  labs(title = "ΔPIR in Interphase (Treatment vs NT)",
       y = "ΔPIR (mean)", x = "") +
  scale_fill_manual(values = c("CDK1" = "#0a9396", "CLK_DYRK" = "#6FB3A8")) +
  theme_minimal() +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1),
    plot.title   = element_text(hjust = 0.5, size = 16),
    text         = element_text(size = 14)
  ) +
  coord_cartesian(ylim = range(df_inter_delta$delta_PIR, na.rm = TRUE) + c(-label_offset_inter*3, label_offset_inter*3))

p_inter_delta

# Save interphase plot
ggplot2::ggsave("CDK1_CLK_DYRK_interphase_delta.pdf", p_inter_delta, width = 8, height = 9)


###### MEAN VALUES CONSIDERING INDEPENDENT REPLICAS / DATES OF EXPERIMENT ########

# Define PIR columns (mitosis only, 3 treatments, 5 genes)
pir_cols <- c(
  "CENPT_M_NT_total_PIR", "CENPT_M_CDK1_total_PIR", "CENPT_M_Clk_DYRK_total_PIR",
  "METTL3_M_NT_total_PIR", "METTL3_M_CDK1_total_PIR", "METTL3_M_Clk_DYRK_total_PIR",
  "TERT_M_NT_total_PIR", "TERT_M_CDK1_total_PIR", "TERT_M_Clk_DYRK_total_PIR",
  "Brd8_M_NT_total_PIR", "Brd8_M_CDK1_total_PIR", "Brd8_M_Clk_DYRK_total_PIR",
  "FANCA_M_NT_total_PIR", "FANCA_M_CDK1_total_PIR", "FANCA_M_Clk_DYRK_total_PIR"
)

# Derive corresponding replica column names
replica_cols <- gsub("total_PIR", "replica", pir_cols)

# Create long-format dataframe with PIR and replica values
df_long <- tibble()

for (i in seq_along(pir_cols)) {
  pir_col <- pir_cols[i]
  replica_col <- replica_cols[i]
  
  temp <- df_CLK_CDK_DYRK %>%
    select(all_of(c(pir_col, replica_col))) %>%
    rename(PIR = !!pir_col, replica = !!replica_col) %>%
    mutate(condition = pir_col)
  
  df_long <- bind_rows(df_long, temp)
}

# Extract Gene, Compartment (M), and Treatment from condition name
df_long <- df_long %>%
  separate(condition, into = c("Gene", "Compartment", "Treatment", "Measure"),
           sep = "_", extra = "merge") %>%
  select(Gene, Compartment, Treatment, PIR, replica)

df_long <- df_long %>%
  filter(!is.na(replica), !is.na(PIR))

# Calculate mean PIR per replica
df_replica_mean <- df_long %>%
  group_by(Gene, Treatment, Compartment, replica) %>%
  summarise(mean_PIR = mean(PIR, na.rm = TRUE), .groups = "drop")

# Calculate mean and standard error across replicates
df_plot <- df_replica_mean %>%
  group_by(Gene, Treatment, Compartment) %>%
  summarise(
    PIR_mean = mean(mean_PIR),
    PIR_sd = sd(mean_PIR),
    n = n(),
    PIR_se = PIR_sd / sqrt(n),
    .groups = "drop"
  )

# Set treatment order manually
df_plot$Treatment <- factor(df_plot$Treatment, levels = c("NT", "CDK1", "Clk"))

# Set treatment order manually
df_replica_mean$Treatment <- factor(df_replica_mean$Treatment, 
                                    levels = c("NT", "CDK1", "Clk"))

# Plot: individual replicate means + overall mean ± SE
ggplot(df_replica_mean %>% filter(Compartment == "M"),
       aes(x = Treatment, y = mean_PIR, color = Treatment)) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +   # individual replicate points
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4, 
               color = "black", fill = "black") +      # black diamond = mean
  stat_summary(fun.data = mean_se, geom = "errorbar", 
               width = 0.2, color = "black") +         # error bars = SE
  facet_wrap(~ Gene, nrow = 3, ncol = 2) +
  labs(x = "Treatment", y = "PIR (%)", title = "PIR in Mitosis across Treatments") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"))


###### DELTA PIR for matching replicas only ####
# Check matching replicas for mitosis
replica_check <- df_replica_mean %>%
  filter(Compartment == "M") %>%
  group_by(Gene, replica) %>%
  summarise(treatments = n_distinct(Treatment),
            list_treatments = paste(sort(unique(Treatment)), collapse = ", "),
            .groups = "drop") %>%
  arrange(Gene, replica)

replica_check

# --- up to df_replica_mean you keep your code exactly the same ---

# Calculate delta PIR relative to NT for each gene/replica
df_delta <- df_replica_mean %>%
  filter(Compartment == "M") %>%
  pivot_wider(names_from = Treatment, values_from = mean_PIR) %>%
  mutate(
    delta_CDK1 = CDK1 - NT,
    delta_Clk  = Clk  - NT
  ) %>%
  select(Gene, replica, delta_CDK1, delta_Clk) %>%
  pivot_longer(cols = starts_with("delta"),
               names_to = "Treatment",
               values_to = "delta_PIR") %>%
  mutate(Treatment = recode(Treatment,
                            "delta_CDK1" = "CDK1",
                            "delta_Clk"  = "Clk"))

# Plot delta PIR (points + mean ± SE)
ggplot(df_delta, aes(x = Treatment, y = delta_PIR, color = Treatment)) +
  geom_jitter(width = 0.15, size = 2, alpha = 0.8) +
  stat_summary(fun = mean, geom = "point", shape = 18, size = 4,
               color = "black", fill = "black") +
  stat_summary(fun.data = mean_se, geom = "errorbar",
               width = 0.2, color = "black") +
  facet_wrap(~ Gene, nrow = 3, ncol = 2) +
  labs(x = "Treatment", y = "ΔPIR vs NT (%)",
       title = "ΔPIR in Mitosis relative to NT") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"))


###### SOME REPLICAS ARE NOT MATCHING SO COMPARING TO MEAN NT ######

## 1) Get the NT reference (per gene, mitosis only)
nt_ref <- df_replica_mean %>%
  filter(Compartment == "M", Treatment == "NT") %>%
  group_by(Gene) %>%
  summarise(
    NT_mean = mean(mean_PIR, na.rm = TRUE),
    n_NT = n(),
    .groups = "drop"
  )

## 2) Compute ΔPIR for each treatment replicate vs NT_mean
df_delta_ntmean <- df_replica_mean %>%
  filter(Compartment == "M", Treatment != "NT") %>%
  left_join(nt_ref, by = "Gene") %>%
  mutate(delta_PIR = mean_PIR - NT_mean) %>%
  filter(!is.na(NT_mean))  # drop genes with no NT reference

## 3) Order treatments as requested
df_delta_ntmean$Treatment <- factor(df_delta_ntmean$Treatment,
                                    levels = c("CDK1", "Clk"))

## 4) Summary for error bars (mean ± SE of ΔPIR)
df_delta_summary <- df_delta_ntmean %>%
  group_by(Gene, Treatment) %>%
  summarise(
    delta_mean = mean(delta_PIR, na.rm = TRUE),
    delta_sd   = sd(delta_PIR, na.rm = TRUE),
    n          = dplyr::n(),
    delta_se   = delta_sd / sqrt(n),
    .groups = "drop"
  )

## 5) Plot: replicate points + overall mean ± SE (per gene), NT not shown
ggplot(df_delta_ntmean,
       aes(x = Treatment, y = delta_PIR, color = Treatment)) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  geom_jitter(width = 0.15, size = 2, alpha = 0.85) +
  # mean points (black diamonds)
  geom_point(data = df_delta_summary,
             aes(x = Treatment, y = delta_mean),
             shape = 18, size = 4, color = "black", inherit.aes = FALSE) +
  # error bars (SE)
  geom_errorbar(data = df_delta_summary,
                aes(x = Treatment,
                    ymin = delta_mean - delta_se,
                    ymax = delta_mean + delta_se),
                width = 0.2, color = "black", inherit.aes = FALSE) +
  facet_wrap(~ Gene, nrow = 3, ncol = 2) +
  labs(x = "Treatment",
       y = "ΔPIR vs NT mean (%)",
       title = "Mitosis: ΔPIR per replicate vs across-replicate NT mean") +
  theme_minimal() +
  theme(strip.text = element_text(size = 12, face = "bold"),
        legend.position = "none")

####### The same for INETRPHASE comparing each replica to NT mean ###########

# Define PIR columns (interphase only, 3 treatments, 5 genes)
pir_cols_inter <- c(
  "CENPT_Inter_NT_total_PIR", "CENPT_Inter_CDK1_total_PIR", "CENPT_Inter_Clk_DYRK_total_PIR",
  "METTL3_Inter_NT_total_PIR", "METTL3_Inter_CDK1_total_PIR", "METTL3_Inter_Clk_DYRK_total_PIR",
  "TERT_Inter_NT_total_PIR",   "TERT_Inter_CDK1_total_PIR",   "TERT_Inter_Clk_DYRK_total_PIR",
  "Brd8_Inter_NT_total_PIR",   "Brd8_Inter_CDK1_total_PIR",   "Brd8_Inter_Clk_DYRK_total_PIR",
  "FANCA_Inter_NT_total_PIR",  "FANCA_Inter_CDK1_total_PIR",  "FANCA_Inter_Clk_DYRK_total_PIR"
)

# Replica cols
replica_cols_inter <- gsub("total_PIR", "replica", pir_cols_inter)

# Build long-format dataframe for interphase
df_long_inter <- tibble()
for (i in seq_along(pir_cols_inter)) {
  pir_col <- pir_cols_inter[i]
  replica_col <- replica_cols_inter[i]
  
  temp <- df_CLK_CDK_DYRK %>%
    select(all_of(c(pir_col, replica_col))) %>%
    rename(PIR = !!pir_col, replica = !!replica_col) %>%
    mutate(condition = pir_col)
  
  df_long_inter <- bind_rows(df_long_inter, temp)
}

# Parse condition into Gene, Compartment, Treatment
df_long_inter <- df_long_inter %>%
  separate(condition, into = c("Gene", "Compartment", "Treatment", "Measure"),
           sep = "_", extra = "merge") %>%
  select(Gene, Compartment, Treatment, PIR, replica) %>%
  filter(!is.na(replica), !is.na(PIR))

# Mean PIR per replica
df_replica_mean_inter <- df_long_inter %>%
  group_by(Gene, Treatment, Compartment, replica) %>%
  summarise(mean_PIR = mean(PIR, na.rm = TRUE), .groups = "drop")

## 1) Get the NT reference (per gene, interphase only)
nt_ref_inter <- df_replica_mean_inter %>%
  filter(Compartment == "Inter", Treatment == "NT") %>%
  group_by(Gene) %>%
  summarise(
    NT_mean = mean(mean_PIR, na.rm = TRUE),
    n_NT = n(),
    .groups = "drop"
  )

## 2) Compute ΔPIR for each treatment replicate vs NT_mean
df_delta_ntmean_inter <- df_replica_mean_inter %>%
  filter(Compartment == "Inter", Treatment != "NT") %>%
  left_join(nt_ref_inter, by = "Gene") %>%
  mutate(delta_PIR = mean_PIR - NT_mean) %>%
  filter(!is.na(NT_mean))  # drop genes with no NT reference

## 3) Order treatments
df_delta_ntmean_inter$Treatment <- factor(df_delta_ntmean_inter$Treatment,
                                          levels = c("CDK1", "Clk"))

## 4) Summary for error bars (mean ± SE of ΔPIR)
df_delta_summary_inter <- df_delta_ntmean_inter %>%
  group_by(Gene, Treatment) %>%
  summarise(
    delta_mean = mean(delta_PIR, na.rm = TRUE),
    delta_sd   = sd(delta_PIR, na.rm = TRUE),
    n          = dplyr::n(),
    delta_se   = delta_sd / sqrt(n),
    .groups = "drop"
  )

## 5) Plot: replicate points + overall mean ± SE (per gene), NT not shown
library(ggplot2)
library(patchwork)

# Function to make one plot per gene
plot_gene <- function(gene_name) {
  df_gene <- df_delta_ntmean_inter %>% filter(Gene == gene_name)
  df_gene_summary <- df_delta_summary_inter %>% filter(Gene == gene_name)
  
  ggplot(df_gene, aes(x = Treatment, y = delta_PIR, color = Treatment)) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_jitter(width = 0.15, size = 2, alpha = 0.85) +
    geom_point(data = df_gene_summary,
               aes(x = Treatment, y = delta_mean),
               shape = 18, size = 4, color = "black", inherit.aes = FALSE) +
    geom_errorbar(data = df_gene_summary,
                  aes(x = Treatment,
                      ymin = delta_mean - delta_se,
                      ymax = delta_mean + delta_se),
                  width = 0.2, color = "black", inherit.aes = FALSE) +
    scale_y_continuous(limits = c(-10, 100),
                       breaks = seq(-10, 100, by = 20)) +
    labs(x = "Treatment", y = "ΔPIR vs NT mean (%)", title = gene_name) +
    theme_classic() +
    theme(
      legend.position = "none",
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black", linewidth = 0.6),
      panel.grid = element_blank(),
      plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
    )
}

# Get list of genes
genes <- unique(df_delta_ntmean_inter$Gene)

# Generate plots
plots <- lapply(genes, plot_gene)

# Arrange into 3 rows x 2 columns
final_plot_interphase <- wrap_plots(plots, nrow = 3, ncol = 2)
final_plot_interphase

###### NUCLEAR CYTOPLASMIC ENRICHMENT ######



## 2) Reshape to long format
df_enrich <- tibble()

for (gene in c("CENPT","METTL3","TERT","Brd8","FANCA")) {
  for (treat in c("NT","CDK1","Clk_DYRK")) {
    nuc_col  <- paste0(gene, "_Inter_", treat, "_Exon_nuc")
    cyto_col <- paste0(gene, "_Inter_", treat, "_Exon_cyto")
    rep_col  <- paste0(gene, "_Inter_", treat, "_replica")  # <-- FIXED
    
    temp <- df_CLK_CDK_DYRK %>%
      select(all_of(c(nuc_col, cyto_col, rep_col))) %>%
      rename(
        nuc = !!sym(nuc_col),
        cyto = !!sym(cyto_col),
        replica = !!sym(rep_col)
      ) %>%
      mutate(
        Gene = gene,
        Treatment = gsub("Clk_DYRK", "Clk", treat),
        Enrichment = (nuc / (nuc + cyto)) * 100
      ) %>%
      filter(!is.na(Enrichment), !is.na(replica))
    
    df_enrich <- bind_rows(df_enrich, temp)
  }
}



## 3) Plot nuclear enrichment (% nuclear RNA)
ggplot(df_enrich, aes(x = Treatment, y = Enrichment, color = Treatment)) +
  geom_boxplot(alpha = 0.2, width = 0.5, outlier.shape = NA) +  # boxplots
  geom_jitter(width = 0.15, size = 2, alpha = 0.7) +            # replicate points
  facet_wrap(~ Gene, nrow = 3, ncol = 2) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 10)) +
  labs(
    x = "Treatment",
    y = "Nuclear enrichment (%)",
    title = "Interphase nuclear vs cytoplasmic exon enrichment"
  ) +
  theme_classic() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12, face = "bold"),
    legend.position = "none",
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.ticks = element_line(color = "black", linewidth = 0.6),
    panel.grid = element_blank()
  )

