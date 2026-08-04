# NGlcNAc_score_comparison_plot.R
# Violin + boxplot comparing Hyperscore and DeltaScore between O-GlcNAc and N-GlcNAc
# For the same spectra assigned differently in two searches

library(tidyverse)

# Load the statistical comparison data
stats <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HEK293T/OGlyco/NGlcNAc_comparison/statistical_comparison.csv',
  show_col_types = FALSE
)

figure_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/'
dir.create(paste0(figure_file_path, "NGlcNAc_comparison"), showWarnings = FALSE)

# Need to reconstruct paired data from the PSM files
# Load both searches
original <- read_tsv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HEK293T/OGlyco/EThcD_OPair_TMT_Search/OGlyco/psm.tsv',
  show_col_types = FALSE
)
combined <- read_tsv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/OGlycoTM_HEK293T/OGlyco/EThcD_OPair_NGlcNAc_TMT_Search/OGlyco/psm.tsv',
  show_col_types = FALSE
)

# Parse scan keys
parse_scan <- function(spectrum) {
  parts <- str_split(spectrum, "\\.", simplify = TRUE)
  scan <- as.integer(parts[, ncol(parts) - 2])
  base <- parts[, 1]
  paste0(scan, "_", base)
}

original <- original |> mutate(scan_key = parse_scan(Spectrum))
combined <- combined |> mutate(scan_key = parse_scan(Spectrum))

# Load bonafide scans
bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/filtered/OGlyco_HEK293T_bonafide.csv',
  show_col_types = FALSE
)
oglcnac_bf <- bonafide |>
  filter(str_detect(Assigned.Modifications, "528\\.2859|299\\.1230")) |>
  filter(!str_detect(Assigned.Modifications, "555\\.2968|326\\.1339")) |>
  mutate(scan_key = parse_scan(Spectrum))
bf_keys <- unique(oglcnac_bf$scan_key)

# N-GlcNAc from combined, in bonafide scans
n_glcnac <- combined |>
  filter(str_detect(`Assigned Modifications`, "N\\(528\\.2859\\)|N\\(299\\.1230\\)")) |>
  filter(scan_key %in% bf_keys)

# Merge
merged <- n_glcnac |>
  select(scan_key, Hyperscore_n = Hyperscore, Nextscore_n = Nextscore, Peptide_n = Peptide) |>
  inner_join(
    original |> select(scan_key, Hyperscore_o = Hyperscore, Nextscore_o = Nextscore, Peptide_o = Peptide),
    by = "scan_key"
  ) |>
  filter(Peptide_n == Peptide_o) |>
  mutate(
    DeltaScore_o = Hyperscore_o - Nextscore_o,
    DeltaScore_n = Hyperscore_n - Nextscore_n
  )

cat("Paired PSMs:", nrow(merged), "\n")

# Reshape for plotting
hyper_long <- merged |>
  select(scan_key, `O-GlcNAc` = Hyperscore_o, `N-GlcNAc` = Hyperscore_n) |>
  pivot_longer(cols = c(`O-GlcNAc`, `N-GlcNAc`), names_to = "Assignment", values_to = "Hyperscore") |>
  mutate(Assignment = factor(Assignment, levels = c("O-GlcNAc", "N-GlcNAc")))

delta_long <- merged |>
  select(scan_key, `O-GlcNAc` = DeltaScore_o, `N-GlcNAc` = DeltaScore_n) |>
  pivot_longer(cols = c(`O-GlcNAc`, `N-GlcNAc`), names_to = "Assignment", values_to = "DeltaScore") |>
  mutate(Assignment = factor(Assignment, levels = c("O-GlcNAc", "N-GlcNAc")))

# Colors
colors_assign <- c("O-GlcNAc" = "#F39B7F", "N-GlcNAc" = "#8491B4")

library(ggpubr)
library(rstatix)

# Wilcoxon test for Hyperscore
hyper_stat <- hyper_long |>
  wilcox_test(Hyperscore ~ Assignment, paired = TRUE) |>
  add_significance() |>
  add_xy_position(x = "Assignment")

# Wilcoxon test for DeltaScore
delta_stat <- delta_long |>
  wilcox_test(DeltaScore ~ Assignment, paired = TRUE) |>
  add_significance() |>
  add_xy_position(x = "Assignment")

# Hyperscore plot
p_hyper <- hyper_long |>
  ggplot(aes(x = Assignment, y = Hyperscore, fill = Assignment)) +
  geom_violin(alpha = 0.5, width = 0.8, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, linewidth = 0.3, alpha = 0.8) +
  stat_pvalue_manual(
    hyper_stat, label = "p.signif",
    y.position = max(hyper_long$Hyperscore, na.rm = TRUE) * 1.05,
    size = 3, bracket.size = 0.3, tip.length = 0.01
  ) +
  scale_fill_manual(values = colors_assign) +
  labs(x = "", y = "Hyperscore") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave(
  paste0(figure_file_path, "NGlcNAc_comparison/Hyperscore_violin.pdf"),
  p_hyper, width = 2, height = 2.2, units = "in"
)

# DeltaScore plot
p_delta <- delta_long |>
  ggplot(aes(x = Assignment, y = DeltaScore, fill = Assignment)) +
  geom_violin(alpha = 0.5, width = 0.8, linewidth = 0.3) +
  geom_boxplot(width = 0.15, outlier.size = 0.3, linewidth = 0.3, alpha = 0.8) +
  stat_pvalue_manual(
    delta_stat, label = "p.signif",
    y.position = max(delta_long$DeltaScore, na.rm = TRUE) * 1.05,
    size = 3, bracket.size = 0.3, tip.length = 0.01
  ) +
  scale_fill_manual(values = colors_assign) +
  labs(x = "", y = "DeltaScore") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 9, color = "black"),
    axis.text.x = element_text(size = 9, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 9, color = "black"),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave(
  paste0(figure_file_path, "NGlcNAc_comparison/DeltaScore_violin.pdf"),
  p_delta, width = 2, height = 2.2, units = "in"
)

# Combined plot
library(patchwork)

p_combined <- p_hyper + p_delta +
  plot_layout(ncol = 2)

ggsave(
  paste0(figure_file_path, "NGlcNAc_comparison/Score_comparison_combined.pdf"),
  p_combined, width = 3.5, height = 2.2, units = "in"
)

cat("Hyperscore median: O-GlcNAc =", median(merged$Hyperscore_o), ", N-GlcNAc =", median(merged$Hyperscore_n), "\n")
cat("DeltaScore median: O-GlcNAc =", median(merged$DeltaScore_o), ", N-GlcNAc =", median(merged$DeltaScore_n), "\n")
cat("Plots saved.\n")
