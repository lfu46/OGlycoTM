# OGalNAc_reported_site_examples_plot.R
# Plot representative previously reported O-GalNAc site examples for Figure 6.
#
# Default panel uses the three most defensible examples from the current data:
#   1. TFRC T104 (Jurkat)
#   2. IGSF8 T169 (Jurkat)
#   3. CD99 T41 (HEK293T)
#
# Optional HepG2 examples are listed below. Add them to `candidate_sites`
# if visual balance across cell types is preferred.

library(tidyverse)

source_file_path <- "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/"
figure_file_path <- "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/"
fallback_output_dir <- "/Users/longpingfu/Downloads/OGlycoTM/data_analysis/Figure6_OGalNAc_preview"

colors_cell <- c(
  "HEK293T" = "#4DBBD5",
  "HepG2" = "#F39B7F",
  "Jurkat" = "#00A087"
)

default_output_dir <- paste0(figure_file_path, "Figure6_OGalNAc")
if (!dir.exists(default_output_dir)) {
  dir.create(default_output_dir, recursive = TRUE, showWarnings = FALSE)
}

output_dir <- if (file.access(default_output_dir, 2) == 0) {
  default_output_dir
} else {
  dir.create(fallback_output_dir, recursive = TRUE, showWarnings = FALSE)
  fallback_output_dir
}

candidate_sites <- tribble(
  ~CellType,  ~Protein.ID, ~Gene,   ~modified_residue, ~site_number, ~SiteLabel,
  "Jurkat",   "P02786",    "TFRC",  "T",               104,          "TFRC T104",
  "Jurkat",   "Q969P0",    "IGSF8", "T",               169,          "IGSF8 T169",
  "HEK293T",  "P14209",    "CD99",  "T",               41,           "CD99 T41"
)

candidate_sites_optional <- tribble(
  ~CellType, ~Protein.ID, ~Gene,   ~modified_residue, ~site_number, ~SiteLabel,
  "HepG2",   "P15529",    "CD46",  "T",               163,          "CD46 T163",
  "HepG2",   "P01130",    "LDLR",  "T",               108,          "LDLR T108"
)

load_site_norm <- function(cell_type) {
  read_csv(
    paste0(source_file_path, "normalization/OGalNAc_site_norm_", cell_type, ".csv"),
    show_col_types = FALSE
  ) |>
    mutate(CellType = cell_type)
}

site_norm_all <- bind_rows(
  load_site_norm("HEK293T"),
  load_site_norm("HepG2"),
  load_site_norm("Jurkat")
)

selected_sites <- site_norm_all |>
  inner_join(
    candidate_sites,
    by = c("CellType", "Protein.ID", "Gene", "modified_residue", "site_number")
  ) |>
  mutate(
    site_key = paste0(Protein.ID, "_", modified_residue, site_number)
  )

missing_sites <- anti_join(
  candidate_sites,
  distinct(
    selected_sites,
    CellType, Protein.ID, Gene, modified_residue, site_number
  ),
  by = c("CellType", "Protein.ID", "Gene", "modified_residue", "site_number")
)

if (nrow(missing_sites) > 0) {
  stop("Some candidate sites were not found:\n", paste(capture.output(print(missing_sites)), collapse = "\n"))
}

source_data <- selected_sites |>
  select(
    SiteLabel, CellType, Protein.ID, Gene, modified_residue, site_number,
    Intensity.Tuni_1_sl_tmm, Intensity.Tuni_2_sl_tmm, Intensity.Tuni_3_sl_tmm,
    Intensity.Ctrl_4_sl_tmm, Intensity.Ctrl_5_sl_tmm, Intensity.Ctrl_6_sl_tmm
  ) |>
  pivot_longer(
    cols = starts_with("Intensity."),
    names_to = "Channel",
    values_to = "NormIntensity"
  ) |>
  mutate(
    Condition = case_when(
      str_detect(Channel, "Tuni") ~ "Tuni",
      str_detect(Channel, "Ctrl") ~ "Ctrl",
      TRUE ~ NA_character_
    ),
    Replicate = case_when(
      str_detect(Channel, "_1_") ~ "1",
      str_detect(Channel, "_2_") ~ "2",
      str_detect(Channel, "_3_") ~ "3",
      str_detect(Channel, "_4_") ~ "4",
      str_detect(Channel, "_5_") ~ "5",
      str_detect(Channel, "_6_") ~ "6",
      TRUE ~ NA_character_
    ),
    Condition = factor(Condition, levels = c("Ctrl", "Tuni")),
    SiteLabel = factor(SiteLabel, levels = candidate_sites$SiteLabel),
    log2Intensity = log2(NormIntensity),
    PlotColor = if_else(Condition == "Ctrl", "grey65", colors_cell[CellType])
  )

summary_table <- source_data |>
  group_by(SiteLabel, CellType, Condition) |>
  summarize(
    mean_log2_intensity = mean(log2Intensity),
    mean_norm_intensity = mean(NormIntensity),
    .groups = "drop"
  ) |>
  pivot_wider(
    names_from = Condition,
    values_from = c(mean_log2_intensity, mean_norm_intensity)
  ) |>
  mutate(log2_Tuni_vs_Ctrl = mean_log2_intensity_Tuni - mean_log2_intensity_Ctrl)

write_csv(
  source_data,
  file.path(output_dir, "Figure6D_reported_sites_source_data.csv")
)
write_csv(
  summary_table,
  file.path(output_dir, "Figure6D_reported_sites_summary.csv")
)

Figure6D <- source_data |>
  ggplot(aes(x = Condition, y = log2Intensity)) +
  geom_boxplot(
    width = 0.45,
    outlier.shape = NA,
    fill = "white",
    color = "black",
    linewidth = 0.3
  ) +
  geom_point(
    aes(color = PlotColor),
    position = position_jitter(width = 0.08, height = 0),
    size = 2
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    shape = 95,
    size = 6,
    color = "black"
  ) +
  facet_wrap(~ SiteLabel, scales = "free_y", nrow = 1) +
  scale_color_identity() +
  labs(
    x = "",
    y = expression(log[2] * "(normalized site abundance)")
  ) +
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    axis.title = element_text(size = 9, color = "black"),
    axis.text = element_text(size = 9, color = "black"),
    strip.text = element_text(size = 9, face = "bold"),
    legend.position = "none"
  )

ggsave(
  file.path(output_dir, "Figure6D_reported_sites_examples.pdf"),
  Figure6D,
  width = 4.6,
  height = 1.9,
  units = "in"
)

cat("Figure 6D example panel saved.\n")
cat("Outputs written to:", output_dir, "\n")
