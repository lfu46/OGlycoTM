# Figure 1: Workflow and Example TMT Reporter Ion Barplot

library(tidyverse)

# Source data paths and colors
source('data_source.R')

# =============================================================================
# Figure 1 Panel: TMT Reporter Ion Barplot (Example Spectrum)
# =============================================================================
# Spectrum: Eclipse_LF_OGlycoTM_HEK293T_OG_7_10062025_ppm.26380.26380.3
# Peptide: SSNDYTSQMYSAK
# Gene: EYA3
# Site: S63 (Serine 63)
# Site Index: Q99504_S63

# Load normalized site data
OGlcNAc_site_norm_HEK293T <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_site_norm_HEK293T.csv'),
  show_col_types = FALSE
)

# Extract EYA3 S63 site data
eya3_s63 <- OGlcNAc_site_norm_HEK293T |>
  filter(site_index == "Q99504_S63")

# Extract SL+TMM normalized intensities for TMT channels
tmt_intensities <- tibble(
  Channel = c("126", "127", "128", "129", "130", "131"),
  Condition = c("Tuni", "Tuni", "Tuni", "Ctrl", "Ctrl", "Ctrl"),
  Intensity = c(
    eya3_s63$Intensity.Tuni_1_sl_tmm,
    eya3_s63$Intensity.Tuni_2_sl_tmm,
    eya3_s63$Intensity.Tuni_3_sl_tmm,
    eya3_s63$Intensity.Ctrl_4_sl_tmm,
    eya3_s63$Intensity.Ctrl_5_sl_tmm,
    eya3_s63$Intensity.Ctrl_6_sl_tmm
  )
) |>
  mutate(
    # Calculate relative abundance (normalized to maximum = 100%)
    Relative_Abundance = Intensity / max(Intensity) * 100,
    # Keep channel order
    Channel = factor(Channel, levels = c("126", "127", "128", "129", "130", "131"))
  )

# Print intensity values
cat("=== EYA3 S63 TMT Intensities (SL+TMM Normalized) ===\n")
print(tmt_intensities)

# Define colors for each TMT channel (matching reference image style)
tmt_colors <- c(
  "126" = "#1F77B4",  # Blue
  "127" = "#FF7F0E",  # Orange
  "128" = "#2CA02C",  # Green
  "129" = "#D62728",  # Red
  "130" = "#9467BD",  # Purple
  "131" = "#8C564B"   # Brown
)

# Create barplot
tmt_barplot <- ggplot(tmt_intensities, aes(x = Channel, y = Relative_Abundance, fill = Channel)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = tmt_colors) +
  scale_y_continuous(
    expand = expansion(mult = c(0, 0.05)),
    limits = c(0, 105)
  ) +
  labs(
    x = NULL,
    y = "Relative Abundance"
  ) +
  theme_classic(base_size = 6) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 6),
    axis.text.y = element_text(size = 6),
    axis.title.y = element_text(size = 6),
    legend.position = "none",
    axis.line = element_line(linewidth = 0.3),
    axis.ticks = element_line(linewidth = 0.3),
    plot.margin = margin(5, 5, 5, 5)
  )

# Display plot
print(tmt_barplot)

# Save plot
ggsave(
  filename = paste0(figure_file_path, "Figure1/TMT_reporter_barplot_EYA3_S63.pdf"),
  plot = tmt_barplot,
  width = 2,
  height = 1,
  units = "in"
)

cat("\nPlot saved to:", figure_file_path, "Figure1/TMT_reporter_barplot_EYA3_S63.pdf\n")
