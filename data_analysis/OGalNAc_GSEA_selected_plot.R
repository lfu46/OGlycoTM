# OGalNAc_GSEA_selected_plot.R
# Figure 6C (revised): GSEA dotplot with 4 negative enriched terms (Jurkat exclusive)

library(tidyverse)
library(devEMF)

source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'
figure_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/'

# Load all Jurkat GSEA results
gsea_all <- bind_rows(
  read_csv(paste0(source_file_path, "enrichment/OGalNAc_exclusive_GSEA_BP_Jurkat.csv"), show_col_types = FALSE) |> mutate(ont = "BP"),
  read_csv(paste0(source_file_path, "enrichment/OGalNAc_exclusive_GSEA_CC_Jurkat.csv"), show_col_types = FALSE) |> mutate(ont = "CC"),
  read_csv(paste0(source_file_path, "enrichment/OGalNAc_exclusive_GSEA_MF_Jurkat.csv"), show_col_types = FALSE) |> mutate(ont = "MF")
)

# Select 4 downregulated terms only (upregulated removed per advisor revision)
selected_descriptions <- c(
  "endoplasmic reticulum membrane",
  "endomembrane system",
  "cell periphery",
  "carbohydrate derivative metabolic process"
)

# Display labels with abbreviations
display_labels <- c(
  "endoplasmic reticulum membrane" = "ER membrane",
  "endomembrane system"            = "Endomembrane system",
  "cell periphery"                 = "Cell periphery",
  "carbohydrate derivative metabolic process" = "Carbohydrate deriv.\nmetabolism"
)

selected <- gsea_all |>
  filter(Description %in% selected_descriptions) |>
  mutate(
    setSize = as.numeric(setSize),
    label = display_labels[Description],
    label = factor(label, levels = rev(unname(display_labels)))
  )

cat("Selected terms:\n")
selected |> dplyr::select(ont, Description, NES, p.adjust, setSize) |> print(n = 10)

# Jurkat green gradient
Figure6C <- selected |>
  ggplot(aes(x = NES, y = label, size = setSize, color = p.adjust)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_color_gradient(low = "#004D40", high = "#80CBC4",
                       name = expression(adj.~italic(p)~value)) +
  scale_size_continuous(range = c(2, 6), name = "Set size", breaks = c(10, 20, 40), limits = c(5, 45)) +
  guides(size = guide_legend(override.aes = list(shape = 1, color = "black"))) +
  scale_x_continuous(breaks = c(-2, -1, 0)) +
  labs(
    x = "NES",
    y = ""
  ) +
  theme_bw(base_size = 6) +
  theme(
    axis.title = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.text.x = element_text(size = 7, color = "black"),
    legend.position = "right",
    legend.title = element_text(size = 6),
    legend.text = element_text(size = 6),
    legend.key.size = unit(0.15, "cm"),
    legend.spacing.y = unit(0.05, "cm"),
    legend.margin = margin(0, 0, 0, 0),
    legend.box.margin = margin(0, 0, 0, -3),
    panel.grid.minor = element_blank(),
    plot.margin = margin(2, 2, 2, 2)
  )

# Save PDF
ggsave(
  paste0(figure_file_path, "Figure6_OGalNAc/Figure6C_GSEA_Jurkat.pdf"),
  Figure6C, width = 2.2, height = 1.5, units = "in"
)
cat("Figure 6C PDF saved.\n")

# Save EMF
emf(paste0(figure_file_path, "Figure6_OGalNAc/Figure6C_GSEA_Jurkat.emf"),
    width = 2.2, height = 1.5)
print(Figure6C)
dev.off()
cat("Figure 6C EMF saved.\n")
