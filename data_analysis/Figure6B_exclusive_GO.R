# Figure6B_exclusive_GO.R
# Figure 6B: GO enrichment dotplot for cell-type exclusive O-GalNAc proteins

library(tidyverse)
library(devEMF)

source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'
figure_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/Figure6_OGalNAc/'

colors_cell <- c("HepG2" = "#F39B7F", "Jurkat" = "#00A087")

# Load cell-type exclusive GO results
go_data <- bind_rows(
  read_csv(paste0(source_file_path, "enrichment/OGalNAc_HepG2_exclusive_GO.csv"),
           show_col_types = FALSE) |> mutate(cell = "HepG2"),
  read_csv(paste0(source_file_path, "enrichment/OGalNAc_Jurkat_exclusive_GO.csv"),
           show_col_types = FALSE) |> mutate(cell = "Jurkat")
)

# Select terms shown in Figure 6B
selected_terms <- c(
  "extracellular space",
  "endoplasmic reticulum lumen",
  "glycoprotein metabolic process",
  "carbohydrate derivative biosynthetic process",
  "Golgi membrane"
)

# Display labels with abbreviations
display_labels <- c(
  "extracellular space"                       = "Extracellular space",
  "endoplasmic reticulum lumen"               = "ER lumen",
  "glycoprotein metabolic process"            = "Glycoprotein\nmetabolism",
  "carbohydrate derivative biosynthetic process" = "Carbohydrate derivative\nbiosynthesis",
  "Golgi membrane"                            = "Golgi membrane"
)

plot_data <- go_data |>
  filter(Description %in% selected_terms) |>
  mutate(
    neg_log10p = -log10(pvalue),
    label = display_labels[Description],
    label = factor(label, levels = rev(unname(display_labels))),
    cell = factor(cell, levels = c("HepG2", "Jurkat"))
  ) |>
  # Only keep dots with meaningful enrichment
  filter(neg_log10p > 1)

cat("Plot data:\n")
plot_data |> select(cell, Description, Count, neg_log10p) |> print(n = 20)

Figure6B <- ggplot(plot_data,
                   aes(x = neg_log10p, y = label, size = Count, color = cell)) +
  geom_point() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.3) +
  scale_color_manual(values = colors_cell, name = NULL) +
  scale_size_continuous(range = c(2, 6), name = "Count", breaks = c(10, 20, 30), limits = c(5, 45)) +
  guides(size = guide_legend(override.aes = list(shape = 1, color = "black"))) +
  labs(
    x = expression(-log[10]~"("*italic(p)~value*")"),
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

ggsave(paste0(figure_file_path, "Figure6B_exclusive_GO_combined.pdf"),
       Figure6B, width = 2.2, height = 1.5, units = "in")
cat("Figure 6B PDF saved.\n")

emf(paste0(figure_file_path, "Figure6B_exclusive_GO_combined.emf"),
    width = 2.2, height = 1.5)
print(Figure6B)
dev.off()
cat("Figure 6B EMF saved.\n")
