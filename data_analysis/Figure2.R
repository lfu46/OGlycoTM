# import packages
library(tidyverse)

# Figure 2. Identification of O-GlcNAc and O-GalNAc sites and proteins in HEK293T, HepG2 and Jurkat cells

# Figure 2A ---------------------------------------------------------------
# Indentification of O-GlcNAc and O-GalNAc proteins in three type of cells

# Number of unique O-GlcNAcylated protein
# HEK293T
OGlcNAc_protein_unique_number_HEK293T <- OGlyco_HEK293T_bonafide |> 
  filter(
    Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859')
  ) |> 
  distinct(Protein.ID) |> 
  nrow()

# HepG2
OGlcNAc_protein_unique_number_HepG2 <- OGlyco_HepG2_bonafide |> 
  filter(
    Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859')
  ) |> 
  distinct(Protein.ID) |> 
  nrow()

# Jurkat
OGlcNAc_protein_unique_number_Jurkat <- OGlyco_Jurkat_bonafide |> 
  filter(
    Total.Glycan.Composition %in% c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859')
  ) |> 
  distinct(Protein.ID) |> 
  nrow()

# Number of unique O-GalNAcylated protein
# HEK293T
OGalNAc_protein_unique_number_HEK293T <- OGlyco_HEK293T_bonafide |> 
  filter(
    Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1) % 326.1339', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968')
  ) |> 
  distinct(Protein.ID) |> 
  nrow()

# HepG2
OGalNAc_protein_unique_number_HepG2 <- OGlyco_HepG2_bonafide |> 
  filter(
    Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1) % 326.1339', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968')
  ) |> 
  distinct(Protein.ID) |> 
  nrow()

# Jurkat
OGalNAc_protein_unique_number_Jurkat <- OGlyco_Jurkat_bonafide |> 
  filter(
    Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1) % 326.1339', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968')
  ) |> 
  distinct(Protein.ID) |> 
  nrow()

# summary of the identification results
Figure2A_df <- data.frame(
  cell_type = rep(c("HEK293T", "HepG2", "Jurkat"), 2),
  glycan_type = factor(
    c(rep("O-GlcNAc", 3), rep("O-GalNAc", 3)),
    levels = c("O-GlcNAc", "O-GalNAc")
  ),
  protein_count = c(
    OGlcNAc_protein_unique_number_HEK293T,
    OGlcNAc_protein_unique_number_HepG2,
    OGlcNAc_protein_unique_number_Jurkat,
    OGalNAc_protein_unique_number_HEK293T,
    OGalNAc_protein_unique_number_HepG2,
    OGalNAc_protein_unique_number_Jurkat
  )
)

# barplot
Figure2A <- ggplot(Figure2A_df, aes(x = cell_type, y = protein_count, fill = glycan_type)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(
    aes(label = protein_count),
    position = position_dodge(width = 0.9),
    vjust = -0.2,
    size = 2,
    color = "black"
  ) +
  scale_fill_manual(values = colors_glycan) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  labs(
    x = NULL,
    y = "No. of unique glycoproteins",
    fill = NULL
  ) +
  theme_classic() +
  theme(
    title = element_text(family = "Helvetica", color = "black", size = 6),
    text = element_text(family = "Helvetica", color = "black"),
    axis.text = element_text(color = "black", size = 6),
    axis.text.x = element_text(angle = 30, hjust = 1),
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size = 6),
    legend.title = element_text(size = 6)
  )

ggsave(
  filename = paste0(figure_file_path, "Figure2/Figure2A.pdf"),
  plot = Figure2A,
  width = 2.5, height = 1.5, units = "in"
)

# Figure 2B ---------------------------------------------------------------
# 

