# import packages
library(tidyverse)

# import HEK293T Preliminary Results
Preliminary_Results <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlyco_HEK293T_Preliminary_Results/OGlycoTM_EThcD_OPair_TMT_Search_1/OGlyco_HEK293T_Preliminary_Results/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(!is.na(O.Pair.Score)) |> 
  mutate(
    PPM = Delta.Mass/Calculated.Peptide.Mass
  ) |>
  select(
    Spectrum,
    Peptide, Peptide.Length, Charge, 
    Observed.Mass, Calculated.Peptide.Mass, Observed.M.Z, Delta.Mass, PPM, Hyperscore, Nextscore,
    Number.of.Enzymatic.Termini, Number.of.Missed.Cleavages, Protein.Start, Protein.End, 
    Assigned.Modifications, O.Pair.Score:Site.Probabilities, 
    Ratio.138.144 = ..138.144.Ratio,
    Has.N.Glyc.Sequon, Paired.Scan.Num, Parent.Scan.Number,
    Protein.ID:Protein.Description
  )

# data filtering
# exclude glycopeptides with localized sites on cysteine resiudes
# exclude glycopeptides that contains cysteine in their sequences but lack localized sites
Preliminary_Results_localized <- Preliminary_Results |> 
  filter(Confidence.Level %in% c('Level1', 'Level1b')) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(555\\.2968\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(528\\.2859\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(326\\.1339\\)')
  ) |> 
  filter(
    !str_detect(Assigned.Modifications, 'C\\(299\\.1230\\)')
  )

Preliminary_Results_nolocalized <- Preliminary_Results |> 
  filter(!Confidence.Level %in% c('Level1', 'Level1b')) |> 
  filter(!str_detect(Peptide, 'C'))

Preliminary_Results_bonafide <- bind_rows(
  Preliminary_Results_localized,
  Preliminary_Results_nolocalized
)

write_csv(
  Preliminary_Results_bonafide,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlyco_HEK293T_Preliminary_Results/Preliminary_Results_bonafide.csv'
)

# import filtered results
library(tidyverse)

OGlyco_HEK293T_bonafide_reassigned <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/data_source/OGlyco_HEK293T_bonafide_reassigned.csv'
)

# check number of total glycoPSM and unique glycoprotein
# O-GlcNAc
Total_GlycoPSM_OGlcNAc <- OGlyco_HEK293T_bonafide_reassigned |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)', 'HexNAt(1)TMT6plex(1)')) |> 
  nrow()

Unique_Glycoprotein_OGlcNAc <- OGlyco_HEK293T_bonafide_reassigned |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)', 'HexNAt(1)TMT6plex(1)')) |> 
  distinct(Protein.ID) |> 
  nrow()

GlycoPSM_Glycoprotein_Results_OGlcNAc <- tibble(
  category = c('Total GlycoPSM', 'Unique Glycoprotein'),
  number = c(9242, 738)
)

# O-GalNAc
Total_GlycoPSM_OGalNAc <- OGlyco_HEK293T_bonafide_reassigned |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1)', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1)')) |> 
  nrow()

Unique_Glycoprotein_OGalNAc <- OGlyco_HEK293T_bonafide_reassigned |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1)', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1)')) |> 
  distinct(Protein.ID) |> 
  nrow()

GlycoPSM_Glycoprotein_Results_OGalNAc <- tibble(
  category = c('Total GlycoPSM', 'Unique Glycoprotein'),
  number = c(1143, 235)
)

# check number of localized and non-localized glycoPSM
# O-GlcNAc
Localized_GlycoPSM_OGlcNAc <- Preliminary_Results_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)', 'HexNAt(1)TMT6plex(1)')) |> 
  filter(Confidence.Level %in% c('Level1', 'Level1b')) |> 
  nrow()

Nonlocalized_GlycoPSM_OGlcNAc <- Preliminary_Results_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)', 'HexNAt(1)TMT6plex(1)')) |> 
  filter(!Confidence.Level %in% c('Level1', 'Level1b')) |> 
  nrow()

Localized_Nonlocalized_Results_OGlcNAc <- tibble(
  category = c('Localized', 'Nonlocalized'),
  number = c(383, 569)
)

# O-GalNAc
Localized_GlycoPSM_OGalNAc <- Preliminary_Results_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1)', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1)')) |> 
  filter(Confidence.Level %in% c('Level1', 'Level1b')) |> 
  nrow()

Nonlocalized_GlycoPSM_OGalNAc <- Preliminary_Results_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1)', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1)')) |> 
  filter(!Confidence.Level %in% c('Level1', 'Level1b')) |> 
  nrow()

Localized_Nonlocalized_Results_OGalNAc <- tibble(
  category = c('Localized', 'Nonlocalized'),
  number = c(320, 213)
)

# GO Enrichment Analysis
library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

Total_Glycoprotein <- Preliminary_Results_bonafide |> 
  dplyr::distinct(Protein.ID) |> 
  pull()

OGlcNAc_GLycoprotein <- Preliminary_Results_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)', 'HexNAt(1)TMT6plex(1)')) |> 
  dplyr::distinct(Protein.ID) |> 
  pull()

OGalNAc_GLycoprotein <- Preliminary_Results_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1)', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1)')) |> 
  dplyr::distinct(Protein.ID) |> 
  pull()

OGlcNAc_GO <- enrichGO(
  gene = OGlcNAc_GLycoprotein,
  OrgDb = org.Hs.eg.db,
  universe = Total_Glycoprotein,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OGlcNAc_GO@result, file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlyco_HEK293T_Preliminary_Results/OGlcNAc_GO.csv'
)

OGalNAc_GO <- enrichGO(
  gene = OGalNAc_GLycoprotein,
  OrgDb = org.Hs.eg.db,
  universe = Total_Glycoprotein,
  keyType = 'UNIPROT',
  ont = 'ALL',
  pvalueCutoff = 1,
  qvalueCutoff = 1
)

write_csv(
  OGalNAc_GO@result, file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlyco_HEK293T_Preliminary_Results/OGalNAc_GO.csv'
)

# barplot
library(tidyverse)

GlycoPSM_Glycoprotein_combined <- bind_rows(
  GlycoPSM_Glycoprotein_Results_OGlcNAc %>% mutate(modification = "O-GlcNAc"),
  GlycoPSM_Glycoprotein_Results_OGalNAc %>% mutate(modification = "O-GalNAc")
)

GlycoPSM_Evidence_combined <- bind_rows(
  Localized_Nonlocalized_Results_OGlcNAc %>% mutate(modification = "O-GlcNAc"),
  Localized_Nonlocalized_Results_OGalNAc %>% mutate(modification = "O-GalNAc")
)

# Number of total glycoPSM and unique glycoprotein
GlycoPSM_Glycoprotein_barplot <- GlycoPSM_Glycoprotein_combined |> 
  ggplot(aes(x = category, y = number, fill = factor(modification, c('O-GlcNAc', 'O-GalNAc')))) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("O-GlcNAc" = "#1f77b4", "O-GalNAc" = "#ff7f0e")) +
  ylim(0, 1000) +
  labs(
    x = NULL,
    y = "Count",
    fill = ""
  ) +
  geom_text(
    aes(label = number), color = 'black', size = 3, position = position_dodge(width = 0.7), vjust = -0.1
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 8, color = 'black', angle = 30, hjust = 1),
    axis.text.y = element_text(size = 8, color = 'black'),
    axis.title.y = element_text(size = 8, color = 'black'),
    legend.position = "none"
  )

ggsave(
  filename = '/Volumes/cos-lab-rwu60/Longping/Retreat2025/Figures/GlycoPSM_Glycoprotein_barplot.tif',
  plot = GlycoPSM_Glycoprotein_barplot,
  height = 2, width = 2, units = 'in', dpi = 300
)

# Number of glycoPSM with different levels of evidence
GlycoPSM_Evidence_barplot <- GlycoPSM_Evidence_combined |> 
  ggplot(aes(x = category, y = number, fill = factor(modification, c('O-GlcNAc', 'O-GalNAc')))) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = c("O-GlcNAc" = "#1f77b4", "O-GalNAc" = "#ff7f0e")) +
  ylim(0, 600) +
  labs(
    x = NULL,
    y = "# of glycoPSM",
    fill = ""
  ) +
  geom_text(
    aes(label = number), color = 'black', size = 3, position = position_dodge(width = 0.7), vjust = -0.1
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(size = 8, color = 'black', angle = 30, hjust = 1),
    axis.text.y = element_text(size = 8, color = 'black'),
    axis.title.y = element_text(size = 8, color = 'black'),
    legend.position = "none"
  )

ggsave(
  filename = '/Volumes/cos-lab-rwu60/Longping/Retreat2025/Figures/GlycoPSM_Evidence_barplot.tif',
  plot = GlycoPSM_Evidence_barplot,
  height = 2, width = 2, units = 'in', dpi = 300
)

# GO enrichment analysis of O-GalNAc
library(tidyverse)

OGalNAc_GO <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlyco_HEK293T_Preliminary_Results/OGalNAc_GO.csv'
)

barplot_GO_OGalNAc <- OGalNAc_GO |> 
  filter(Description %in% c(
    'endoplasmic reticulum membrane',
    'protein folding',
    'protein maturation',
    'lysosome',
    'establishment of protein localization',
    'extracellular organelle'
  )) |> 
  ggplot() +
  geom_bar(
    aes(
      x = fct_reorder(Description, -log10(pvalue)), 
      y = -log10(pvalue)
    ),
    fill = '#ff7f0e', color = 'transparent', stat = 'identity'
  ) +
  labs(x = '', y = '') +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = '/Volumes/cos-lab-rwu60/Longping/Retreat2025/Figures/barplot_GO_OGalNAc.tif',
  plot = barplot_GO_OGalNAc,
  height = 2, width =  3, units = 'in', dpi = 300
)

# GO enrichment analysis of O-GlcNAc
library(tidyverse)

OGlcNAc_GO <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlyco_HEK293T_Preliminary_Results/OGlcNAc_GO.csv'
)

barplot_GO_OGlcNAc <- OGlcNAc_GO |> 
  filter(Description %in% c(
    'regulation of transcription by RNA polymerase II',
    'positive regulation of metabolic process',
    'regulation of DNA-templated transcription',
    'nucleus',
    'nucleoplasm',
    'chromatin'
  )) |> 
  ggplot() +
  geom_bar(
    aes(
      x = fct_reorder(Description, -log10(pvalue)), 
      y = -log10(pvalue)
    ),
    fill = '#1f77b4', color = 'transparent', stat = 'identity'
  ) +
  labs(x = '', y = '') +
  coord_flip() +
  theme_classic() +
  theme(
    axis.text = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = '/Volumes/cos-lab-rwu60/Longping/Retreat2025/Figures/barplot_GO_OGlcNAc.tif',
  plot = barplot_GO_OGlcNAc,
  height = 2, width =  3.5, units = 'in', dpi = 300
)

# O-GlcNAc on tyrosine residue
library(tidyverse)

Preliminary_Results_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlyco_HEK293T_Preliminary_Results/Preliminary_Results_bonafide.csv'
)

OGlcNAc_Tyrosine_Level_1 <- Preliminary_Results_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)', 'HexNAt(1)TMT6plex(1)')) |> 
  filter(str_detect(Assigned.Modifications, 'Y')) |> 
  filter(Confidence.Level == 'Level1') |> 
  mutate(site_number_insequence = as.numeric(str_extract(Site.Probabilities, '(?<=\\[)\\d+'))) |> 
  mutate(site_number_inprotein = Protein.Start + site_number_insequence -1) |> 
  mutate(site_index = paste(Protein.ID, site_number_inprotein, sep = "_"))

write_csv(
  OGlcNAc_Tyrosine_Level_1,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlyco_HEK293T_Preliminary_Results/OGlcNAc_Tyrosine_Level_1.csv'
)

OGlcNAc_Tyrosine_site_index <- OGlcNAc_Tyrosine_Level_1 |> 
  distinct(site_index) |> 
  pull()

# import O-GlcNAc on tyrosine residue from literature
# https://www.pnas.org/doi/10.1073/pnas.2409501121, Dataset S03 (XLSX)
library(readxl)
library(tidyverse)

OGlcNAc_Tyrosine_site_literature <- read_xlsx(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlyco_HEK293T_Preliminary_Results/pnas.2409501121.sd03.xlsx',
  range = "A5:D2835",
  col_names = c('UniProt_Accession', 'Protein.Description', 'Site.Position', 'Residue')
) |> 
  filter(Residue == 'Y') |> 
  mutate(site_index = paste(UniProt_Accession, Site.Position, sep = '_'))

write_csv(
  OGlcNAc_Tyrosine_site_literature,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlyco_HEK293T_Preliminary_Results/OGlcNAc_Tyrosine_site_literature.csv'
)

OGlcNAc_Tyrosine_site_literature_index <- OGlcNAc_Tyrosine_site_literature |> 
  distinct(site_index) |> 
  pull()

# venn diagram
library(VennDiagram)

venn.diagram(
  x = list(
    "My Data" = OGlcNAc_Tyrosine_site_index,
    "Literature" = OGlcNAc_Tyrosine_site_literature_index
  ),
  filename = "/Volumes/cos-lab-rwu60/Longping/Retreat2025/Figures/OGlcNAc_Tyrosine_Venn.tiff",
  imagetype = "tiff",
  resolution = 300,
  cex = 4,
  cat.cex = 4,
  cat.pos = c(0, 180),
  cat.dist = c(0.01, 0.01),
  
  lty = "blank",
  fill = c("blue", "red"),
  alpha = 0.3,
  scaled = TRUE
)
