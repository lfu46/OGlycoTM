# import packages
library(tidyverse)

# import results from psm.tsv file
OGlyco_Jurkat_raw <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Jurkat/EThcD_OPair_Search_1/OGlycoTM_Jurkat/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(!is.na(Total.Glycan.Composition)) |> 
  select(
    Spectrum, 
    Peptide, Peptide.Length, Charge, 
    Observed.Mass, Observed.M.Z, Delta.Mass, Hyperscore, Nextscore,
    Number.of.Enzymatic.Termini, Number.of.Missed.Cleavages, Protein.Start, Protein.End, 
    Assigned.Modifications, O.Pair.Score:Site.Probabilities,
    Ratio.138.144 = ..138.144.Ratio, 
    Has.N.Glyc.Sequon, Paired.Scan.Num, Parent.Scan.Number,
    Protein.ID:Protein.Description
  )

## data filtering
## exclude glycopeptides with localized sites on cysteine residues
## exclude glycopeptides that contain cysteine in their sequences but lack localized sites
OGlyco_Jurkat_localized <- OGlyco_Jurkat_raw |> 
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

OGlyco_Jurkat_nolocalized <- OGlyco_Jurkat_raw |> 
  filter(!Confidence.Level %in% c('Level1', 'Level1b')) |> 
  filter(!str_detect(Peptide, 'C'))

OGlyco_Jurkat_bonafide <- bind_rows(
  OGlyco_Jurkat_localized,
  OGlyco_Jurkat_nolocalized
)

write_csv(
  OGlyco_Jurkat_bonafide,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Jurkat/EThcD_OPair_Search_1/OGlyco_Jurkat_bonafide.csv'
)

# import filtered results
library(tidyverse)

OGlyco_Jurkat_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Jurkat/EThcD_OPair_Search_1/OGlyco_Jurkat_bonafide.csv'
)

# check preliminary results
# O-GlcNAc
Total_GlycoPSM_OGlcNAc <- OGlyco_Jurkat_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)', 'HexNAt(1)TMT6plex(1)')) |> 
  nrow()

Unique_Glycoprotein_OGlcNAc <- OGlyco_Jurkat_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)', 'HexNAt(1)TMT6plex(1)')) |> 
  distinct(Protein.ID) |> 
  nrow()

GlycoPSM_Glycoprotein_Results_OGlcNAc <- tibble(
  category = c('Total GlycoPSM', 'Unique Glycoprotein'),
  number = c(6634, 691)
)

# O-GalNAc
Total_GlycoPSM_OGalNAc <- OGlyco_Jurkat_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1)', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1)')) |> 
  nrow()

Unique_Glycoprotein_OGalNAc <- OGlyco_Jurkat_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1)', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1)')) |> 
  distinct(Protein.ID) |> 
  nrow()

GlycoPSM_Glycoprotein_Results_OGalNAc <- tibble(
  category = c('Total GlycoPSM', 'Unique Glycoprotein'),
  number = c(1233, 246)
)

# check the total glycoPSM group by raw file
# O-GlcNAc
Total_GlycoPSM_OGlcNAc_by_raw_file <- OGlyco_Jurkat_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)', 'HexNAt(1)TMT6plex(1)')) |> 
  mutate(
    raw_file = sub("\\..*", "", Spectrum)
  ) |> 
  group_by(raw_file) |> 
  count() |> 
  ungroup()

# O-GalNAc
Total_GlycoPSM_OGalNAc_by_raw_file <- OGlyco_Jurkat_bonafide |> 
  filter(Total.Glycan.Composition %in% c('HexNAt(1)GAO_Methoxylamine(1)', 'HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1)')) |> 
  mutate(
    raw_file = sub("\\..*", "", Spectrum)
  ) |> 
  group_by(raw_file) |> 
  count() |> 
  ungroup()

## plot for preliminary results
# Total glycoPSMs and unique glycoproteins
Total_GlycoPSMs_Unique_Glycoproteins <- bind_rows(
  GlycoPSM_Glycoprotein_Results_OGlcNAc %>% mutate(glycosylation = 'O-GlcNAc'),
  GlycoPSM_Glycoprotein_Results_OGalNAc %>% mutate(glycosylation = 'O-GalNAc')
)

Barplot_Total_GlycoPSMs_Unique_Glycoproteins <- Total_GlycoPSMs_Unique_Glycoproteins |> 
  mutate(glycosylation = factor(glycosylation, levels = c("O-GlcNAc", "O-GalNAc"))) |> 
  ggplot(aes(x = glycosylation, y = number, fill = category)) +
  geom_bar(stat = 'identity', position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = number), 
            position = position_dodge(width = 0.8), 
            vjust = -0.1,
            size = 2.5) +
  scale_fill_manual(values = c('Total GlycoPSM' = '#2C7BB6', 
                               'Unique Glycoprotein' = '#D7191C')) +
  labs(x = NULL,
       y = 'No. of O-glycan identifications',
       fill = NULL) +
  theme_bw() +
  theme(
    legend.position = 'right',
    axis.text = element_text(color = 'black', size = 7, family = 'arial'),
    axis.title = element_text(color = 'black', size = 7, family = 'arial'),
    legend.text = element_text(color = 'black', size = 5, family = 'arial'),
    legend.key.size = unit(0.1, 'in')
  )

ggsave(
  filename = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Jurkat/EThcD_OPair_Search_1/Barplot_Total_GlycoPSMs_Unique_Glycoproteins.tiff',
  plot = Barplot_Total_GlycoPSMs_Unique_Glycoproteins,
  width = 3, height = 2, units = 'in', dpi = 300
)

# Total glycoPSMs by raw file
# O-GlcNAc
barplot_Total_GlycoPSM_OGlcNAc_by_raw_file <- Total_GlycoPSM_OGlcNAc_by_raw_file |> 
  mutate(
    # Extract the fraction number from raw_file
    fraction = as.numeric(str_extract(raw_file, "(?<=_OG_)\\d+(?=_)")),
    # Create new label
    x_label = paste0("F", fraction)
  ) |> 
  # Arrange by fraction number to ensure correct order
  arrange(fraction) |> 
  # Convert x_label to factor with levels in the arranged order
  mutate(x_label = factor(x_label, levels = x_label)) |>
  ggplot(aes(x = x_label, y = n)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Fraction",
    y = "No. of total O-GlcNAc PSMs"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(color = 'black', size = 7, family = 'arial', angle = 30),
    axis.text.y = element_text(color = 'black', size = 7, family = 'arial'),
    axis.title = element_text(color = 'black', size = 7, family = 'arial')
  )

ggsave(
  filename = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Jurkat/EThcD_OPair_Search_1/barplot_Total_GlycoPSM_OGlcNAc_by_raw_file.tiff',
  plot = barplot_Total_GlycoPSM_OGlcNAc_by_raw_file,
  width = 3, height = 2.5, units = 'in', dpi = 300
)

# O-GalNAc
barplot_Total_GlycoPSM_OGalNAc_by_raw_file <- Total_GlycoPSM_OGalNAc_by_raw_file |> 
  mutate(
    # Extract the fraction number from raw_file
    fraction = as.numeric(str_extract(raw_file, "(?<=_OG_)\\d+(?=_)")),
    # Create new label
    x_label = paste0("F", fraction)
  ) |> 
  # Arrange by fraction number to ensure correct order
  arrange(fraction) |> 
  # Convert x_label to factor with levels in the arranged order
  mutate(x_label = factor(x_label, levels = x_label)) |>
  ggplot(aes(x = x_label, y = n)) +
  geom_bar(stat = "identity") +
  labs(
    x = "Fraction",
    y = "No. of total O-GalNAc PSMs"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(color = 'black', size = 7, family = 'arial', angle = 30),
    axis.text.y = element_text(color = 'black', size = 7, family = 'arial'),
    axis.title = element_text(color = 'black', size = 7, family = 'arial')
  )

ggsave(
  filename = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Jurkat/EThcD_OPair_Search_1/barplot_Total_GlycoPSM_OGalNAc_by_raw_file.tiff',
  plot = barplot_Total_GlycoPSM_OGalNAc_by_raw_file,
  width = 3, height = 2.5, units = 'in', dpi = 300
)
