# import packages
library(tidyverse)

# import EThcD OPair results
EThcD_OPair_Results <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlyco_HEK293T/OGlycoTM_EThcD_OPair_TMT_Search_1/Eclipse_LF_OGlycoTM_HEK293T_OG_4_10062025_HCD_pd_EThcD_Test_ppm/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Entry.Name, 'HUMAN')) |> 
  filter(!is.na(O.Pair.Score)) |> 
  select(
    Spectrum,
    Peptide, Peptide.Length, Charge, Observed.Mass, Observed.M.Z, Delta.Mass, Hyperscore, Nextscore,
    Number.of.Enzymatic.Termini, Number.of.Missed.Cleavages, Protein.Start, Protein.End, 
    Assigned.Modifications, O.Pair.Score:Paired.Scan.Num,
    Protein.ID:Protein.Description
  )

# data filtering
# exclude glycopeptides with localized sites on cysteine resiudes
# exclude glycopeptides that contains cysteine in their sequences but lack localized sites
EThcD_OPair_Results_localized <- EThcD_OPair_Results |> 
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

EThcD_OPair_Results_nolocalized <- EThcD_OPair_Results |> 
  filter(!Confidence.Level %in% c('Level1', 'Level1b')) |> 
  filter(!str_detect(Peptide, 'C'))

EThcD_OPair_Results_bonafide <- bind_rows(
  EThcD_OPair_Results_localized,
  EThcD_OPair_Results_nolocalized
)

write_csv(
  EThcD_OPair_Results_bonafide,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlyco_HEK293T/OGlycoTM_EThcD_OPair_TMT_Search_1/Eclipse_LF_OGlycoTM_HEK293T_OG_4_10062025_HCD_pd_EThcD_Test_ppm/EThcD_OPair_Results_bonafide.csv'
)

# import filtered results
library(tidyverse)

EThcD_OPair_Results_bonafide <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlyco_HEK293T/OGlycoTM_EThcD_OPair_TMT_Search_1/Eclipse_LF_OGlycoTM_HEK293T_OG_4_10062025_HCD_pd_EThcD_Test_ppm/EThcD_OPair_Results_bonafide.csv'
)

# check number of total glycoPSM and unique glycoprotein
Total_GlycoPSM <- EThcD_OPair_Results_bonafide |> 
  nrow()

Unique_Glycoprotein <- EThcD_OPair_Results_bonafide |> 
  distinct(Protein.ID) |> 
  nrow()

GlycoPSM_Glycoprotein_Results <- tibble(
  category = c('Total GlycoPSM', 'Unique Glycoprotein'),
  number = c(440, 172)
)

# check number of localized and non-localized glycoPSM
Localized_GlycoPSM <- EThcD_OPair_Results_bonafide |> 
  filter(Confidence.Level %in% c('Level1', 'Level1b')) |> 
  nrow()

Nonlocalized_GlycoPSM <- EThcD_OPair_Results_bonafide |> 
  filter(!Confidence.Level %in% c('Level1', 'Level1b')) |> 
  nrow()

Localized_Nonlocalized_Results <- tibble(
  category = c('Localized', 'Nonlocalized'),
  number = c(72, 368)
)

# barplot
# glycoPSM & glycoprotein
GlycoPSM_Glycoprotein_Fraction_4 <- GlycoPSM_Glycoprotein_Results |> 
  ggplot() +
  geom_bar(
    aes(x = category, y = number),
    stat = 'identity'
  ) +
  geom_text(
    aes(x = category, y = number, label = number),
    size = 3, vjust = -0.1
  ) +
  labs(x = '', y = 'Number') +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 7, angle = 30, hjust = 1, color = 'black'),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/GlycoPSM_Glycoprotein_Fraction_4.tif',
  plot = GlycoPSM_Glycoprotein_Fraction_4,
  height = 3, width = 2,
  dpi = 300
)

# glycosylation site
Localized_Nonlocalized_Fraction_4 <- Localized_Nonlocalized_Results |> 
  ggplot() +
  geom_bar(
    aes(x = category, y = number),
    stat = 'identity'
  ) +
  geom_text(
    aes(x = category, y = number, label = number),
    size = 3, vjust = -0.1
  ) +
  labs(x = '', y = 'Number of PSM') +
  theme_bw() +
  theme(
    axis.text.x = element_text(size = 7, angle = 30, hjust = 1, color = 'black'),
    axis.text.y = element_text(size = 7, color = "black"),
    axis.title = element_text(size = 8, color = 'black')
  )

ggsave(
  filename = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/Localized_Nonlocalized_Fraction_4.tif',
  plot = Localized_Nonlocalized_Fraction_4,
  height = 3, width = 2,
  dpi = 300
)
