# import packages
library(tidyverse)

## O-Pair results
# O-GlcNAc
OPair_GlycoPSM_Glycoprotein_Results_OGlcNAc <- tibble(
  category = c('Total GlycoPSM', 'Unique Glycoprotein'),
  number = c(9242, 738)
)

# O-GalNAc
OPair_GlycoPSM_Glycoprotein_Results_OGalNAc <- tibble(
  category = c('Total GlycoPSM', 'Unique Glycoprotein'),
  number = c(1143, 235)
)

## MSFragger Glyco results
# O-GlcNAc
GlycoPSM_Glycoprotein_Results_OGlcNAc <- tibble(
  category = c('Total GlycoPSM', 'Unique Glycoprotein'),
  number = c(9059, 699)
)

# O-GalNAc
GlycoPSM_Glycoprotein_Results_OGalNAc <- tibble(
  category = c('Total GlycoPSM', 'Unique Glycoprotein'),
  number = c(878, 139)
)

# Combine all data
plot_data <- bind_rows(
  OPair_GlycoPSM_Glycoprotein_Results_OGlcNAc %>% mutate(glycan = 'O-GlcNAc', method = 'O-Pair'),
  OPair_GlycoPSM_Glycoprotein_Results_OGalNAc %>% mutate(glycan = 'O-GalNAc', method = 'O-Pair'),
  GlycoPSM_Glycoprotein_Results_OGlcNAc %>% mutate(glycan = 'O-GlcNAc', method = 'MSFragger Glyco'),
  GlycoPSM_Glycoprotein_Results_OGalNAc %>% mutate(glycan = 'O-GalNAc', method = 'MSFragger Glyco')
) %>%
  mutate(
    glycan = factor(glycan, levels = c('O-GlcNAc', 'O-GalNAc')),
    category = factor(category, levels = c('Total GlycoPSM', 'Unique Glycoprotein')),
    method = factor(method, levels = c('O-Pair', 'MSFragger Glyco'))
  )

# Publication-ready color palette
colors <- c('O-Pair' = '#2E8B9E', 'MSFragger Glyco' = '#C85A54')

# Create plot with broken y-axis
OPair_vs_MSFraggerGlyco <- ggplot(plot_data, aes(x = category, y = number, fill = method)) +
  geom_col(position = 'dodge', width = 0.7) +
  geom_text(aes(label = number), position = position_dodge(width = 0.7), 
            vjust = -0.1, size = 2.5) +
  facet_grid(. ~ glycan) +
  scale_fill_manual(values = colors) +
  ylim(0, 9500) +
  labs(
    x = NULL,
    y = 'Count',
    fill = 'Method',
    title = 'Glycoproteomics Results: O-Pair vs MSFragger Glyco'
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13),
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 7, color = 'black'),
    axis.text.y = element_text(size = 8, color = 'black'),
    axis.title.y = element_text(size = 8,vjust = 0.5),
    legend.position = 'bottom',
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    legend.key.size = unit(0.15, 'in'),
    strip.text = element_text(size = 8),
    plot.background = element_rect(fill = "white", color = NA)
  )

ggsave(
  filename = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OPair_vs_MSFraggerGlyco.png',
  plot = OPair_vs_MSFraggerGlyco,
  height = 3, width = 5, dpi = 300
)
