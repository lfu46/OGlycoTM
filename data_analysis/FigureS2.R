# Figure S2: HEK293T GO Enrichment Barplot
# Originally Figure 4B - Shows top 5 GO terms enriched in HEK293T upregulated O-GlcNAc proteins

library(tidyverse)
library(devEMF)  # For EMF output

# Source data paths and colors
source('data_source.R')

# Create output directory
dir.create(paste0(figure_file_path, "FigureS2"), showWarnings = FALSE)

# Define lighter pastel color for HEK293T gradient bars
color_HEK293T_light <- "#7DCDE5"

# Function to create gradient bar data
create_gradient_data_S2 <- function(df, n_segments = 50) {
  df <- df |> mutate(y_num = as.numeric(fct_reorder(Description, log_pvalue)))

  gradient_df <- do.call(rbind, lapply(1:nrow(df), function(i) {
    row <- df[i, ]
    segments <- data.frame(
      Description = row$Description,
      y_num = row$y_num,
      xmin = seq(0, row$log_pvalue, length.out = n_segments + 1)[-(n_segments + 1)],
      xmax = seq(0, row$log_pvalue, length.out = n_segments + 1)[-1],
      segment = 1:n_segments
    )
    segments$alpha_val <- seq(0.9, 0.3, length.out = n_segments)
    segments
  }))
  gradient_df
}

# Load HEK293T GO enrichment results
OGlcNAc_up_HEK293T_GO_result <- read_csv(
  paste0(source_file_path, 'enrichment/OGlcNAc_up_HEK293T_GO.csv')
)

# Figure S2 - HEK293T Barplot with gradient
figureS2_data <- OGlcNAc_up_HEK293T_GO_result |>
  filter(
    Description %in% c(
      'ribonucleotide metabolic process',
      'response to retinoic acid',
      'ribosome biogenesis',
      'forebrain generation of neurons',
      'protein folding'
    )
  ) |>
  mutate(
    Description = case_when(
      Description == "ribonucleotide metabolic process" ~ "Ribonucleotide metabolism",
      Description == "response to retinoic acid" ~ "Response to retinoic acid",
      Description == "ribosome biogenesis" ~ "Ribosome biogenesis",
      Description == "forebrain generation of neurons" ~ "Forebrain neuron generation",
      Description == "protein folding" ~ "Protein folding",
      TRUE ~ Description
    ),
    log_pvalue = -log10(pvalue)
  )

figureS2_gradient <- create_gradient_data_S2(figureS2_data)

FigureS2 <- ggplot() +
  geom_rect(
    data = figureS2_gradient,
    aes(xmin = xmin, xmax = xmax,
        ymin = y_num - 0.35, ymax = y_num + 0.35,
        alpha = alpha_val),
    fill = color_HEK293T_light
  ) +
  geom_text(
    data = figureS2_data |> mutate(y_num = as.numeric(fct_reorder(Description, log_pvalue))),
    aes(label = Description, x = 0.05, y = y_num),
    hjust = 0, size = 2.8, color = "black"
  ) +
  scale_alpha_identity() +
  scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
  scale_y_continuous(breaks = 1:nrow(figureS2_data), labels = NULL) +
  labs(x = bquote("-"*log[10]*"("*italic(P)~Value*")"), y = "") +
  theme_classic() +
  theme(
    axis.title.x = element_text(size = 9),
    axis.text.x = element_text(color = "black", size = 9),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank()
  )

# Save as EMF format
emf(
  file = paste0(figure_file_path, 'FigureS2/FigureS2.emf'),
  height = 1.5, width = 2
)
print(FigureS2)
dev.off()

cat("\nFigure S2 saved to:", figure_file_path, "FigureS2/FigureS2.emf\n")
