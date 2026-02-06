# Figure S1: Reproducibility plots for O-GlcNAc protein quantification
# Shows correlation between biological replicates using log2 fold change
# Fold change calculated as: Tuni_X / mean(Ctrl_4, Ctrl_5, Ctrl_6)

library(tidyverse)
library(GGally)
library(grid)
library(devEMF)  # For EMF output

# Source data paths and colors
source('data_source.R')

# Load normalized protein data
OGlcNAc_protein_norm_HEK293T <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_HEK293T.csv')
)
OGlcNAc_protein_norm_HepG2 <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_HepG2.csv')
)
OGlcNAc_protein_norm_Jurkat <- read_csv(
  paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_Jurkat.csv')
)

# Figure S1A: HEK293T reproducibility
OGlcNAc_protein_repro_HEK293T <- OGlcNAc_protein_norm_HEK293T |>
  mutate(
    Ctrl_mean = (Intensity.Ctrl_4_sl_tmm + Intensity.Ctrl_5_sl_tmm + Intensity.Ctrl_6_sl_tmm) / 3,
    Rep1 = log2(Intensity.Tuni_1_sl_tmm / Ctrl_mean),
    Rep2 = log2(Intensity.Tuni_2_sl_tmm / Ctrl_mean),
    Rep3 = log2(Intensity.Tuni_3_sl_tmm / Ctrl_mean)
  ) |>
  select(Rep1, Rep2, Rep3)

# Figure S1B: HepG2 reproducibility
OGlcNAc_protein_repro_HepG2 <- OGlcNAc_protein_norm_HepG2 |>
  mutate(
    Ctrl_mean = (Intensity.Ctrl_4_sl_tmm + Intensity.Ctrl_5_sl_tmm + Intensity.Ctrl_6_sl_tmm) / 3,
    Rep1 = log2(Intensity.Tuni_1_sl_tmm / Ctrl_mean),
    Rep2 = log2(Intensity.Tuni_2_sl_tmm / Ctrl_mean),
    Rep3 = log2(Intensity.Tuni_3_sl_tmm / Ctrl_mean)
  ) |>
  select(Rep1, Rep2, Rep3)

# Figure S1C: Jurkat reproducibility
OGlcNAc_protein_repro_Jurkat <- OGlcNAc_protein_norm_Jurkat |>
  mutate(
    Ctrl_mean = (Intensity.Ctrl_4_sl_tmm + Intensity.Ctrl_5_sl_tmm + Intensity.Ctrl_6_sl_tmm) / 3,
    Rep1 = log2(Intensity.Tuni_1_sl_tmm / Ctrl_mean),
    Rep2 = log2(Intensity.Tuni_2_sl_tmm / Ctrl_mean),
    Rep3 = log2(Intensity.Tuni_3_sl_tmm / Ctrl_mean)
  ) |>
  select(Rep1, Rep2, Rep3)

# Create ggpairs plots for each cell type
FigureS1A <- ggpairs(
  OGlcNAc_protein_repro_HEK293T,
  upper = list(continuous = wrap("cor", size = 2, color = "black")),
  lower = list(continuous = wrap("points", alpha = 0.5, size = 0.3)),
  diag = list(continuous = wrap("densityDiag", fill = "#4DBBD5", alpha = 0.7)),
  title = "HEK293T"
) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 6, hjust = 0.5, color = "black"),
    axis.text = element_text(size = 5, color = "black"),
    axis.title = element_text(size = 5, color = "black"),
    strip.text = element_text(size = 5, color = "black")
  )

FigureS1B <- ggpairs(
  OGlcNAc_protein_repro_HepG2,
  upper = list(continuous = wrap("cor", size = 2, color = "black")),
  lower = list(continuous = wrap("points", alpha = 0.5, size = 0.3)),
  diag = list(continuous = wrap("densityDiag", fill = "#F39B7F", alpha = 0.7)),
  title = "HepG2"
) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 6, hjust = 0.5, color = "black"),
    axis.text = element_text(size = 5, color = "black"),
    axis.title = element_text(size = 5, color = "black"),
    strip.text = element_text(size = 5, color = "black")
  )

FigureS1C <- ggpairs(
  OGlcNAc_protein_repro_Jurkat,
  upper = list(continuous = wrap("cor", size = 2, color = "black")),
  lower = list(continuous = wrap("points", alpha = 0.5, size = 0.3)),
  diag = list(continuous = wrap("densityDiag", fill = "#00A087", alpha = 0.7)),
  title = "Jurkat"
) +
  theme_bw() +
  theme(
    plot.title = element_text(size = 6, hjust = 0.5, color = "black"),
    axis.text = element_text(size = 5, color = "black"),
    axis.title = element_text(size = 5, color = "black"),
    strip.text = element_text(size = 5, color = "black")
  )

# Save combined plot using grid viewports
pdf(paste0(figure_file_path, "FigureS1/FigureS1_reproducibility.pdf"),
    width = 6, height = 2)

# Create layout with 3 columns
pushViewport(viewport(layout = grid.layout(1, 3)))

# Plot HEK293T in first column
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
print(FigureS1A, newpage = FALSE)
popViewport()

# Plot HepG2 in second column
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
print(FigureS1B, newpage = FALSE)
popViewport()

# Plot Jurkat in third column
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
print(FigureS1C, newpage = FALSE)
popViewport()

popViewport()
dev.off()

# Save as EMF format
emf(paste0(figure_file_path, "FigureS1/FigureS1_reproducibility.emf"),
    width = 6, height = 2)

# Create layout with 3 columns
pushViewport(viewport(layout = grid.layout(1, 3)))

# Plot HEK293T in first column
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
print(FigureS1A, newpage = FALSE)
popViewport()

# Plot HepG2 in second column
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
print(FigureS1B, newpage = FALSE)
popViewport()

# Plot Jurkat in third column
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
print(FigureS1C, newpage = FALSE)
popViewport()

popViewport()
dev.off()

# Save as TIFF format (600 DPI)
tiff(paste0(figure_file_path, "FigureS1/FigureS1_reproducibility.tiff"),
     width = 6, height = 2, units = "in", res = 600, compression = "lzw")

# Create layout with 3 columns
pushViewport(viewport(layout = grid.layout(1, 3)))

# Plot HEK293T in first column
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
print(FigureS1A, newpage = FALSE)
popViewport()

# Plot HepG2 in second column
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
print(FigureS1B, newpage = FALSE)
popViewport()

# Plot Jurkat in third column
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
print(FigureS1C, newpage = FALSE)
popViewport()

popViewport()
dev.off()

# Display in RStudio
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 3)))
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 1))
print(FigureS1A, newpage = FALSE)
popViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 2))
print(FigureS1B, newpage = FALSE)
popViewport()
pushViewport(viewport(layout.pos.row = 1, layout.pos.col = 3))
print(FigureS1C, newpage = FALSE)
popViewport()
popViewport()

cat("\nFigure S1 (Reproducibility) saved to:", figure_file_path, "FigureS1/")
cat("\n  - FigureS1_reproducibility.pdf")
cat("\n  - FigureS1_reproducibility.emf")
cat("\n  - FigureS1_reproducibility.tiff\n")
