# Glycoprotein Classification Analysis
# Classify glycoproteins as exclusive O-GlcNAc, exclusive O-GalNAc, or dual-labeled

# import packages
library(tidyverse)

# source data
source("data_analysis/data_source.R")

# Define glycan composition filters
OGlcNAc_compositions <- c(
 "HexNAt(1) % 299.1230",
 "HexNAt(1)TMT6plex(1) % 528.2859"
)

OGalNAc_compositions <- c(
 "HexNAt(1)GAO_Methoxylamine(1) % 326.1339",
 "HexNAt(1)GAO_Methoxylamine(1)TMT6plex(1) % 555.2968"
)


# Function to classify proteins and separate peptides -------------------------

classify_glycoproteins <- function(df, cell_type) {
  # Identify O-GlcNAc and O-GalNAc proteins
 OGlcNAc_proteins <- df |>
   filter(Total.Glycan.Composition %in% OGlcNAc_compositions) |>
   pull(Protein.ID) |>
   unique()

 OGalNAc_proteins <- df |>
   filter(Total.Glycan.Composition %in% OGalNAc_compositions) |>
   pull(Protein.ID) |>
   unique()

 # Classify proteins
 dual_labeled <- intersect(OGlcNAc_proteins, OGalNAc_proteins)
 exclusive_OGlcNAc <- setdiff(OGlcNAc_proteins, OGalNAc_proteins)
 exclusive_OGalNAc <- setdiff(OGalNAc_proteins, OGlcNAc_proteins)

 # Separate peptides based on protein classification
 # Exclusive O-GlcNAc: peptides from proteins that only have O-GlcNAc
 exclusive_OGlcNAc_peptides <- df |>
   filter(Protein.ID %in% exclusive_OGlcNAc)

 # Exclusive O-GalNAc: peptides from proteins that only have O-GalNAc
 exclusive_OGalNAc_peptides <- df |>
   filter(Protein.ID %in% exclusive_OGalNAc)

 # Dual-labeled: all peptides from proteins that have both modifications
 dual_labeled_peptides <- df |>
   filter(Protein.ID %in% dual_labeled)

 # Save to CSV files
 write_csv(
   exclusive_OGlcNAc_peptides,
   paste0(source_file_path, "classification/exclusive_OGlcNAc_", cell_type, ".csv")
 )

 write_csv(
   exclusive_OGalNAc_peptides,
   paste0(source_file_path, "classification/exclusive_OGalNAc_", cell_type, ".csv")
 )

 write_csv(
   dual_labeled_peptides,
   paste0(source_file_path, "classification/dual_labeled_", cell_type, ".csv")
 )

 # Return summary statistics
 data.frame(
   Cell.Type = cell_type,
   Total.OGlcNAc = length(OGlcNAc_proteins),
   Total.OGalNAc = length(OGalNAc_proteins),
   Exclusive.OGlcNAc = length(exclusive_OGlcNAc),
   Exclusive.OGalNAc = length(exclusive_OGalNAc),
   Dual.Labeled = length(dual_labeled),
   Pct.Exclusive.OGlcNAc = round(length(exclusive_OGlcNAc) / length(OGlcNAc_proteins) * 100, 1),
   Pct.Exclusive.OGalNAc = round(length(exclusive_OGalNAc) / length(OGalNAc_proteins) * 100, 1)
 )
}


# Classify glycoproteins for each cell type -----------------------------------

summary_HEK293T <- classify_glycoproteins(OGlyco_HEK293T_bonafide, "HEK293T")
summary_HepG2 <- classify_glycoproteins(OGlyco_HepG2_bonafide, "HepG2")
summary_Jurkat <- classify_glycoproteins(OGlyco_Jurkat_bonafide, "Jurkat")

# Combine summary statistics
summary_table <- bind_rows(
 summary_HEK293T,
 summary_HepG2,
 summary_Jurkat
)

print(summary_table)


# Generate summary table figure -----------------------------------------------

# Reshape data for visualization
summary_long <- summary_table |>
 select(Cell.Type, Exclusive.OGlcNAc, Exclusive.OGalNAc, Dual.Labeled) |>
 pivot_longer(
   cols = -Cell.Type,
   names_to = "Category",
   values_to = "Count"
 ) |>
 mutate(
   Category = factor(
     Category,
     levels = c("Exclusive.OGlcNAc", "Dual.Labeled", "Exclusive.OGalNAc"),
     labels = c("Exclusive\nO-GlcNAc", "Dual\nLabeled", "Exclusive\nO-GalNAc")
   ),
   Cell.Type = factor(Cell.Type, levels = c("HEK293T", "HepG2", "Jurkat"))
 )

# Create table visualization
summary_figure <- ggplot(summary_long, aes(x = Category, y = Cell.Type)) +
 geom_tile(aes(fill = Category), color = "white", linewidth = 1) +
 geom_text(aes(label = Count), size = 4, fontface = "bold") +
 scale_fill_manual(
   values = c(
     "Exclusive\nO-GlcNAc" = "#F39B7F",
     "Dual\nLabeled" = "#8491B4",
     "Exclusive\nO-GalNAc" = "#4DBBD5"
   )
 ) +
 labs(
   title = "Glycoprotein Classification by Cell Type",
   x = NULL,
   y = NULL
 ) +
 theme_minimal() +
 theme(
   plot.title = element_text(hjust = 0.5, size = 10, face = "bold"),
   axis.text = element_text(color = "black", size = 9),
   axis.text.x = element_text(lineheight = 0.9),
   legend.position = "none",
   panel.grid = element_blank()
 )

ggsave(
 filename = paste0(figure_file_path, "glycoprotein_classification_summary.pdf"),
 plot = summary_figure,
 width = 4, height = 2.5, units = "in"
)

# Also create a more detailed table with percentages
summary_detail <- summary_table |>
 mutate(
   `O-GlcNAc (Exclusive/Total)` = paste0(Exclusive.OGlcNAc, "/", Total.OGlcNAc, " (", Pct.Exclusive.OGlcNAc, "%)"),
   `O-GalNAc (Exclusive/Total)` = paste0(Exclusive.OGalNAc, "/", Total.OGalNAc, " (", Pct.Exclusive.OGalNAc, "%)"),
   `Dual Labeled` = as.character(Dual.Labeled)
 ) |>
 select(Cell.Type, `O-GlcNAc (Exclusive/Total)`, `O-GalNAc (Exclusive/Total)`, `Dual Labeled`)

# Reshape for detailed table visualization
summary_detail_long <- summary_detail |>
 pivot_longer(
   cols = -Cell.Type,
   names_to = "Metric",
   values_to = "Value"
 ) |>
 mutate(
   Metric = factor(Metric, levels = c("O-GlcNAc (Exclusive/Total)", "Dual Labeled", "O-GalNAc (Exclusive/Total)")),
   Cell.Type = factor(Cell.Type, levels = c("HEK293T", "HepG2", "Jurkat")),
   Value = as.character(Value)
 )

# Create detailed table visualization
summary_detail_figure <- ggplot(summary_detail_long, aes(x = Metric, y = Cell.Type)) +
 geom_tile(fill = "grey95", color = "grey70", linewidth = 0.5) +
 geom_text(aes(label = Value), size = 2.5) +
 labs(
   title = "Glycoprotein Classification Summary",
   x = NULL,
   y = NULL
 ) +
 theme_minimal() +
 theme(
   plot.title = element_text(hjust = 0.5, size = 9, face = "bold"),
   axis.text = element_text(color = "black", size = 7),
   axis.text.x = element_text(angle = 0, hjust = 0.5),
   panel.grid = element_blank()
 )

ggsave(
 filename = paste0(figure_file_path, "glycoprotein_classification_detail.pdf"),
 plot = summary_detail_figure,
 width = 5, height = 2, units = "in"
)

cat("\n=== Output files generated ===\n")
cat("CSV files saved to:", source_file_path, "\n")
cat("- exclusive_OGlcNAc_HEK293T.csv\n")
cat("- exclusive_OGalNAc_HEK293T.csv\n")
cat("- dual_labeled_HEK293T.csv\n")
cat("- exclusive_OGlcNAc_HepG2.csv\n")
cat("- exclusive_OGalNAc_HepG2.csv\n")
cat("- dual_labeled_HepG2.csv\n")
cat("- exclusive_OGlcNAc_Jurkat.csv\n")
cat("- exclusive_OGalNAc_Jurkat.csv\n")
cat("- dual_labeled_Jurkat.csv\n")
cat("\nFigures saved to:", figure_file_path, "\n")
cat("- glycoprotein_classification_summary.pdf\n")
cat("- glycoprotein_classification_detail.pdf\n")
