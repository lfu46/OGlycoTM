# Table of Contents: GO Dot Plot for Commonly Regulated O-GlcNAc Proteins

library(tidyverse)

# Source data paths and colors
source('data_source.R')

# =============================================================================
# Read GO enrichment results for commonly regulated proteins
# =============================================================================

# Load GO enrichment results
commonly_up_GO <- read_csv(
  paste0(source_file_path, 'enrichment/Figure3E_commonly_up_GO_total.csv'),
  show_col_types = FALSE
)

commonly_down_GO <- read_csv(
  paste0(source_file_path, 'enrichment/Figure3E_commonly_down_GO_total.csv'),
  show_col_types = FALSE
)

# =============================================================================
# Extract selected GO terms
# =============================================================================

# Selected terms:
# - 'regulation of translation' from commonly upregulated proteins
# - 'cytoplasmic stress granule' from commonly downregulated proteins

up_term <- commonly_up_GO |>
  filter(Description == "regulation of translation") |>
  mutate(Regulation = "up") |>
  dplyr::select(Description, pvalue, Count, Regulation)

down_term <- commonly_down_GO |>
  filter(Description == "cytoplasmic stress granule") |>
  mutate(Regulation = "down") |>
  dplyr::select(Description, pvalue, Count, Regulation)

# Print the selected terms
cat("Selected GO terms:\n")
cat("\nUpregulated:\n")
print(up_term)
cat("\nDownregulated:\n")
print(down_term)

# =============================================================================
# Prepare data for dot plot
# =============================================================================

# Combine the two terms and create full grid for dot plot
# We want to show both terms on y-axis, and up/down on x-axis
# Each term should only have a dot in its relevant column

GO_dotplot_data <- bind_rows(up_term, down_term) |>
  mutate(
    Regulation = factor(Regulation, levels = c("up", "down")),
    log_pvalue = -log10(pvalue),
    # Abbreviate long term names
    Description = case_when(
      Description == "regulation of translation" ~ "Translation reg.",
      Description == "cytoplasmic stress granule" ~ "Stress granule",
      TRUE ~ Description
    )
  )

# =============================================================================
# Create GO dot plot
# =============================================================================

# Color palette
dot_color <- "#8491B4"

TOC_GO_dotplot <- ggplot(GO_dotplot_data, aes(x = Regulation, y = Description)) +
  geom_point(aes(size = Count), color = dot_color) +
  scale_size_continuous(range = c(2, 4)) +
  scale_x_discrete(position = "bottom") +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.3, color = "gray85"),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 7, color = "black"),
    axis.text.y = element_text(size = 7, color = "black"),
    legend.position = "none",
    plot.margin = margin(5, 5, 5, 5)
  )

print(TOC_GO_dotplot)

# Save figure
ggsave(
  filename = paste0(figure_file_path, "Table_of_Contents/TOC_GO_dotplot.pdf"),
  plot = TOC_GO_dotplot,
  width = 1.3, height = 0.6, units = "in"
)

cat("\nTOC GO dot plot saved to:", figure_file_path, "Table_of_Contents/\n")

# =============================================================================
# Print summary statistics
# =============================================================================

cat("\n", paste(rep("=", 60), collapse = ""), "\n")
cat("Table of Contents - GO Terms Summary\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("\nCommon Regulation of O-GlcNAc Proteins:\n")
cat("\n1. Upregulated proteins enriched for:\n")
cat("   - ", up_term$Description, "\n")
cat("     Count: ", up_term$Count, ", p-value: ", format(up_term$pvalue, digits = 4), "\n")
cat("\n2. Downregulated proteins enriched for:\n")
cat("   - ", down_term$Description, "\n")
cat("     Count: ", down_term$Count, ", p-value: ", format(down_term$pvalue, digits = 4), "\n")
