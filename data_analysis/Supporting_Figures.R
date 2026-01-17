# Supporting Figures for OGlycoTM Manuscript

library(tidyverse)
library(ggpubr)
library(rstatix)

# Source data paths, colors, and differential analysis results
source('data_source_DE.R')

# Figure S2 ---------------------------------------------------------------
# Violin boxplot comparing logFC of O-GlcNAc proteins in specific functional categories
# Categories: DNA binding, RNA binding, Transcription
# Statistical test: Kolmogorov-Smirnov test (compares distribution shapes)
# (Originally Figure 3C)

# Read GO enrichment results for common O-GlcNAc proteins
common_OGlcNAc_GO <- read_csv(
  paste0(source_file_path, 'enrichment/Figure2C_common_OGlcNAc_GO.csv')
)

# Define GO terms for each functional category
DNA_binding_terms <- c(
  "chromatin DNA binding",
  "damaged DNA binding",
  "telomeric DNA binding",
  "core promoter sequence-specific DNA binding",
  "single-stranded DNA binding"
)

RNA_binding_terms <- c(
  "mRNA binding",
  "single-stranded RNA binding",
  "double-stranded RNA binding",
  "pre-mRNA binding",
  "telomerase RNA binding",
  "snRNA binding"
)

# Extract all transcription-related terms using pattern matching
transcription_terms <- common_OGlcNAc_GO %>%
  filter(str_detect(Description, regex("transcription", ignore_case = TRUE))) %>%
  pull(Description)

# Extract proteins for each category
extract_proteins <- function(go_df, terms) {
  go_df %>%
    filter(Description %in% terms) %>%
    pull(geneID) %>%
    paste(collapse = "/") %>%
    str_split("/") %>%
    unlist() %>%
    unique()
}

DNA_binding_proteins <- extract_proteins(common_OGlcNAc_GO, DNA_binding_terms)
RNA_binding_proteins <- extract_proteins(common_OGlcNAc_GO, RNA_binding_terms)
transcription_proteins <- extract_proteins(common_OGlcNAc_GO, transcription_terms)

cat("Figure S2 - DNA binding proteins:", length(DNA_binding_proteins), "\n")
cat("Figure S2 - RNA binding proteins:", length(RNA_binding_proteins), "\n")
cat("Figure S2 - Transcription proteins:", length(transcription_proteins), "\n")

# Create combined logFC data for each category
create_category_df <- function(proteins, category_name) {
  bind_rows(
    OGlcNAc_protein_DE_HEK293T %>%
      filter(Protein.ID %in% proteins) %>%
      select(Protein.ID, logFC) %>%
      mutate(CellType = "HEK293T"),
    OGlcNAc_protein_DE_HepG2 %>%
      filter(Protein.ID %in% proteins) %>%
      select(Protein.ID, logFC) %>%
      mutate(CellType = "HepG2"),
    OGlcNAc_protein_DE_Jurkat %>%
      filter(Protein.ID %in% proteins) %>%
      select(Protein.ID, logFC) %>%
      mutate(CellType = "Jurkat")
  ) %>%
    mutate(Category = category_name)
}

FigureS2_df <- bind_rows(
  create_category_df(DNA_binding_proteins, "DNA binding"),
  create_category_df(RNA_binding_proteins, "RNA binding"),
  create_category_df(transcription_proteins, "Transcription")
)

# Set factor levels
FigureS2_df <- FigureS2_df %>%
  mutate(
    CellType = factor(CellType, levels = c("HEK293T", "HepG2", "Jurkat")),
    Category = factor(Category, levels = c("DNA binding", "RNA binding", "Transcription"))
  )

# Perform KS tests for each category (between cell types)
perform_ks_tests <- function(df, category_name) {
  df_cat <- df %>% filter(Category == category_name)

  HEK293T_logFC <- df_cat %>% filter(CellType == "HEK293T") %>% pull(logFC)
  HepG2_logFC <- df_cat %>% filter(CellType == "HepG2") %>% pull(logFC)
  Jurkat_logFC <- df_cat %>% filter(CellType == "Jurkat") %>% pull(logFC)

  tribble(
    ~Category, ~group1, ~group2, ~p,
    category_name, "HEK293T", "HepG2",
    if(length(HEK293T_logFC) > 0 & length(HepG2_logFC) > 0) ks.test(HEK293T_logFC, HepG2_logFC)$p.value else NA,
    category_name, "HEK293T", "Jurkat",
    if(length(HEK293T_logFC) > 0 & length(Jurkat_logFC) > 0) ks.test(HEK293T_logFC, Jurkat_logFC)$p.value else NA,
    category_name, "HepG2", "Jurkat",
    if(length(HepG2_logFC) > 0 & length(Jurkat_logFC) > 0) ks.test(HepG2_logFC, Jurkat_logFC)$p.value else NA
  )
}

FigureS2_ks_results <- bind_rows(
  perform_ks_tests(FigureS2_df, "DNA binding"),
  perform_ks_tests(FigureS2_df, "RNA binding"),
  perform_ks_tests(FigureS2_df, "Transcription")
) %>%
  filter(!is.na(p)) %>%
  add_significance("p") %>%
  filter(p.signif != "ns")  # Remove non-significant comparisons

cat("\nFigure S2 - KS Test Results:\n")
print(FigureS2_ks_results)

# Add y positions for significance bars
FigureS2_ks_results <- FigureS2_ks_results %>%
  mutate(
    y.position = case_when(
      group1 == "HEK293T" & group2 == "HepG2" ~ 1.0,
      group1 == "HEK293T" & group2 == "Jurkat" ~ 1.4,
      group1 == "HepG2" & group2 == "Jurkat" ~ 1.2,
      TRUE ~ 1.0
    )
  )

# Create faceted violin boxplot
FigureS2 <- FigureS2_df %>%
  ggplot(aes(x = CellType, y = logFC)) +
  geom_violin(aes(fill = CellType), color = "transparent") +
  geom_boxplot(color = "black", outliers = FALSE, width = 0.2, linewidth = 0.3) +
  facet_wrap(~ Category, nrow = 1) +
  scale_fill_manual(values = colors_cell) +
  labs(
    x = "",
    y = expression(log[2]*"(Tuni/Ctrl)")
  ) +
  stat_pvalue_manual(
    data = FigureS2_ks_results,
    label = "p.signif",
    tip.length = 0,
    size = 3
  ) +
  coord_cartesian(ylim = c(-1.5, 1.7)) +
  theme_bw() +
  theme(
    panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
    panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
    axis.title = element_text(size = 7),
    axis.text.x = element_text(size = 6, color = "black", angle = 30, hjust = 1),
    axis.text.y = element_text(size = 7, color = "black"),
    strip.text = element_text(size = 6),
    strip.background = element_blank(),
    legend.position = "none"
  )

print(FigureS2)

ggsave(
  filename = paste0(figure_file_path, "Supporting_Figures/FigureS2.pdf"),
  plot = FigureS2,
  width = 2, height = 1.5, units = "in"
)

cat("\nFigure S2 saved to:", figure_file_path, "Supporting_Figures/\n")
