library(tidyverse)
library(eulerr)

source_file_path <- "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/"
default_output_dir <- "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/Figure6_OGalNAc"
fallback_output_dir <- "/Users/longpingfu/Downloads/OGlycoTM/data_analysis/Figure6_OGalNAc_preview/secretory_downregulated"

colors_cell <- c(
  "HEK293T" = "#4DBBD5",
  "HepG2" = "#F39B7F",
  "Jurkat" = "#00A087"
)

fc_cutoff <- log2(1.5)
cell_types <- c("HEK293T", "HepG2", "Jurkat")

if (!dir.exists(default_output_dir)) {
  dir.create(default_output_dir, recursive = TRUE, showWarnings = FALSE)
}

output_dir <- if (file.access(default_output_dir, 2) == 0) {
  default_output_dir
} else {
  dir.create(fallback_output_dir, recursive = TRUE, showWarnings = FALSE)
  fallback_output_dir
}

message("Writing outputs to: ", output_dir)

feature_summary <- read_tsv(
  paste0(source_file_path, "reference/uniprot_positional_features_OGalNAc.tsv"),
  show_col_types = FALSE
) |>
  mutate(description = coalesce(description, "")) |>
  group_by(Protein.ID) |>
  summarize(
    has_signal = any(feature_type == "Signal"),
    has_tm = any(feature_type == "Transmembrane"),
    has_extracellular = any(
      feature_type == "Topological domain" &
        str_detect(str_to_lower(description), "extracellular|lumen|luminal")
    ),
    secretory = has_signal | has_tm | has_extracellular,
    .groups = "drop"
  )

load_de <- function(cell_type) {
  read_csv(
    paste0(source_file_path, "differential_analysis/OGalNAc_protein_DE_", cell_type, ".csv"),
    show_col_types = FALSE
  ) |>
    mutate(CellType = cell_type) |>
    left_join(feature_summary, by = "Protein.ID")
}

de_list <- set_names(map(cell_types, load_de), cell_types)

count_summary <- map_dfr(cell_types, function(cell_type) {
  de <- de_list[[cell_type]]

  tibble(
    CellType = cell_type,
    total_down = sum(de$logFC < -fc_cutoff & de$adj.P.Val < 0.05, na.rm = TRUE),
    secretory_down = sum(
      de$logFC < -fc_cutoff & de$adj.P.Val < 0.05 & de$secretory,
      na.rm = TRUE
    )
  )
})

write_csv(count_summary, file.path(output_dir, "OGalNAc_secretory_down_count_summary.csv"))

sec_down_tables <- map(
  de_list,
  ~ .x |>
    filter(logFC < -fc_cutoff, adj.P.Val < 0.05, secretory) |>
    distinct(Protein.ID, Gene, .keep_all = TRUE) |>
    arrange(logFC)
)

sec_down_list <- map(sec_down_tables, ~ pull(.x, Gene))

euler_fit <- euler(sec_down_list)

euler_plot <- plot(
  euler_fit,
  fills = list(fill = colors_cell, alpha = 0.5),
  edges = list(col = "white", lwd = 2),
  labels = list(font = 2, cex = 0.9),
  quantities = list(font = 1, cex = 0.85)
)

pdf(file.path(output_dir, "OGalNAc_secretory_down_euler.pdf"), width = 3.5, height = 3)
print(euler_plot)
dev.off()

png(
  file.path(output_dir, "OGalNAc_secretory_down_euler.png"),
  width = 3.5,
  height = 3,
  units = "in",
  res = 300
)
print(euler_plot)
dev.off()

all_three <- reduce(sec_down_list, intersect)
hek_hep_only <- setdiff(intersect(sec_down_list$HEK293T, sec_down_list$HepG2), sec_down_list$Jurkat)
hek_jur_only <- setdiff(intersect(sec_down_list$HEK293T, sec_down_list$Jurkat), sec_down_list$HepG2)
hep_jur_only <- setdiff(intersect(sec_down_list$HepG2, sec_down_list$Jurkat), sec_down_list$HEK293T)

overlap_summary <- tribble(
  ~Overlap, ~GeneCount, ~Genes,
  "All 3 cell types", length(all_three), paste(sort(all_three), collapse = ", "),
  "HEK293T + HepG2 only", length(hek_hep_only), paste(sort(hek_hep_only), collapse = ", "),
  "HEK293T + Jurkat only", length(hek_jur_only), paste(sort(hek_jur_only), collapse = ", "),
  "HepG2 + Jurkat only", length(hep_jur_only), paste(sort(hep_jur_only), collapse = ", ")
)

write_csv(overlap_summary, file.path(output_dir, "OGalNAc_secretory_down_overlap_summary.csv"))

shared_genes <- sort(names(which(table(unlist(sec_down_list)) >= 2)))

heatmap_data <- map_dfr(cell_types, function(cell_type) {
  de_list[[cell_type]] |>
    filter(Gene %in% shared_genes) |>
    distinct(Protein.ID, Gene, .keep_all = TRUE) |>
    transmute(
      CellType = cell_type,
      Protein.ID,
      Gene,
      logFC,
      adj.P.Val,
      is_down = logFC < -fc_cutoff & adj.P.Val < 0.05,
      sig = case_when(
        logFC < -fc_cutoff & adj.P.Val < 0.001 ~ "***",
        logFC < -fc_cutoff & adj.P.Val < 0.01 ~ "**",
        logFC < -fc_cutoff & adj.P.Val < 0.05 ~ "*",
        TRUE ~ ""
      )
    )
})

gene_order <- heatmap_data |>
  group_by(Gene) |>
  summarize(
    n_down = sum(is_down),
    mean_logFC = mean(logFC),
    .groups = "drop"
  ) |>
  arrange(desc(n_down), mean_logFC) |>
  pull(Gene)

heatmap_data <- heatmap_data |>
  mutate(
    Gene = factor(Gene, levels = rev(gene_order)),
    CellType = factor(CellType, levels = cell_types)
  )

write_csv(heatmap_data, file.path(output_dir, "OGalNAc_secretory_down_overlap_heatmap_source.csv"))

p_heatmap <- ggplot(heatmap_data, aes(x = CellType, y = Gene, fill = logFC)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = sig), size = 3, vjust = 0.75) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,
    limits = c(-4, 1),
    name = expression(log[2] * "(Tuni/Ctrl)")
  ) +
  scale_x_discrete(position = "top") +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x = element_text(face = "bold", color = "black"),
    axis.text.y = element_text(face = "italic", color = "black"),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.key.height = unit(0.8, "cm"),
    legend.key.width = unit(0.3, "cm")
  )

ggsave(
  file.path(output_dir, "OGalNAc_secretory_down_logFC_heatmap.pdf"),
  p_heatmap,
  width = 4.2,
  height = 4
)
ggsave(
  file.path(output_dir, "OGalNAc_secretory_down_logFC_heatmap.png"),
  p_heatmap,
  width = 4.2,
  height = 4,
  dpi = 300
)

category_levels <- c(
  "ER quality control / ERAD",
  "ER-Golgi transport / organelle homeostasis",
  "Glycan processing / Golgi enzymes",
  "Cell surface receptors / adhesion",
  "Solute transport",
  "Other secretory proteins"
)

category_map <- tribble(
  ~Gene, ~Category,
  "EDEM1", "ER quality control / ERAD",
  "OS9", "ER quality control / ERAD",
  "PRKCSH", "ER quality control / ERAD",
  "SEL1L", "ER quality control / ERAD",
  "TXNDC5", "ER quality control / ERAD",
  "MIA3", "ER-Golgi transport / organelle homeostasis",
  "TMEM165", "ER-Golgi transport / organelle homeostasis",
  "TGOLN2", "ER-Golgi transport / organelle homeostasis",
  "GPR108", "ER-Golgi transport / organelle homeostasis",
  "GALNT2", "Glycan processing / Golgi enzymes",
  "B4GALT4", "Glycan processing / Golgi enzymes",
  "B3GALT6", "Glycan processing / Golgi enzymes",
  "XXYLT1", "Glycan processing / Golgi enzymes",
  "MAN1A2", "Glycan processing / Golgi enzymes",
  "MANEA", "Glycan processing / Golgi enzymes",
  "CPD", "Glycan processing / Golgi enzymes",
  "TFRC", "Cell surface receptors / adhesion",
  "PTPRC", "Cell surface receptors / adhesion",
  "SPN", "Cell surface receptors / adhesion",
  "IGSF8", "Cell surface receptors / adhesion",
  "CD99", "Cell surface receptors / adhesion",
  "LDLR", "Cell surface receptors / adhesion",
  "EFNB1", "Cell surface receptors / adhesion",
  "RTN4RL2", "Cell surface receptors / adhesion",
  "TRGC2", "Cell surface receptors / adhesion",
  "SLC29A1", "Solute transport",
  "SLC39A10", "Solute transport",
  "OGN", "Other secretory proteins",
  "QSOX2", "Other secretory proteins",
  "TXNDC15", "Other secretory proteins"
)

jurkat_sec_down <- sec_down_tables$Jurkat |>
  left_join(category_map, by = "Gene") |>
  mutate(
    Category = replace_na(Category, "Other secretory proteins"),
    Category = factor(Category, levels = category_levels)
  ) |>
  arrange(Category, logFC) |>
  select(
    Protein.ID,
    Gene,
    Protein.Description,
    logFC,
    adj.P.Val,
    has_signal,
    has_tm,
    has_extracellular,
    Category
  )

write_csv(
  jurkat_sec_down,
  file.path(output_dir, "OGalNAc_secretory_down_categorized_Jurkat.csv")
)

category_summary <- jurkat_sec_down |>
  count(Category, name = "ProteinCount")

write_csv(
  category_summary,
  file.path(output_dir, "OGalNAc_secretory_down_category_counts.csv")
)

gene_levels <- jurkat_sec_down |>
  arrange(Category, logFC) |>
  pull(Gene)

jurkat_plot_data <- jurkat_sec_down |>
  mutate(Gene = factor(Gene, levels = rev(gene_levels)))

p_category <- ggplot(jurkat_plot_data, aes(x = logFC, y = Gene, color = Category)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey70", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = logFC, yend = Gene), linewidth = 0.8, show.legend = FALSE) +
  geom_point(size = 2.4, show.legend = FALSE) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_color_manual(
    values = c(
      "ER quality control / ERAD" = "#1B9E77",
      "ER-Golgi transport / organelle homeostasis" = "#7570B3",
      "Glycan processing / Golgi enzymes" = "#D95F02",
      "Cell surface receptors / adhesion" = "#E7298A",
      "Solute transport" = "#66A61E",
      "Other secretory proteins" = "#A6761D"
    )
  ) +
  labs(
    x = expression(log[2] * "(Tuni/Ctrl)"),
    y = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(
    strip.placement = "outside",
    strip.background = element_blank(),
    strip.text.y.left = element_text(angle = 0, face = "bold"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = element_text(face = "italic", color = "black"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(output_dir, "OGalNAc_secretory_down_category_plot.pdf"),
  p_category,
  width = 5.4,
  height = 6.1
)
ggsave(
  file.path(output_dir, "OGalNAc_secretory_down_category_plot.png"),
  p_category,
  width = 5.4,
  height = 6.1,
  dpi = 300
)

gsea_lookup <- bind_rows(
  read_csv(
    paste0(source_file_path, "enrichment/OGalNAc_exclusive_GSEA_CC_Jurkat.csv"),
    show_col_types = FALSE
  ) |>
    mutate(Source = "Exclusive OGalNAc", Ontology = "CC"),
  read_csv(
    paste0(source_file_path, "enrichment/OGalNAc_GSEA_CC_Jurkat.csv"),
    show_col_types = FALSE
  ) |>
    mutate(Source = "All OGalNAc", Ontology = "CC"),
  read_csv(
    paste0(source_file_path, "enrichment/OGalNAc_exclusive_GSEA_BP_Jurkat.csv"),
    show_col_types = FALSE
  ) |>
    mutate(Source = "Exclusive OGalNAc", Ontology = "BP")
)

gsea_key_terms <- tribble(
  ~Source, ~Ontology, ~ID, ~Description,
  "Exclusive OGalNAc", "CC", "GO:0031984", "organelle subcompartment",
  "Exclusive OGalNAc", "CC", "GO:0005789", "endoplasmic reticulum membrane",
  "Exclusive OGalNAc", "CC", "GO:0005783", "endoplasmic reticulum",
  "All OGalNAc", "CC", "GO:0000139", "Golgi membrane",
  "Exclusive OGalNAc", "CC", "GO:0012505", "endomembrane system",
  "Exclusive OGalNAc", "BP", "GO:1901135", "carbohydrate derivative metabolic process"
) |>
  left_join(
    gsea_lookup |>
      select(Source, Ontology, ID, Description, NES, p.adjust, setSize),
    by = c("Source", "Ontology", "ID", "Description")
  ) |>
  mutate(Significant = if_else(!is.na(p.adjust) & p.adjust < 0.05, "yes", "no")) |>
  arrange(Source, Ontology, NES)

write_csv(
  gsea_key_terms,
  file.path(output_dir, "OGalNAc_GSEA_secretory_key_terms_Jurkat.csv")
)

message("Count summary:")
print(count_summary)

message("Overlap summary:")
print(overlap_summary, n = Inf)

message("Jurkat category counts:")
print(category_summary, n = Inf)

message("Key GSEA terms:")
print(gsea_key_terms, n = Inf)
