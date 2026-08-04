library(tidyverse)

source_file_path <- "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/"
default_output_dir <- "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/Figure6_OGalNAc"
fallback_output_dir <- "/Users/longpingfu/Downloads/OGlycoTM/data_analysis/Figure6_OGalNAc_preview/secretory_downregulated"

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

features <- read_tsv(
  paste0(source_file_path, "reference/uniprot_positional_features_OGalNAc.tsv"),
  show_col_types = FALSE
) |>
  mutate(description = coalesce(description, ""))

site_de_jurkat <- read_csv(
  paste0(source_file_path, "differential_analysis/OGalNAc_site_DE_Jurkat.csv"),
  show_col_types = FALSE
)

glycodomain_annotation <- read_csv(
  paste0(source_file_path, "annotation/OGalNAc_site_glycodomain_annotation.csv"),
  show_col_types = FALSE
)

example_sites <- tribble(
  ~Example, ~CellType, ~Protein.ID, ~Gene, ~modified_residue, ~site_number, ~Literature_PMID, ~Evidence_note, ~Interpretation_note,
  "TFRC T104", "Jurkat", "P02786", "TFRC", "T", 104, "1421757; 8286753", "Exact site-level literature support for Thr104 on transferrin receptor; mutation of Thr104 disrupts protection from proteolytic cleavage.", "Current site falls in the extracellular stem region immediately after the transmembrane helix, consistent with a cleavage-sensitive mucin-like context.",
  "PTPRC T139/S146", "Jurkat", "P08575", "PTPRC", "T", 139, "7501023; 19920154; 30582698", "CD45 O-glycans are functionally established at the protein-region level, but not specifically at Thr139.", "Current site lies in the extracellular disordered N-terminal region that carries dense mucin-type O-glycans on Jurkat CD45.",
  "PTPRC T139/S146", "Jurkat", "P08575", "PTPRC", "S", 146, "7501023; 19920154; 30582698", "CD45 O-glycans are functionally established at the protein-region level, but not specifically at Ser146.", "Current site lies in the extracellular disordered N-terminal region that carries dense mucin-type O-glycans on Jurkat CD45."
) |>
  mutate(Example = factor(Example, levels = c("TFRC T104", "PTPRC T139/S146")))

example_summary <- example_sites |>
  left_join(
    site_de_jurkat,
    by = c("Protein.ID", "Gene", "modified_residue", "site_number")
  ) |>
  left_join(
    glycodomain_annotation |>
      filter(CellType == "Jurkat") |>
      select(
        Protein.ID,
        Gene,
        modified_residue,
        site_number,
        topology_class,
        domain_context,
        region_class,
        tm_distance
      ),
    by = c("Protein.ID", "Gene", "modified_residue", "site_number")
  ) |>
  select(
    Example,
    CellType,
    Protein.ID,
    Gene,
    Site = modified_residue,
    site_number,
    logFC,
    adj.P.Val,
    topology_class,
    domain_context,
    region_class,
    tm_distance,
    Literature_PMID,
    Evidence_note,
    Interpretation_note
  ) |>
  mutate(Example = factor(Example, levels = levels(example_sites$Example)))

if (any(is.na(example_summary$logFC))) {
  stop("One or more example sites were not found in the Jurkat O-GalNAc site DE table.")
}

write_csv(
  example_summary,
  file.path(output_dir, "OGalNAc_secretory_example_site_summary.csv")
)

protein_lengths <- features |>
  filter(feature_type == "Protein_length") |>
  group_by(Protein.ID) |>
  summarize(protein_length = max(end, na.rm = TRUE), .groups = "drop")

plot_proteins <- example_sites |>
  distinct(Example, Protein.ID, Gene) |>
  left_join(protein_lengths, by = "Protein.ID") |>
  mutate(Example = factor(Example, levels = c("TFRC T104", "PTPRC T139/S146")))

feature_plot_df <- features |>
  filter(Protein.ID %in% plot_proteins$Protein.ID) |>
  mutate(
    track = case_when(
      feature_type == "Topological domain" ~ "Topology",
      feature_type %in% c("Signal", "Transmembrane") ~ "Membrane features",
      feature_type %in% c("Domain", "Repeat") ~ "Domains / regions",
      feature_type == "Region" &
        description %in% c("Disordered", "Ligand-binding", "Mediates interaction with SH3BP4") ~ "Domains / regions",
      TRUE ~ NA_character_
    )
  ) |>
  filter(!is.na(track)) |>
  left_join(plot_proteins, by = "Protein.ID") |>
  mutate(
    Example = factor(Example, levels = levels(example_sites$Example)),
    feature_label = case_when(
      feature_type == "Topological domain" &
        str_detect(str_to_lower(description), "extracellular|luminal|lumen") ~ "Extracellular / luminal",
      feature_type == "Topological domain" &
        str_detect(str_to_lower(description), "cytoplasmic|cytosolic") ~ "Cytoplasmic",
      feature_type == "Signal" ~ "Signal peptide",
      feature_type == "Transmembrane" ~ "TM helix",
      description == "Fibronectin type-III 1" ~ "FNIII-1",
      description == "Fibronectin type-III 2" ~ "FNIII-2",
      description == "Tyrosine-protein phosphatase 1" ~ "PTP-1",
      description == "Tyrosine-protein phosphatase 2" ~ "PTP-2",
      description == "Mediates interaction with SH3BP4" ~ "SH3BP4 bind.",
      TRUE ~ description
    ),
    feature_class = case_when(
      feature_type == "Topological domain" &
        str_detect(str_to_lower(description), "extracellular|luminal|lumen") ~ "Extracellular / luminal",
      feature_type == "Topological domain" &
        str_detect(str_to_lower(description), "cytoplasmic|cytosolic") ~ "Cytoplasmic",
      feature_type == "Signal" ~ "Signal peptide",
      feature_type == "Transmembrane" ~ "TM helix",
      feature_type == "Region" & description == "Disordered" ~ "Disordered region",
      TRUE ~ "Annotated domain / region"
    ),
    track_id = case_when(
      track == "Topology" ~ 1,
      track == "Membrane features" ~ 2,
      track == "Domains / regions" ~ 3
    ),
    label_x = (start + end) / 2
  )

feature_labels <- feature_plot_df |>
  filter(track == "Domains / regions") |>
  mutate(feature_label = str_trunc(feature_label, 18))

site_plot_df <- example_summary |>
  mutate(
    SiteLabel = paste0(Site, site_number),
    SiteText = paste0(Site, site_number, "\nlogFC=", sprintf("%.2f", logFC))
  ) |>
  left_join(plot_proteins, by = c("Example", "Protein.ID", "Gene")) |>
  group_by(Example) |>
  arrange(site_number, .by_group = TRUE) |>
  mutate(
    label_rank = row_number(),
    label_x = case_when(
      n() == 1 ~ site_number,
      label_rank %% 2 == 1 ~ site_number - 32,
      TRUE ~ site_number + 32
    ),
    label_y = 3.62 + (label_rank - 1) * 0.28,
    label_hjust = case_when(
      n() == 1 ~ 0.5,
      label_rank %% 2 == 1 ~ 1,
      TRUE ~ 0
    )
  ) |>
  ungroup()

protein_backbone <- plot_proteins |>
  mutate(y = 0.45)

domain_map <- ggplot() +
  geom_segment(
    data = protein_backbone,
    aes(x = 1, xend = protein_length, y = y, yend = y),
    linewidth = 1.2,
    color = "grey55"
  ) +
  geom_rect(
    data = feature_plot_df,
    aes(
      xmin = start,
      xmax = end,
      ymin = track_id - 0.22,
      ymax = track_id + 0.22,
      fill = feature_class
    ),
    color = "grey30",
    linewidth = 0.2
  ) +
  geom_text(
    data = feature_labels,
    aes(x = label_x, y = track_id, label = feature_label),
    size = 2.4,
    color = "black",
    check_overlap = TRUE
  ) +
  geom_segment(
    data = site_plot_df,
    aes(x = site_number, xend = site_number, y = 0.45, yend = 3.45),
    linewidth = 0.35,
    color = "grey45",
    linetype = "dashed"
  ) +
  geom_segment(
    data = site_plot_df,
    aes(x = site_number, xend = label_x, y = 3.45, yend = label_y - 0.08),
    linewidth = 0.3,
    color = "#D95F02"
  ) +
  geom_point(
    data = site_plot_df,
    aes(x = site_number, y = 3.45),
    size = 2.8,
    color = "#D95F02"
  ) +
  geom_text(
    data = site_plot_df,
    aes(x = label_x, y = label_y, label = SiteText, hjust = label_hjust),
    size = 2.7,
    fontface = "bold",
    color = "#D95F02",
    vjust = 0
  ) +
  geom_text(
    data = plot_proteins,
    aes(x = protein_length, y = 0.1, label = paste0("aa ", protein_length)),
    hjust = 1,
    size = 2.6
  ) +
  facet_wrap(~ Example, ncol = 1, scales = "free_x") +
  scale_fill_manual(
    values = c(
      "Extracellular / luminal" = "#A6CEE3",
      "Cytoplasmic" = "#CAB2D6",
      "Signal peptide" = "#FDBF6F",
      "TM helix" = "#FB9A99",
      "Disordered region" = "#B2DF8A",
      "Annotated domain / region" = "#FDD0A2"
    )
  ) +
  scale_y_continuous(
    breaks = c(1, 2, 3),
    labels = c("Topology", "Membrane", "Domains"),
    limits = c(0, 4.1),
    expand = expansion(mult = c(0.02, 0.04))
  ) +
  labs(
    x = "Residue position",
    y = NULL,
    fill = NULL
  ) +
  theme_bw(base_size = 10) +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    strip.text = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    panel.grid.minor = element_blank()
  )

ggsave(
  file.path(output_dir, "OGalNAc_secretory_example_domain_maps.pdf"),
  domain_map,
  width = 8,
  height = 5.4
)
ggsave(
  file.path(output_dir, "OGalNAc_secretory_example_domain_maps.png"),
  domain_map,
  width = 8,
  height = 5.4,
  dpi = 300
)

write_csv(
  feature_plot_df |>
    select(
      Example,
      Protein.ID,
      Gene,
      feature_type,
      start,
      end,
      description,
      track,
      feature_class,
      feature_label
    ),
  file.path(output_dir, "OGalNAc_secretory_example_domain_map_source.csv")
)

message("Example summary:")
print(example_summary, n = Inf)
