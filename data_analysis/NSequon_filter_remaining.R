# Generate remaining filtered subfigures: Figure 4B, 4C + copy unchanged panels
# Run AFTER NSequon_filter_comparison.R

library(tidyverse)
library(introdataviz)

source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'
output_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/NSequon_filtered_figures/'
colors_cell <- c("HEK293T" = "#4DBBD5", "HepG2" = "#F39B7F", "Jurkat" = "#00A087")

# Load filtered protein lists
OGlcNAc_compositions <- c('HexNAt(1) % 299.1230', 'HexNAt(1)TMT6plex(1) % 528.2859')
apply_nsequon_filter <- function(bonafide_df) {
  bonafide_df |>
    filter(Total.Glycan.Composition %in% OGlcNAc_compositions) |>
    filter(!(Confidence.Level == "Level3" & Has.N.Glyc.Sequon == TRUE))
}

OGlyco_HEK293T_bonafide <- read_csv(paste0(source_file_path, 'filtered/OGlyco_HEK293T_bonafide.csv'), show_col_types = FALSE)
OGlyco_HepG2_bonafide <- read_csv(paste0(source_file_path, 'filtered/OGlyco_HepG2_bonafide.csv'), show_col_types = FALSE)
OGlyco_Jurkat_bonafide <- read_csv(paste0(source_file_path, 'filtered/OGlyco_Jurkat_bonafide.csv'), show_col_types = FALSE)

valid_HEK <- unique(apply_nsequon_filter(OGlyco_HEK293T_bonafide)$Protein.ID)
valid_HepG2 <- unique(apply_nsequon_filter(OGlyco_HepG2_bonafide)$Protein.ID)
valid_Jurkat <- unique(apply_nsequon_filter(OGlyco_Jurkat_bonafide)$Protein.ID)

# Load and filter DE results
OGlcNAc_protein_DE_HEK293T <- read_csv(paste0(source_file_path, 'differential_analysis/OGlcNAc_protein_DE_HEK293T.csv'), show_col_types = FALSE) |> filter(Protein.ID %in% valid_HEK)
OGlcNAc_protein_DE_HepG2 <- read_csv(paste0(source_file_path, 'differential_analysis/OGlcNAc_protein_DE_HepG2.csv'), show_col_types = FALSE) |> filter(Protein.ID %in% valid_HepG2)
OGlcNAc_protein_DE_Jurkat <- read_csv(paste0(source_file_path, 'differential_analysis/OGlcNAc_protein_DE_Jurkat.csv'), show_col_types = FALSE) |> filter(Protein.ID %in% valid_Jurkat)

# =============================================================================
# Figure 4B: Jurkat GO enrichment barplot (upregulated)
# =============================================================================
cat("=== Figure 4B ===\n")

library(clusterProfiler)
library(org.Hs.eg.db)

up_Jurkat <- OGlcNAc_protein_DE_Jurkat |> filter(logFC > 0.5, adj.P.Val < 0.05) |> pull(Protein.ID)
bg <- Reduce(union, list(OGlcNAc_protein_DE_HEK293T$Protein.ID, OGlcNAc_protein_DE_HepG2$Protein.ID, OGlcNAc_protein_DE_Jurkat$Protein.ID))

cat("Jurkat up:", length(up_Jurkat), " background:", length(bg), "\n")

up_Jurkat_GO <- enrichGO(gene = up_Jurkat, OrgDb = org.Hs.eg.db, universe = bg,
                          keyType = 'UNIPROT', ont = 'ALL', pvalueCutoff = 1, qvalueCutoff = 1)

figure4B_data <- up_Jurkat_GO@result |>
  filter(Description %in% c(
    'regulation of leukocyte proliferation', 'cell activation',
    'lymphocyte activation', 'leukocyte migration', 'regulation of T cell activation'
  )) |>
  mutate(
    Description = case_when(
      Description == "regulation of leukocyte proliferation" ~ "Leukocyte proliferation",
      Description == "cell activation" ~ "Cell activation",
      Description == "lymphocyte activation" ~ "Lymphocyte activation",
      Description == "leukocyte migration" ~ "Leukocyte migration",
      Description == "regulation of T cell activation" ~ "T cell activation",
      TRUE ~ Description
    ),
    log_pvalue = -log10(pvalue)
  )

cat("Figure 4B GO terms found:", nrow(figure4B_data), "\n")
if (nrow(figure4B_data) > 0) print(figure4B_data |> dplyr::select(Description, pvalue, Count, log_pvalue))

# Gradient bar function
create_gradient_data <- function(df, n_segments = 50) {
  df <- df |> mutate(y_num = as.numeric(fct_reorder(Description, log_pvalue)))
  do.call(rbind, lapply(1:nrow(df), function(i) {
    row <- df[i, ]
    segs <- data.frame(
      Description = row$Description, y_num = row$y_num,
      xmin = seq(0, row$log_pvalue, length.out = n_segments + 1)[-(n_segments + 1)],
      xmax = seq(0, row$log_pvalue, length.out = n_segments + 1)[-1],
      segment = 1:n_segments
    )
    segs$alpha_val <- seq(0.9, 0.3, length.out = n_segments)
    segs
  }))
}

if (nrow(figure4B_data) > 0) {
  gradient <- create_gradient_data(figure4B_data)
  figure4B <- ggplot() +
    geom_rect(data = gradient, aes(xmin = xmin, xmax = xmax, ymin = y_num - 0.35, ymax = y_num + 0.35, alpha = alpha_val), fill = "#4DC4B0") +
    geom_text(data = figure4B_data |> mutate(y_num = as.numeric(fct_reorder(Description, log_pvalue))),
              aes(label = Description, x = 0.05, y = y_num), hjust = 0, size = 2.8, color = "black") +
    scale_alpha_identity() +
    scale_x_continuous(expand = expansion(mult = c(0, 0.1))) +
    scale_y_continuous(breaks = 1:nrow(figure4B_data), labels = NULL) +
    labs(x = expression(-log[10]*"("*paste(italic(P), " Value")*")"), y = "") +
    theme_classic() +
    theme(axis.title.x = element_text(size = 9), axis.text.x = element_text(color = "black", size = 9),
          axis.text.y = element_blank(), axis.ticks.y = element_blank())
  ggsave(paste0(output_path, "Figure4B.pdf"), figure4B, height = 1.5, width = 2)
  cat("Figure 4B saved\n")
}

# =============================================================================
# Figure 4C: Jurkat example proteins (NFATC2, CTTN, SEC31A)
# =============================================================================
cat("\n=== Figure 4C ===\n")

OGlcNAc_protein_norm_HEK293T <- read_csv(paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_HEK293T.csv'), show_col_types = FALSE) |> filter(Protein.ID %in% valid_HEK)
OGlcNAc_protein_norm_HepG2 <- read_csv(paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_HepG2.csv'), show_col_types = FALSE) |> filter(Protein.ID %in% valid_HepG2)
OGlcNAc_protein_norm_Jurkat <- read_csv(paste0(source_file_path, 'normalization/OGlcNAc_protein_norm_Jurkat.csv'), show_col_types = FALSE) |> filter(Protein.ID %in% valid_Jurkat)

# NFATC2=Q13469, CTTN=Q14247, SEC31A=O94979
proteins_4C <- c("Q13469", "Q14247", "O94979")

# Check which survive the filter
for (p in proteins_4C) {
  gene <- c("NFATC2", "CTTN", "SEC31A")[match(p, proteins_4C)]
  cat(gene, "- HEK:", p %in% valid_HEK, " HepG2:", p %in% valid_HepG2, " Jurkat:", p %in% valid_Jurkat, "\n")
}

calculate_fc <- function(norm_data, cell_name, proteins) {
  norm_data |>
    filter(Protein.ID %in% proteins) |>
    dplyr::select(Protein.ID, Gene,
                  Intensity.Tuni_1_sl_tmm, Intensity.Tuni_2_sl_tmm, Intensity.Tuni_3_sl_tmm,
                  Intensity.Ctrl_4_sl_tmm, Intensity.Ctrl_5_sl_tmm, Intensity.Ctrl_6_sl_tmm) |>
    rowwise() |>
    mutate(
      mean_Ctrl = mean(c(Intensity.Ctrl_4_sl_tmm, Intensity.Ctrl_5_sl_tmm, Intensity.Ctrl_6_sl_tmm), na.rm = TRUE),
      log2FC_rep1 = log2(Intensity.Tuni_1_sl_tmm / mean_Ctrl),
      log2FC_rep2 = log2(Intensity.Tuni_2_sl_tmm / mean_Ctrl),
      log2FC_rep3 = log2(Intensity.Tuni_3_sl_tmm / mean_Ctrl)
    ) |>
    ungroup() |>
    dplyr::select(Protein.ID, Gene, log2FC_rep1, log2FC_rep2, log2FC_rep3) |>
    pivot_longer(cols = starts_with("log2FC"), names_to = "replicate", values_to = "log2FC") |>
    mutate(cell = cell_name)
}

fc_H <- calculate_fc(OGlcNAc_protein_norm_HEK293T, "HEK293T", proteins_4C)
fc_G <- calculate_fc(OGlcNAc_protein_norm_HepG2, "HepG2", proteins_4C)
fc_J <- calculate_fc(OGlcNAc_protein_norm_Jurkat, "Jurkat", proteins_4C)

fc_combined_4C <- bind_rows(fc_H, fc_G, fc_J) |>
  mutate(cell = factor(cell, levels = c("HEK293T", "HepG2", "Jurkat")),
         Gene = factor(Gene, levels = c("NFATC2", "CTTN", "SEC31A")))

# Remove HEK293T NFATC2 (not identified) and any empty combos
fc_combined_4C <- fc_combined_4C |>
  filter(!(Gene == "NFATC2" & cell == "HEK293T")) |>
  mutate(cell = droplevels(cell))

cat("Proteins with data after filter:\n")
print(fc_combined_4C |> group_by(Gene) |> summarise(cells = paste(unique(cell), collapse=", "), n = n()))

figure4C <- fc_combined_4C |>
  ggplot(aes(x = cell, y = log2FC, color = cell)) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  geom_hline(yintercept = 0.5, color = "black", linetype = "dashed", linewidth = 0.5) +
  geom_point(size = 2, position = position_jitter(width = 0.1, seed = 42)) +
  scale_color_manual(values = colors_cell) +
  scale_x_discrete(drop = TRUE) +
  facet_grid(. ~ Gene, scales = "free_x", space = "free_x") +
  labs(x = "", y = expression(log[2]*"(Tuni/Ctrl)")) +
  theme_classic() +
  theme(axis.title.y = element_text(size = 9),
        axis.text.x = element_text(color = "black", size = 9, angle = 90, hjust = 1),
        axis.text.y = element_text(color = "black", size = 9),
        strip.text = element_text(size = 9, face = "bold"),
        strip.background = element_blank(),
        legend.position = "none", panel.spacing = unit(0.3, "lines"))

ggsave(paste0(output_path, "Figure4C.pdf"), figure4C, height = 2, width = 2.5)
cat("Figure 4C saved (note: CTTN missing from Jurkat)\n")

# =============================================================================
# Figure 2B/2C: Copy unchanged O-GalNAc panels
# =============================================================================
cat("\n=== Copying unchanged O-GalNAc panels ===\n")
fig2_src <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/Figure2_new/'
for (f in c("Figure2B.pdf", "Figure2C.pdf")) {
  src <- paste0(fig2_src, f)
  dst <- paste0(output_path, f)
  if (file.exists(src)) {
    file.copy(src, dst, overwrite = TRUE)
    cat("Copied", f, "(unchanged - O-GalNAc)\n")
  }
}

# =============================================================================
# Figure 3C: UMAP (simplified using R uwot)
# =============================================================================
cat("\n=== Figure 3C (UMAP) ===\n")

library(uwot)

overlapping <- Reduce(intersect, list(
  OGlcNAc_protein_DE_HEK293T$Protein.ID,
  OGlcNAc_protein_DE_HepG2$Protein.ID,
  OGlcNAc_protein_DE_Jurkat$Protein.ID
))

cat("Overlapping proteins:", length(overlapping), "\n")

logFC_mat <- tibble(Protein.ID = overlapping) |>
  left_join(OGlcNAc_protein_DE_HEK293T |> dplyr::select(Protein.ID, logFC_HEK293T = logFC), by = "Protein.ID") |>
  left_join(OGlcNAc_protein_DE_HepG2 |> dplyr::select(Protein.ID, logFC_HepG2 = logFC), by = "Protein.ID") |>
  left_join(OGlcNAc_protein_DE_Jurkat |> dplyr::select(Protein.ID, logFC_Jurkat = logFC), by = "Protein.ID")

# Re-derive commonly regulated
get_sig <- function(de, dir) {
  if (dir == "up") de |> filter(logFC > 0.5, adj.P.Val < 0.05) |> pull(Protein.ID)
  else de |> filter(logFC < -0.5, adj.P.Val < 0.05) |> pull(Protein.ID)
}
up_H <- get_sig(OGlcNAc_protein_DE_HEK293T, "up")
up_G <- get_sig(OGlcNAc_protein_DE_HepG2, "up")
up_J <- get_sig(OGlcNAc_protein_DE_Jurkat, "up")
down_H <- get_sig(OGlcNAc_protein_DE_HEK293T, "down")
down_G <- get_sig(OGlcNAc_protein_DE_HepG2, "down")
down_J <- get_sig(OGlcNAc_protein_DE_Jurkat, "down")

commonly_up <- setdiff(unique(c(intersect(up_H, up_G), intersect(up_H, up_J), intersect(up_G, up_J))), c(down_H, down_G, down_J))
commonly_down <- setdiff(unique(c(intersect(down_H, down_G), intersect(down_H, down_J), intersect(down_G, down_J))), c(up_H, up_G, up_J))

# Run UMAP
set.seed(42)
umap_result <- umap(as.matrix(logFC_mat |> dplyr::select(starts_with("logFC_"))),
                     n_neighbors = 15, min_dist = 0.1, n_components = 2)

umap_df <- tibble(
  Protein.ID = logFC_mat$Protein.ID,
  UMAP1 = umap_result[,1],
  UMAP2 = umap_result[,2],
  Regulation = case_when(
    Protein.ID %in% commonly_up ~ "Generally Up",
    Protein.ID %in% commonly_down ~ "Generally Down",
    TRUE ~ "Other"
  )
) |> mutate(Regulation = factor(Regulation, levels = c("Generally Up", "Generally Down", "Other")))

cat("UMAP classification:", table(umap_df$Regulation), "\n")

Figure3C <- ggplot() +
  geom_point(data = umap_df |> filter(Regulation == "Other"),
             aes(x = UMAP1, y = UMAP2, color = Regulation, shape = Regulation), size = 1, alpha = 0.5) +
  geom_point(data = umap_df |> filter(Regulation != "Other"),
             aes(x = UMAP1, y = UMAP2, color = Regulation, shape = Regulation), size = 2, alpha = 0.8) +
  scale_color_manual(values = c("Generally Up" = "#F39B7F", "Generally Down" = "#4DBBD5", "Other" = "grey70"),
                     breaks = c("Generally Up", "Generally Down")) +
  scale_shape_manual(values = c("Generally Up" = 17, "Generally Down" = 15, "Other" = 16),
                     breaks = c("Generally Up", "Generally Down")) +
  labs(x = "UMAP1", y = "UMAP2", color = "", shape = "") +
  theme_bw() +
  theme(panel.grid.major = element_line(linewidth = 0.2, color = "gray"),
        panel.grid.minor = element_line(linewidth = 0.1, color = "gray"),
        axis.title = element_text(size = 9), axis.text = element_text(size = 9, color = "black"),
        legend.position = "bottom", legend.text = element_text(size = 7), legend.key.size = unit(0.2, "cm"))

ggsave(paste0(output_path, "Figure3C.pdf"), Figure3C, width = 1.8, height = 2.2)
cat("Figure 3C saved\n")

cat("\n=== All subfigures generated ===\n")
cat("Files in output:\n")
print(list.files(output_path))
