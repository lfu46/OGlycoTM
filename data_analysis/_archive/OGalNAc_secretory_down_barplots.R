library(tidyverse)

source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'
figure_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/Figure6_OGalNAc/'

# Load normalized Jurkat OGalNAc protein data
norm <- read_csv(paste0(source_file_path, 'normalization/OGalNAc_protein_norm_Jurkat.csv'),
                 show_col_types = FALSE)

# Compute per-replicate log2FC (Tuni/mean_Ctrl)
fc_data <- norm |>
  rowwise() |>
  mutate(
    mean_Ctrl = mean(c(Intensity.Ctrl_4_sl_tmm, Intensity.Ctrl_5_sl_tmm, Intensity.Ctrl_6_sl_tmm)),
    log2FC_rep1 = log2(Intensity.Tuni_1_sl_tmm / mean_Ctrl),
    log2FC_rep2 = log2(Intensity.Tuni_2_sl_tmm / mean_Ctrl),
    log2FC_rep3 = log2(Intensity.Tuni_3_sl_tmm / mean_Ctrl)
  ) |>
  ungroup() |>
  select(Protein.ID, Gene, log2FC_rep1, log2FC_rep2, log2FC_rep3) |>
  pivot_longer(cols = starts_with("log2FC"), names_to = "replicate", values_to = "log2FC")

# Define categories
categories <- list(
  "ER quality control" = c("SEL1L", "PRKCSH", "TXNDC5", "EDEM1", "OS9"),
  "ER-Golgi transport" = c("TMEM165", "MIA3"),
  "Glycosyltransferases" = c("MANEA", "MAN1A2", "GALNT2", "B4GALT4", "B3GALT6", "XXYLT1"),
  "Cell surface receptors" = c("TFRC", "IGSF8", "CD99", "PTPRC", "LDLR", "SPN"),
  "Transporters" = c("SLC29A1", "SLC39A10", "GPR108")
)

# Plot function
make_barplot <- function(genes, title, width, height) {
  plot_data <- fc_data |>
    filter(Gene %in% genes) |>
    mutate(Gene = factor(Gene, levels = genes))

  mean_data <- plot_data |>
    group_by(Gene) |>
    summarize(mean_log2FC = mean(log2FC), .groups = "drop")

  p <- ggplot() +
    geom_hline(yintercept = 0, color = "grey40", linewidth = 0.4) +
    geom_col(data = mean_data, aes(x = Gene, y = mean_log2FC),
             fill = "#4DBBD5", alpha = 0.6, width = 0.6) +
    geom_point(data = plot_data, aes(x = Gene, y = log2FC),
               size = 1.5, color = "black",
               position = position_jitter(width = 0.1, seed = 42)) +
    labs(x = NULL, y = expression(log[2]*"(Tuni/Ctrl)"), title = title) +
    theme_classic(base_size = 9) +
    theme(
      plot.title = element_text(size = 9, face = "bold", hjust = 0.5),
      axis.text.x = element_text(color = "black", size = 8, angle = 45, hjust = 1, face = "italic"),
      axis.text.y = element_text(color = "black", size = 8),
      axis.title.y = element_text(size = 9),
      plot.margin = margin(5, 10, 5, 10)
    )

  ggsave(paste0(figure_file_path, "OGalNAc_secretory_down_",
                str_replace_all(title, " ", "_"), ".pdf"),
         p, width = width, height = height)
  ggsave(paste0(figure_file_path, "OGalNAc_secretory_down_",
                str_replace_all(title, " ", "_"), ".png"),
         p, width = width, height = height, dpi = 300)
  cat("Saved:", title, "\n")
  p
}

# Generate plots — order genes by mean logFC (most downregulated first)
for (cat_name in names(categories)) {
  genes <- categories[[cat_name]]
  # Sort by mean logFC
  gene_means <- fc_data |>
    filter(Gene %in% genes) |>
    group_by(Gene) |>
    summarize(m = mean(log2FC), .groups = "drop") |>
    arrange(m)
  sorted_genes <- gene_means$Gene

  n <- length(sorted_genes)
  w <- max(1.5, 0.6 * n + 0.8)
  h <- 2.2

  make_barplot(sorted_genes, cat_name, w, h)
}

cat("\nAll barplots saved to:", figure_file_path, "\n")
