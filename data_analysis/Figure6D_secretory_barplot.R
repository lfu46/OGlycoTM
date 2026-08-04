library(tidyverse)

source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'
figure_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/Figure6_OGalNAc/'

# Load normalized Jurkat OGalNAc protein data
norm <- read_csv(paste0(source_file_path, 'normalization/OGalNAc_protein_norm_Jurkat.csv'),
                 show_col_types = FALSE)

# Compute per-replicate log2FC
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

# Define proteins and categories (order within each category by logFC)
protein_category <- tribble(
  ~Gene, ~Category,
  "TFRC",    "Cell-surface receptors",
  "PTPRC",   "Cell-surface receptors",
  "IGSF8",   "Cell-surface receptors",
  "CPD",     "ER / Golgi processing",
  "MAN1A2",  "ER / Golgi processing",
  "XXYLT1",  "ER / Golgi processing",
  "OS9",     "ER quality control",
  "EDEM1",   "ER quality control",
  "SEL1L",   "ER quality control",
)

category_order <- c("ER quality control", "ER / Golgi processing", "Cell-surface receptors")

# Filter and join
plot_data <- fc_data |>
  inner_join(protein_category, by = "Gene")

mean_data <- plot_data |>
  group_by(Gene, Category) |>
  summarize(mean_log2FC = mean(log2FC),
            sd_log2FC = sd(log2FC), .groups = "drop")

# Order genes within each category by mean logFC
gene_order <- mean_data |>
  mutate(Category = factor(Category, levels = category_order)) |>
  arrange(Category, mean_log2FC) |>
  pull(Gene)

plot_data <- plot_data |>
  mutate(
    Gene = factor(Gene, levels = gene_order),
    Category = factor(Category, levels = category_order)
  )

mean_data <- mean_data |>
  mutate(
    Gene = factor(Gene, levels = gene_order),
    Category = factor(Category, levels = category_order)
  )

# Build gradient tiles for each bar using geom_tile (facet-compatible)
n_slices <- 50
grad_data <- mean_data |>
  rowwise() |>
  do({
    row <- .
    y_centers <- seq(row$mean_log2FC, 0, length.out = n_slices + 1)
    y_centers <- (y_centers[1:n_slices] + y_centers[2:(n_slices + 1)]) / 2
    slice_h <- abs(row$mean_log2FC) / n_slices
    frac <- seq(0, 1, length.out = n_slices)
    tibble(Gene = row$Gene, Category = row$Category,
           y = y_centers, h = slice_h, frac = frac)
  }) |>
  ungroup() |>
  mutate(
    Gene = factor(Gene, levels = levels(mean_data$Gene)),
    Category = factor(Category, levels = levels(mean_data$Category))
  )

p <- ggplot() +
  geom_tile(data = grad_data,
            aes(x = Gene, y = y, height = h, fill = frac),
            width = 0.6, show.legend = FALSE) +
  geom_hline(yintercept = 0, color = "grey40", linewidth = 0.5) +
  geom_errorbar(data = mean_data,
                aes(x = Gene,
                    ymin = mean_log2FC - sd_log2FC,
                    ymax = mean_log2FC + sd_log2FC),
                width = 0.25, linewidth = 0.4, color = "black") +
  scale_fill_gradient(low = "#004D40", high = "#80CBC4") +
  facet_grid(. ~ Category, scales = "free_x", space = "free_x") +
  labs(x = NULL, y = expression(log[2]*"(Tuni/Ctrl)")) +
  theme_classic(base_size = 7) +
  theme(
    axis.text.x = element_text(color = "black", size = 7, angle = 45, hjust = 1),
    axis.text.y = element_text(color = "black", size = 6),
    axis.title.y = element_text(size = 7),
    strip.text = element_text(size = 5.5, face = "bold"),
    strip.background = element_blank(),
    panel.spacing = unit(0.3, "lines")
  )

ggsave(paste0(figure_file_path, "Figure6D_secretory_barplot.pdf"),
       p, width = 3.2, height = 1.5)
ggsave(paste0(figure_file_path, "Figure6D_secretory_barplot.png"),
       p, width = 3.2, height = 1.5, dpi = 300)
cat("Saved Figure6D\n")
