# Figure: O-GlcNAc sites in ER/Golgi/PM-only proteins — N-glyc sequon analysis
library(tidyverse)
library(devEMF)

figure_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/'
colors_cell <- c("HEK293T" = "#4DBBD5", "HepG2" = "#F39B7F", "Jurkat" = "#00A087")

# Data from analysis (unique site level)
df <- tibble(
  Cell = factor(c("HEK293T", "HepG2", "Jurkat"), levels = c("HEK293T", "HepG2", "Jurkat")),
  ER = c(7, 5, 4),
  Golgi = c(1, 5, 6),
  PM = c(1, 0, 3),
  total_sites = c(9, 10, 13),
  n_sequon = c(1, 4, 0),
  n_proteins = c(7, 9, 13)
)
df <- df |> mutate(n_no_sequon = total_sites - n_sequon)

# --- Panel A: Stacked bar by subcellular location ---
loc_long <- df |>
  pivot_longer(cols = c(ER, Golgi, PM), names_to = "Location", values_to = "n_sites") |>
  mutate(Location = factor(Location, levels = c("PM", "Golgi", "ER")))

colors_loc <- c("ER" = "#E64B35", "Golgi" = "#3C5488", "PM" = "#8491B4")

pA <- ggplot(loc_long, aes(x = Cell, y = n_sites, fill = Location)) +
  geom_bar(stat = "identity", width = 0.65, color = "black", linewidth = 0.3) +
  geom_text(aes(label = after_stat(y), group = Cell),
            stat = "summary", fun = sum, vjust = -0.5, size = 3.5, fontface = "bold") +
  scale_fill_manual(values = colors_loc, name = "Subcellular\nlocation") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Number of O-GlcNAc sites",
       title = "O-GlcNAc sites in ER/Golgi/PM-only proteins") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )

# --- Panel B: Sequon stacked bar ---
seq_long <- df |>
  dplyr::select(Cell, `With N-sequon` = n_sequon, `Without N-sequon` = n_no_sequon) |>
  pivot_longer(-Cell, names_to = "Sequon", values_to = "n") |>
  mutate(Sequon = factor(Sequon, levels = c("Without N-sequon", "With N-sequon")))

colors_seq <- c("With N-sequon" = "#E64B35", "Without N-sequon" = "#4DBBD5")

pB <- ggplot(seq_long, aes(x = Cell, y = n, fill = Sequon)) +
  geom_bar(stat = "identity", width = 0.65, color = "black", linewidth = 0.3) +
  geom_text(
    data = df |> filter(n_sequon > 0),
    aes(x = Cell, y = total_sites, fill = NULL,
        label = paste0(round(100 * n_sequon / total_sites, 1), "%")),
    vjust = -0.5, size = 3.5, fontface = "bold"
  ) +
  geom_text(
    data = df |> filter(n_sequon == 0),
    aes(x = Cell, y = total_sites, fill = NULL, label = "0%"),
    vjust = -0.5, size = 3.5, fontface = "bold", color = "grey40"
  ) +
  scale_fill_manual(values = colors_seq, name = "N-X-S/T\nsequon") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(x = NULL, y = "Number of O-GlcNAc sites",
       title = "N-glycosylation sequon in ER/Golgi/PM sites") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 11, face = "bold"),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11),
    legend.position = "right",
    legend.text = element_text(size = 10),
    legend.title = element_text(size = 10, face = "bold"),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5)
  )

# --- Combined figure ---
library(patchwork)
combined <- pA + pB + plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))

# Save PDF
out_dir <- paste0(figure_file_path, "OGlcNAc_ER_Golgi_PM_spectra/")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

ggsave(paste0(out_dir, "OGlcNAc_ER_Golgi_PM_sequon_analysis.pdf"),
       combined, width = 9, height = 4, dpi = 300)
cat("Saved PDF\n")

# Save EMF
emf(paste0(out_dir, "OGlcNAc_ER_Golgi_PM_sequon_analysis.emf"),
    width = 9, height = 4)
print(combined)
dev.off()
cat("Saved EMF\n")
