# Statistical power analysis for the OGlycoTM TMT-limma differential design
#
# For each (level, cell type, log2FC effect size, true-positive fraction):
#   1. Take per-feature mean and pooled within-group SD from the normalized matrix.
#   2. Simulate N replicate datasets with those mean/SD parameters.
#   3. For a random subset (true positives), add the chosen log2FC to the Tuni columns.
#   4. Fit the same limma pipeline used in differential_analysis.R.
#   5. Score detection power (TP detected / total TP) and empirical FDR (FP / called).
#
# Output: power_analysis_table.csv + power_analysis_plot.pdf/png

suppressMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
  library(limma)
  library(ggplot2)
})

set.seed(42)

# ---- Configuration ----------------------------------------------------------

NORM_DIR <- "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/normalization"
OUT_DIR  <- "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/Power_analysis"
LOCAL_OUT <- "/tmp/power_analysis"
dir.create(LOCAL_OUT, showWarnings = FALSE, recursive = TRUE)
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

INTENSITY_COLS <- c(
  "Intensity.Tuni_1_sl_tmm", "Intensity.Tuni_2_sl_tmm", "Intensity.Tuni_3_sl_tmm",
  "Intensity.Ctrl_4_sl_tmm", "Intensity.Ctrl_5_sl_tmm", "Intensity.Ctrl_6_sl_tmm"
)

LEVELS <- list(
  list(level = "OGlcNAc_protein", template = "OGlcNAc_protein_norm_%s.csv"),
  list(level = "OGlcNAc_site",    template = "OGlcNAc_site_norm_%s.csv"),
  list(level = "OGalNAc_protein", template = "OGalNAc_protein_norm_%s.csv"),
  list(level = "OGalNAc_site",    template = "OGalNAc_site_norm_%s.csv"),
  list(level = "WP_protein",      template = "WP_protein_norm_%s.csv")  # reference
)
CELLS         <- c("HEK293T", "HepG2", "Jurkat")
EFFECT_SIZES  <- c(0.5, 0.75, 1.0, 1.5)
TP_FRACTIONS  <- c(0.10, 0.20)
N_ITER        <- 500
ALPHA         <- 0.05
LOGFC_CUTOFF  <- 0.5

# ---- Limma design (matches differential_analysis.R) -------------------------

design <- model.matrix(~ 0 + factor(rep(c("Tuni", "Ctrl"), each = 3),
                                    levels = c("Tuni", "Ctrl")))
colnames(design) <- c("Tuni", "Ctrl")
contrast <- makeContrasts(Tuni_vs_Ctrl = Tuni - Ctrl, levels = design)

# ---- Helpers ----------------------------------------------------------------

load_matrix <- function(path) {
  d <- suppressMessages(read_csv(path, show_col_types = FALSE))
  M <- as.matrix(d[, INTENSITY_COLS])
  keep <- complete.cases(M) & rowSums(M > 0) == 6
  Mlog <- log2(M[keep, , drop = FALSE])
  list(
    Mlog   = Mlog,
    mu     = rowMeans(Mlog),
    sd_pool = {
      sd_v <- sqrt((apply(Mlog[, 1:3], 1, var) + apply(Mlog[, 4:6], 1, var)) / 2)
      med  <- median(sd_v, na.rm = TRUE)
      sd_v[is.na(sd_v) | sd_v <= 0] <- med
      sd_v
    }
  )
}

simulate_iter <- function(mu, sd_pool, delta, tp_frac) {
  n <- length(mu)
  n_tp <- max(1L, round(n * tp_frac))
  tp_idx <- sample.int(n, n_tp)

  Y <- matrix(rnorm(n * 6, mean = rep(mu, 6), sd = rep(sd_pool, 6)),
              nrow = n, ncol = 6)
  Y[tp_idx, 1:3] <- Y[tp_idx, 1:3] + delta

  fit  <- lmFit(Y, design)
  fit2 <- contrasts.fit(fit, contrast)
  fit2 <- eBayes(fit2)
  tt   <- topTable(fit2, number = n, sort.by = "none")

  hit_q  <- tt$adj.P.Val < ALPHA
  hit_qf <- hit_q & abs(tt$logFC) > LOGFC_CUTOFF

  is_tp <- logical(n); is_tp[tp_idx] <- TRUE

  c(
    power_q       = sum(hit_q  &  is_tp) / n_tp,
    power_qf      = sum(hit_qf &  is_tp) / n_tp,
    fdr_q         = if (sum(hit_q)  > 0) sum(hit_q  & !is_tp) / sum(hit_q)  else 0,
    fdr_qf        = if (sum(hit_qf) > 0) sum(hit_qf & !is_tp) / sum(hit_qf) else 0,
    n_called_q    = sum(hit_q),
    n_called_qf   = sum(hit_qf)
  )
}

# ---- Main loop --------------------------------------------------------------

results <- list()

t0 <- Sys.time()
for (lvl in LEVELS) {
  for (cell in CELLS) {
    path <- file.path(NORM_DIR, sprintf(lvl$template, cell))
    if (!file.exists(path)) {
      message("Skipping (no file): ", path); next
    }
    d <- load_matrix(path)
    n_features <- nrow(d$Mlog)
    if (n_features < 20) {
      message("Skipping (n<20): ", path); next
    }

    for (es in EFFECT_SIZES) {
      for (tpf in TP_FRACTIONS) {
        sims <- replicate(N_ITER,
                          simulate_iter(d$mu, d$sd_pool, es, tpf))
        # sims is 6 x N_ITER
        results[[length(results) + 1]] <- data.frame(
          Level         = lvl$level,
          CellType      = cell,
          log2FC        = es,
          TP_fraction   = tpf,
          N_features    = n_features,
          N_TP_per_iter = max(1L, round(n_features * tpf)),
          Power_q       = mean(sims["power_q",  ]),
          Power_qf      = mean(sims["power_qf", ]),
          Empirical_FDR_q  = mean(sims["fdr_q",  ]),
          Empirical_FDR_qf = mean(sims["fdr_qf", ]),
          Mean_called_q   = mean(sims["n_called_q",  ]),
          Mean_called_qf  = mean(sims["n_called_qf", ])
        )
        message(sprintf("[%5.1fs] %-18s %-7s log2FC=%.2f TPfrac=%.2f -> Power_q=%.2f FDR_q=%.3f",
                        as.numeric(difftime(Sys.time(), t0, units = "secs")),
                        lvl$level, cell, es, tpf,
                        mean(sims["power_q", ]), mean(sims["fdr_q", ])))
      }
    }
  }
}

power_table <- bind_rows(results) %>%
  arrange(Level, CellType, TP_fraction, log2FC)

# ---- Save outputs (local first, then network) ------------------------------

local_csv <- file.path(LOCAL_OUT, "power_analysis_table.csv")
write_csv(power_table, local_csv)

# Plot 1: primary (TP=10%, q-only threshold)
p1 <- power_table %>%
  filter(TP_fraction == 0.10) %>%
  ggplot(aes(x = log2FC, y = Power_q, color = Level, group = Level)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ CellType) +
  scale_y_continuous(limits = c(0, 1.02), labels = scales::percent_format(accuracy = 1)) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 11) +
  theme(text = element_text(family = "Arial"),
        legend.position = "bottom",
        legend.title = element_blank()) +
  labs(x = expression(paste("Simulated log"[2], "FC")),
       y = "Estimated power (adj.P.Val < 0.05)",
       title = "Power analysis: triplicate TMT-limma design",
       subtitle = "True-positive fraction = 10%; 500 simulations per condition")

# Plot 2: TP=10% vs 20% sensitivity (linetype)
p2 <- power_table %>%
  ggplot(aes(x = log2FC, y = Power_q, color = Level,
             group = interaction(Level, TP_fraction),
             linetype = factor(TP_fraction))) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  facet_wrap(~ CellType) +
  scale_y_continuous(limits = c(0, 1.02), labels = scales::percent_format(accuracy = 1)) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 11) +
  theme(text = element_text(family = "Arial"),
        legend.position = "bottom") +
  labs(x = expression(paste("Simulated log"[2], "FC")),
       y = "Estimated power (adj.P.Val < 0.05)",
       linetype = "True-positive fraction",
       color = "Feature level",
       title = "Power analysis: sensitivity to true-positive fraction")

# Plot 3: stricter threshold (q + |logFC|>0.5), TP=10%
p3 <- power_table %>%
  filter(TP_fraction == 0.10) %>%
  ggplot(aes(x = log2FC, y = Power_qf, color = Level, group = Level)) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  facet_wrap(~ CellType) +
  scale_y_continuous(limits = c(0, 1.02), labels = scales::percent_format(accuracy = 1)) +
  scale_color_brewer(palette = "Set1") +
  theme_bw(base_size = 11) +
  theme(text = element_text(family = "Arial"),
        legend.position = "bottom",
        legend.title = element_blank()) +
  labs(x = expression(paste("Simulated log"[2], "FC")),
       y = expression(paste("Estimated power (adj.P.Val < 0.05  &  |log"[2], "FC| > 0.5)")),
       title = "Power analysis with combined threshold",
       subtitle = "True-positive fraction = 10%")

ggsave(file.path(LOCAL_OUT, "power_plot_primary.pdf"), p1, width = 9, height = 4.5)
ggsave(file.path(LOCAL_OUT, "power_plot_primary.png"), p1, width = 9, height = 4.5, dpi = 600)
ggsave(file.path(LOCAL_OUT, "power_plot_TPsensitivity.pdf"), p2, width = 9, height = 4.5)
ggsave(file.path(LOCAL_OUT, "power_plot_TPsensitivity.png"), p2, width = 9, height = 4.5, dpi = 600)
ggsave(file.path(LOCAL_OUT, "power_plot_combined_threshold.pdf"), p3, width = 9, height = 4.5)
ggsave(file.path(LOCAL_OUT, "power_plot_combined_threshold.png"), p3, width = 9, height = 4.5, dpi = 600)

# Sync to network drive
file.copy(list.files(LOCAL_OUT, full.names = TRUE),
          OUT_DIR, overwrite = TRUE)

cat("\n=== Power analysis complete ===\n")
cat("Output dir:", OUT_DIR, "\n")
cat("Total runtime:", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 2), "min\n")
print(power_table %>% filter(TP_fraction == 0.10) %>%
        select(Level, CellType, log2FC, N_features, Power_q, Empirical_FDR_q),
      n = Inf)
