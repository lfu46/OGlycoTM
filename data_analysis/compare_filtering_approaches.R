# compare_filtering_approaches.R
# Compare two filtering approaches for O-GlcNAc site quantification:
#   Option 1 (Site-level): Filter sites after quantification (current approach)
#   Option 2 (PSM-level): Filter PSMs before quantification (stricter approach)
#
# This script generates a fair comparison of the two approaches.

library(tidyverse)
library(limma)
library(edgeR)

# =============================================================================
# Configuration
# =============================================================================

source_file_path <- "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/"

# Helper function to parse site probability
parse_site_probability <- function(prob_str) {
  prob <- str_extract(prob_str, "[0-9.]+(?=\\]$)")
  as.numeric(prob)
}

# =============================================================================
# Load Raw Site Data (before any filtering)
# =============================================================================

cat("=== Loading Raw Site Data ===\n\n")

site_files <- list(
  HEK293T = paste0(source_file_path, "site/OGlcNAc_site_HEK293T.csv"),
  HepG2 = paste0(source_file_path, "site/OGlcNAc_site_HepG2.csv"),
  Jurkat = paste0(source_file_path, "site/OGlcNAc_site_Jurkat.csv")
)

raw_site_data <- lapply(site_files, read_csv, show_col_types = FALSE)

for (cell_type in names(raw_site_data)) {
  cat(cell_type, ": ", nrow(raw_site_data[[cell_type]]), " PSMs\n", sep = "")
}

# =============================================================================
# OPTION 1: Site-Level Filtering (Current Approach)
# Quantify ALL PSMs, then filter sites that have at least one high-prob PSM
# =============================================================================

cat("\n=== OPTION 1: Site-Level Filtering ===\n")
cat("Quantify all PSMs, then keep sites with any prob >= 0.75 PSM\n\n")

option1_results <- list()

for (cell_type in names(raw_site_data)) {
  data <- raw_site_data[[cell_type]]

  # Parse probability
  data <- data %>%
    mutate(Site_Prob = parse_site_probability(Site.Probabilities))

  # Quantify ALL PSMs (aggregate by site_index)
  site_quant <- data %>%
    group_by(site_index) %>%
    summarize(
      Protein.ID = first(Protein.ID),
      Entry.Name = first(Entry.Name),
      Gene = first(Gene),
      Protein.Description = first(Protein.Description),
      modified_residue = first(modified_residue),
      site_number = first(site_number),
      n_psms_total = n(),
      Intensity.Tuni_1 = sum(Intensity.Tuni_1),
      Intensity.Tuni_2 = sum(Intensity.Tuni_2),
      Intensity.Tuni_3 = sum(Intensity.Tuni_3),
      Intensity.Ctrl_4 = sum(Intensity.Ctrl_4),
      Intensity.Ctrl_5 = sum(Intensity.Ctrl_5),
      Intensity.Ctrl_6 = sum(Intensity.Ctrl_6),
      .groups = "drop"
    )

  # Find sites that have at least one high-confidence PSM
  high_conf_sites <- data %>%
    filter(
      Confidence.Level == "Level1" |
      (Confidence.Level == "Level1b" & Site_Prob >= 0.75)
    ) %>%
    pull(site_index) %>%
    unique()

  # Filter quantified data to keep only high-confidence sites
  site_quant_filtered <- site_quant %>%
    filter(site_index %in% high_conf_sites)

  option1_results[[cell_type]] <- site_quant_filtered

  cat(cell_type, ":\n")
  cat("  Total sites quantified: ", nrow(site_quant), "\n")
  cat("  Sites with high-conf PSM: ", length(high_conf_sites), "\n")
  cat("  Sites retained: ", nrow(site_quant_filtered), "\n\n")
}

# =============================================================================
# OPTION 2: PSM-Level Filtering (Stricter Approach)
# Filter PSMs first, then quantify only high-confidence PSMs
# =============================================================================

cat("=== OPTION 2: PSM-Level Filtering ===\n")
cat("Filter PSMs first (prob >= 0.75), then quantify\n\n")

option2_results <- list()

for (cell_type in names(raw_site_data)) {
  data <- raw_site_data[[cell_type]]

  # Parse probability
  data <- data %>%
    mutate(Site_Prob = parse_site_probability(Site.Probabilities))

  # Count PSMs before filtering
  n_before <- nrow(data)
  n_level1 <- sum(data$Confidence.Level == "Level1", na.rm = TRUE)
  n_level1b <- sum(data$Confidence.Level == "Level1b", na.rm = TRUE)
  n_level1b_low <- sum(data$Confidence.Level == "Level1b" & data$Site_Prob < 0.75, na.rm = TRUE)

  # Filter PSMs FIRST
  data_filtered <- data %>%
    filter(
      Confidence.Level == "Level1" |
      (Confidence.Level == "Level1b" & Site_Prob >= 0.75)
    )

  n_after <- nrow(data_filtered)

  # Quantify only high-confidence PSMs
  site_quant <- data_filtered %>%
    group_by(site_index) %>%
    summarize(
      Protein.ID = first(Protein.ID),
      Entry.Name = first(Entry.Name),
      Gene = first(Gene),
      Protein.Description = first(Protein.Description),
      modified_residue = first(modified_residue),
      site_number = first(site_number),
      n_psms_highconf = n(),
      Intensity.Tuni_1 = sum(Intensity.Tuni_1),
      Intensity.Tuni_2 = sum(Intensity.Tuni_2),
      Intensity.Tuni_3 = sum(Intensity.Tuni_3),
      Intensity.Ctrl_4 = sum(Intensity.Ctrl_4),
      Intensity.Ctrl_5 = sum(Intensity.Ctrl_5),
      Intensity.Ctrl_6 = sum(Intensity.Ctrl_6),
      .groups = "drop"
    )

  option2_results[[cell_type]] <- site_quant

  cat(cell_type, ":\n")
  cat("  PSMs before filtering: ", n_before, "\n")
  cat("    Level1: ", n_level1, "\n")
  cat("    Level1b: ", n_level1b, " (", n_level1b_low, " with prob < 0.75)\n")
  cat("  PSMs after filtering: ", n_after, "\n")
  cat("  Sites quantified: ", nrow(site_quant), "\n\n")
}

# =============================================================================
# Normalization Function (SL + TMM)
# =============================================================================

normalize_data <- function(data, intensity_cols) {
  # Sample loading normalization
  intensity_matrix <- data %>%
    select(all_of(intensity_cols)) %>%
    as.matrix()

  # Replace zeros with small value to avoid issues
  intensity_matrix[intensity_matrix == 0] <- 1

  # Calculate column sums
  col_sums <- colSums(intensity_matrix, na.rm = TRUE)
  target_sum <- mean(col_sums)

  # SL normalization factors
  sl_factors <- target_sum / col_sums

  # Apply SL normalization
  sl_normalized <- sweep(intensity_matrix, 2, sl_factors, "*")

  # TMM normalization using edgeR
  # Transpose so samples are columns
  dge <- DGEList(counts = sl_normalized)
  dge <- calcNormFactors(dge, method = "TMM")
  tmm_factors <- dge$samples$norm.factors

  # Apply TMM normalization
  tmm_normalized <- sweep(sl_normalized, 2, tmm_factors, "/")

  # Create normalized column names
  norm_col_names <- paste0(gsub("Intensity\\.", "Intensity.", intensity_cols), "_sl_tmm")

  # Add normalized columns to data
  result <- data
  for (i in seq_along(intensity_cols)) {
    result[[norm_col_names[i]]] <- tmm_normalized[, i]
  }

  return(result)
}

# =============================================================================
# Differential Analysis Function
# =============================================================================

run_limma_DE <- function(data, id_col, intensity_cols) {
  # Design matrix
  design <- model.matrix(~ 0 + factor(rep(c("Tuni", "Ctrl"), each = 3),
                                       levels = c("Tuni", "Ctrl")))
  colnames(design) <- c("Tuni", "Ctrl")

  # Contrast
  contrast <- makeContrasts(Tuni_vs_Ctrl = Tuni - Ctrl, levels = design)

  # Log2 transform
  log2_data <- data %>%
    mutate(across(all_of(intensity_cols), ~ log2(.x)))

  # Create matrix
  data_matrix <- log2_data %>%
    select(all_of(intensity_cols)) %>%
    as.matrix()
  rownames(data_matrix) <- data[[id_col]]

  # Fit model
  fit <- lmFit(data_matrix, design)
  fit <- contrasts.fit(fit, contrast)
  fit <- eBayes(fit)

  # Extract results
  top_table <- topTable(fit, number = Inf, adjust.method = "BH")

  result <- as_tibble(top_table)
  result[[id_col]] <- rownames(top_table)

  return(result)
}

# =============================================================================
# Run Normalization and DE Analysis for Both Options
# =============================================================================

cat("=== Running Normalization and Differential Analysis ===\n\n")

intensity_cols <- c("Intensity.Tuni_1", "Intensity.Tuni_2", "Intensity.Tuni_3",
                    "Intensity.Ctrl_4", "Intensity.Ctrl_5", "Intensity.Ctrl_6")
norm_intensity_cols <- paste0(intensity_cols, "_sl_tmm")

option1_de <- list()
option2_de <- list()

for (cell_type in names(raw_site_data)) {
  cat("Processing", cell_type, "...\n")

  # Option 1
  data1 <- normalize_data(option1_results[[cell_type]], intensity_cols)
  de1 <- run_limma_DE(data1, "site_index", norm_intensity_cols)
  de1$cell_type <- cell_type
  option1_de[[cell_type]] <- de1

  # Option 2
  data2 <- normalize_data(option2_results[[cell_type]], intensity_cols)
  de2 <- run_limma_DE(data2, "site_index", norm_intensity_cols)
  de2$cell_type <- cell_type
  option2_de[[cell_type]] <- de2
}

# =============================================================================
# Compare Results
# =============================================================================

cat("\n=== COMPARISON: Option 1 vs Option 2 ===\n\n")

comparison_stats <- data.frame(
  cell_type = character(),
  n_sites_option1 = integer(),
  n_sites_option2 = integer(),
  n_common = integer(),
  n_only_option1 = integer(),
  n_only_option2 = integer(),
  correlation_logFC = numeric(),
  mean_abs_diff_logFC = numeric(),
  max_abs_diff_logFC = numeric(),
  stringsAsFactors = FALSE
)

all_comparisons <- list()

for (cell_type in names(raw_site_data)) {
  de1 <- option1_de[[cell_type]]
  de2 <- option2_de[[cell_type]]

  # Find common sites
  common_sites <- intersect(de1$site_index, de2$site_index)
  only_option1 <- setdiff(de1$site_index, de2$site_index)
  only_option2 <- setdiff(de2$site_index, de1$site_index)

  # Merge for comparison
  merged <- de1 %>%
    filter(site_index %in% common_sites) %>%
    select(site_index, logFC_opt1 = logFC, adj.P.Val_opt1 = adj.P.Val) %>%
    left_join(
      de2 %>%
        filter(site_index %in% common_sites) %>%
        select(site_index, logFC_opt2 = logFC, adj.P.Val_opt2 = adj.P.Val),
      by = "site_index"
    ) %>%
    mutate(
      logFC_diff = logFC_opt1 - logFC_opt2,
      cell_type = cell_type
    )

  all_comparisons[[cell_type]] <- merged

  # Calculate statistics
  corr <- cor(merged$logFC_opt1, merged$logFC_opt2, use = "complete.obs")
  mean_diff <- mean(abs(merged$logFC_diff), na.rm = TRUE)
  max_diff <- max(abs(merged$logFC_diff), na.rm = TRUE)

  comparison_stats <- rbind(comparison_stats, data.frame(
    cell_type = cell_type,
    n_sites_option1 = nrow(de1),
    n_sites_option2 = nrow(de2),
    n_common = length(common_sites),
    n_only_option1 = length(only_option1),
    n_only_option2 = length(only_option2),
    correlation_logFC = round(corr, 4),
    mean_abs_diff_logFC = round(mean_diff, 4),
    max_abs_diff_logFC = round(max_diff, 4)
  ))

  cat(cell_type, ":\n")
  cat("  Sites in Option 1: ", nrow(de1), "\n")
  cat("  Sites in Option 2: ", nrow(de2), "\n")
  cat("  Common sites: ", length(common_sites), "\n")
  cat("  Only in Option 1: ", length(only_option1), "\n")
  cat("  Only in Option 2: ", length(only_option2), "\n")
  cat("  logFC correlation: ", round(corr, 4), "\n")
  cat("  Mean |logFC difference|: ", round(mean_diff, 4), "\n")
  cat("  Max |logFC difference|: ", round(max_diff, 4), "\n\n")
}

# Combine all comparisons
all_comparisons_df <- bind_rows(all_comparisons)

# =============================================================================
# Summary Statistics
# =============================================================================

cat("=== OVERALL SUMMARY ===\n\n")

# Overall correlation
overall_corr <- cor(all_comparisons_df$logFC_opt1, all_comparisons_df$logFC_opt2, use = "complete.obs")
overall_mean_diff <- mean(abs(all_comparisons_df$logFC_diff), na.rm = TRUE)
overall_max_diff <- max(abs(all_comparisons_df$logFC_diff), na.rm = TRUE)

cat("Total sites compared: ", nrow(all_comparisons_df), "\n")
cat("Overall logFC correlation: ", round(overall_corr, 4), "\n")
cat("Overall mean |logFC difference|: ", round(overall_mean_diff, 4), "\n")
cat("Overall max |logFC difference|: ", round(overall_max_diff, 4), "\n\n")

# Distribution of differences
cat("Distribution of logFC differences:\n")
quantiles <- quantile(all_comparisons_df$logFC_diff, probs = c(0, 0.05, 0.25, 0.5, 0.75, 0.95, 1), na.rm = TRUE)
print(round(quantiles, 4))

cat("\nSites with |logFC difference| > 0.5: ", sum(abs(all_comparisons_df$logFC_diff) > 0.5, na.rm = TRUE), "\n")
cat("Sites with |logFC difference| > 1.0: ", sum(abs(all_comparisons_df$logFC_diff) > 1.0, na.rm = TRUE), "\n")

# =============================================================================
# Significance Changes
# =============================================================================

cat("\n=== SIGNIFICANCE CHANGES (adj.P.Val < 0.05) ===\n\n")

sig_comparison <- all_comparisons_df %>%
  mutate(
    sig_opt1 = adj.P.Val_opt1 < 0.05,
    sig_opt2 = adj.P.Val_opt2 < 0.05,
    sig_change = case_when(
      sig_opt1 & sig_opt2 ~ "Both significant",
      sig_opt1 & !sig_opt2 ~ "Lost significance",
      !sig_opt1 & sig_opt2 ~ "Gained significance",
      TRUE ~ "Both non-significant"
    )
  )

sig_summary <- sig_comparison %>%
  group_by(sig_change) %>%
  summarize(n = n(), .groups = "drop") %>%
  mutate(pct = round(n / sum(n) * 100, 1))

print(sig_summary)

# =============================================================================
# Save Results
# =============================================================================

output_dir <- paste0(source_file_path, "filtering_comparison/")
dir.create(output_dir, showWarnings = FALSE)

# Save comparison statistics
write_csv(comparison_stats, paste0(output_dir, "comparison_statistics.csv"))

# Save detailed comparisons
write_csv(all_comparisons_df, paste0(output_dir, "site_level_comparison.csv"))

# Save Option 2 DE results (PSM-level filtered)
for (cell_type in names(option2_de)) {
  write_csv(option2_de[[cell_type]],
            paste0(output_dir, "OGlcNAc_site_DE_", cell_type, "_PSM_filtered.csv"))
}

# Save Option 2 quantification results
for (cell_type in names(option2_results)) {
  write_csv(option2_results[[cell_type]],
            paste0(output_dir, "OGlcNAc_site_quant_", cell_type, "_PSM_filtered.csv"))
}

cat("\n=== Results saved to: ", output_dir, " ===\n")

# =============================================================================
# Generate Comparison Plot
# =============================================================================

library(ggplot2)

# Scatter plot of logFC values
p <- ggplot(all_comparisons_df, aes(x = logFC_opt1, y = logFC_opt2, color = cell_type)) +
  geom_point(alpha = 0.5, size = 1.5) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50") +
  labs(
    x = "logFC (Option 1: Site-level filtering)",
    y = "logFC (Option 2: PSM-level filtering)",
    title = paste0("Comparison of Filtering Approaches (r = ", round(overall_corr, 3), ")"),
    color = "Cell Type"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  ) +
  coord_fixed()

ggsave(paste0(output_dir, "logFC_comparison_scatter.pdf"), p, width = 5, height = 5)

# Histogram of differences
p2 <- ggplot(all_comparisons_df, aes(x = logFC_diff, fill = cell_type)) +
  geom_histogram(binwidth = 0.05, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey30") +
  labs(
    x = "logFC difference (Option 1 - Option 2)",
    y = "Count",
    title = "Distribution of logFC Differences",
    fill = "Cell Type"
  ) +
  theme_bw() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5)
  )

ggsave(paste0(output_dir, "logFC_difference_histogram.pdf"), p2, width = 5, height = 4)

cat("Comparison plots saved.\n")

# =============================================================================
# Final Recommendation
# =============================================================================

cat("\n=== RECOMMENDATION ===\n\n")

if (overall_corr > 0.99 && overall_mean_diff < 0.05) {
  cat("The two approaches produce nearly identical results.\n")
  cat("Either approach is acceptable for publication.\n")
  cat("Option 1 (site-level filtering) is simpler and retains more data.\n")
} else if (overall_corr > 0.95 && overall_mean_diff < 0.1) {
  cat("The two approaches produce highly correlated results with minor differences.\n")
  cat("Option 2 (PSM-level filtering) is more conservative and defensible.\n")
  cat("Consider using Option 2 for maximum rigor.\n")
} else {
  cat("The two approaches show meaningful differences.\n")
  cat("Option 2 (PSM-level filtering) is recommended for publication.\n")
  cat("Review sites with large differences carefully.\n")
}
