# import packages
library(tidyverse)
library(edgeR)

# source data
source("data_source.R")


# Define intensity columns ------------------------------------------------
intensity_cols <- c("Intensity.Tuni_1", "Intensity.Tuni_2", "Intensity.Tuni_3",
                    "Intensity.Ctrl_4", "Intensity.Ctrl_5", "Intensity.Ctrl_6")


# Filtering function ------------------------------------------------------
# Removes entries with 0 or NA in any intensity column

filter_zero_na <- function(data, intensity_cols) {
  data |>
    filter(if_all(all_of(intensity_cols), ~ . > 0 & !is.na(.)))
}


# Normalization function --------------------------------------------------
# Performs sample loading (SL) normalization followed by TMM normalization
# Keeps all original columns plus adds _sl and _sl_tmm columns

normalize_intensity <- function(data,
                                intensity_cols = c("Intensity.Tuni_1", "Intensity.Tuni_2", "Intensity.Tuni_3",
                                                   "Intensity.Ctrl_4", "Intensity.Ctrl_5", "Intensity.Ctrl_6")) {

  # Extract intensity matrix
  intensity_matrix <- data |>
    select(all_of(intensity_cols)) |>
    as.matrix()

  # Step 1: Sample Loading (SL) normalization
  # Assumes the sum of each channel should be the same across samples
  col_sums <- colSums(intensity_matrix)
  target_mean <- mean(col_sums)
  sl_factors <- target_mean / col_sums

  # Apply SL normalization
  intensity_sl <- sweep(intensity_matrix, 2, sl_factors, FUN = "*")

  # Create SL-normalized tibble with renamed columns
  intensity_sl_tb <- as_tibble(intensity_sl)
  colnames(intensity_sl_tb) <- paste0(intensity_cols, "_sl")

  # Step 2: TMM normalization
  # Makes the median value similar across samples
  tmm_factors <- calcNormFactors(intensity_sl)

  # Apply TMM normalization
  intensity_sl_tmm <- sweep(intensity_sl, 2, tmm_factors, FUN = "/")

  # Create SL+TMM-normalized tibble with renamed columns
  intensity_sl_tmm_tb <- as_tibble(intensity_sl_tmm)
  colnames(intensity_sl_tmm_tb) <- paste0(intensity_cols, "_sl_tmm")

  # Combine original data with normalized columns
  result <- bind_cols(data, intensity_sl_tb, intensity_sl_tmm_tb)

  return(result)
}


# O-GlcNAc protein level filtering and normalization ----------------------

# HEK293T
OGlcNAc_protein_filtered_HEK293T <- filter_zero_na(OGlcNAc_protein_quant_HEK293T, intensity_cols)
cat("HEK293T protein: ", nrow(OGlcNAc_protein_quant_HEK293T), " -> ", nrow(OGlcNAc_protein_filtered_HEK293T), " after filtering\n")

write_csv(
  OGlcNAc_protein_filtered_HEK293T,
  paste0(source_file_path, "quantification/OGlcNAc_protein_filtered_HEK293T.csv")
)

OGlcNAc_protein_norm_HEK293T <- normalize_intensity(OGlcNAc_protein_filtered_HEK293T)

write_csv(
  OGlcNAc_protein_norm_HEK293T,
  paste0(source_file_path, "normalization/OGlcNAc_protein_norm_HEK293T.csv")
)

# HepG2
OGlcNAc_protein_filtered_HepG2 <- filter_zero_na(OGlcNAc_protein_quant_HepG2, intensity_cols)
cat("HepG2 protein: ", nrow(OGlcNAc_protein_quant_HepG2), " -> ", nrow(OGlcNAc_protein_filtered_HepG2), " after filtering\n")

write_csv(
  OGlcNAc_protein_filtered_HepG2,
  paste0(source_file_path, "quantification/OGlcNAc_protein_filtered_HepG2.csv")
)

OGlcNAc_protein_norm_HepG2 <- normalize_intensity(OGlcNAc_protein_filtered_HepG2)

write_csv(
  OGlcNAc_protein_norm_HepG2,
  paste0(source_file_path, "normalization/OGlcNAc_protein_norm_HepG2.csv")
)

# Jurkat
OGlcNAc_protein_filtered_Jurkat <- filter_zero_na(OGlcNAc_protein_quant_Jurkat, intensity_cols)
cat("Jurkat protein: ", nrow(OGlcNAc_protein_quant_Jurkat), " -> ", nrow(OGlcNAc_protein_filtered_Jurkat), " after filtering\n")

write_csv(
  OGlcNAc_protein_filtered_Jurkat,
  paste0(source_file_path, "quantification/OGlcNAc_protein_filtered_Jurkat.csv")
)

OGlcNAc_protein_norm_Jurkat <- normalize_intensity(OGlcNAc_protein_filtered_Jurkat)

write_csv(
  OGlcNAc_protein_norm_Jurkat,
  paste0(source_file_path, "normalization/OGlcNAc_protein_norm_Jurkat.csv")
)


# O-GlcNAc site level filtering and normalization -------------------------

# HEK293T
OGlcNAc_site_filtered_HEK293T <- filter_zero_na(OGlcNAc_site_quant_HEK293T, intensity_cols)
cat("HEK293T site: ", nrow(OGlcNAc_site_quant_HEK293T), " -> ", nrow(OGlcNAc_site_filtered_HEK293T), " after filtering\n")

write_csv(
  OGlcNAc_site_filtered_HEK293T,
  paste0(source_file_path, "quantification/OGlcNAc_site_filtered_HEK293T.csv")
)

OGlcNAc_site_norm_HEK293T <- normalize_intensity(OGlcNAc_site_filtered_HEK293T)

write_csv(
  OGlcNAc_site_norm_HEK293T,
  paste0(source_file_path, "normalization/OGlcNAc_site_norm_HEK293T.csv")
)

# HepG2
OGlcNAc_site_filtered_HepG2 <- filter_zero_na(OGlcNAc_site_quant_HepG2, intensity_cols)
cat("HepG2 site: ", nrow(OGlcNAc_site_quant_HepG2), " -> ", nrow(OGlcNAc_site_filtered_HepG2), " after filtering\n")

write_csv(
  OGlcNAc_site_filtered_HepG2,
  paste0(source_file_path, "quantification/OGlcNAc_site_filtered_HepG2.csv")
)

OGlcNAc_site_norm_HepG2 <- normalize_intensity(OGlcNAc_site_filtered_HepG2)

write_csv(
  OGlcNAc_site_norm_HepG2,
  paste0(source_file_path, "normalization/OGlcNAc_site_norm_HepG2.csv")
)

# Jurkat
OGlcNAc_site_filtered_Jurkat <- filter_zero_na(OGlcNAc_site_quant_Jurkat, intensity_cols)
cat("Jurkat site: ", nrow(OGlcNAc_site_quant_Jurkat), " -> ", nrow(OGlcNAc_site_filtered_Jurkat), " after filtering\n")

write_csv(
  OGlcNAc_site_filtered_Jurkat,
  paste0(source_file_path, "quantification/OGlcNAc_site_filtered_Jurkat.csv")
)

OGlcNAc_site_norm_Jurkat <- normalize_intensity(OGlcNAc_site_filtered_Jurkat)

write_csv(
  OGlcNAc_site_norm_Jurkat,
  paste0(source_file_path, "normalization/OGlcNAc_site_norm_Jurkat.csv")
)


# Verification: Check column sums -----------------------------------------
# Uncomment to verify normalization worked correctly

# cat("\n=== HEK293T Protein Level Normalization Check ===\n")
# cat("Original column sums:\n")
# print(OGlcNAc_protein_norm_HEK293T |> select(Intensity.Tuni_1:Intensity.Ctrl_6) |> colSums())
# cat("\nSL-normalized column sums:\n")
# print(OGlcNAc_protein_norm_HEK293T |> select(ends_with("_sl")) |> colSums())
# cat("\nSL+TMM-normalized column sums:\n")
# print(OGlcNAc_protein_norm_HEK293T |> select(ends_with("_sl_tmm")) |> colSums())

