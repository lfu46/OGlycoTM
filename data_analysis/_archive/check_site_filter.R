library(readr)
library(stringr)

source_file_path <- "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/"

for (ct in c("HEK293T", "HepG2", "Jurkat")) {
  site_data <- read_csv(paste0(source_file_path, "site/OGlcNAc_site_", ct, ".csv"),
                        show_col_types = FALSE)
  site_data$Site_Prob <- as.numeric(str_extract(site_data$Site.Probabilities, "[0-9.]+(?=\\]$)"))

  n_total <- nrow(site_data)
  n_l1 <- sum(site_data$Confidence.Level == "Level1", na.rm = TRUE)
  n_l1b <- sum(site_data$Confidence.Level == "Level1b", na.rm = TRUE)
  n_l1b_high <- sum(site_data$Confidence.Level == "Level1b" & site_data$Site_Prob >= 0.75, na.rm = TRUE)
  n_l1b_low <- sum(site_data$Confidence.Level == "Level1b" & site_data$Site_Prob < 0.75, na.rm = TRUE)
  n_filtered <- n_l1 + n_l1b_high

  old_count <- c("HEK293T" = 1046, "HepG2" = 600, "Jurkat" = 731)[ct]

  cat(sprintf("%s: total PSMs = %d (old table = %d, match = %s)\n",
    ct, n_total, old_count, n_total == old_count))
  cat(sprintf("  Level1 = %d, Level1b = %d (>=0.75: %d, <0.75: %d)\n",
    n_l1, n_l1b, n_l1b_high, n_l1b_low))
  cat(sprintf("  After filter (L1 + L1b>=0.75) = %d\n\n", n_filtered))
}
