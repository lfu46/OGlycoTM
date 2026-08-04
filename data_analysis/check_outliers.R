library(tidyverse)

source_file_path <- "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/"

parse_site_probability <- function(prob_str) {
  prob <- str_extract(prob_str, "[0-9.]+(?=\\]$)")
  as.numeric(prob)
}

# Check the outlier sites in Jurkat
site_data <- read_csv(paste0(source_file_path, "site/OGlcNAc_site_Jurkat.csv"), show_col_types = FALSE)
site_data <- site_data %>%
  mutate(Site_Prob = parse_site_probability(Site.Probabilities))

# Look at O75376_S2069
cat("=== O75376_S2069 (Jurkat) - Largest difference (9.58) ===\n")
outlier <- site_data %>%
  filter(site_index == "O75376_S2069") %>%
  select(site_index, Confidence.Level, Site_Prob,
         Intensity.Tuni_1, Intensity.Tuni_2, Intensity.Tuni_3,
         Intensity.Ctrl_4, Intensity.Ctrl_5, Intensity.Ctrl_6)

cat("All PSMs:\n")
print(outlier)

all_tuni <- sum(outlier$Intensity.Tuni_1) + sum(outlier$Intensity.Tuni_2) + sum(outlier$Intensity.Tuni_3)
all_ctrl <- sum(outlier$Intensity.Ctrl_4) + sum(outlier$Intensity.Ctrl_5) + sum(outlier$Intensity.Ctrl_6)
cat("\nTotal intensity (all PSMs): Tuni =", all_tuni, ", Ctrl =", all_ctrl, "\n")
cat("Ratio (all PSMs): ", all_tuni/all_ctrl, "\n")

# Filter high-conf only
high_conf <- outlier %>%
  filter(Confidence.Level == "Level1" | (Confidence.Level == "Level1b" & Site_Prob >= 0.75))

cat("\nHigh-conf PSMs only:\n")
print(high_conf)

if (nrow(high_conf) > 0) {
  hc_tuni <- sum(high_conf$Intensity.Tuni_1) + sum(high_conf$Intensity.Tuni_2) + sum(high_conf$Intensity.Tuni_3)
  hc_ctrl <- sum(high_conf$Intensity.Ctrl_4) + sum(high_conf$Intensity.Ctrl_5) + sum(high_conf$Intensity.Ctrl_6)
  cat("Total intensity (high-conf): Tuni =", hc_tuni, ", Ctrl =", hc_ctrl, "\n")
  cat("Ratio (high-conf): ", hc_tuni/hc_ctrl, "\n")
} else {
  cat("No high-confidence PSMs! All PSMs are low-prob Level1b.\n")
}

# Also check P04080_T16
cat("\n=== P04080_T16 (Jurkat) - Second largest difference (7.43) ===\n")
outlier2 <- site_data %>%
  filter(site_index == "P04080_T16") %>%
  select(site_index, Confidence.Level, Site_Prob,
         Intensity.Tuni_1, Intensity.Tuni_2, Intensity.Tuni_3,
         Intensity.Ctrl_4, Intensity.Ctrl_5, Intensity.Ctrl_6)
print(outlier2)

high_conf2 <- outlier2 %>%
  filter(Confidence.Level == "Level1" | (Confidence.Level == "Level1b" & Site_Prob >= 0.75))
cat("\nHigh-conf PSMs:\n")
print(high_conf2)

if (nrow(high_conf2) > 0) {
  hc2_tuni <- sum(high_conf2$Intensity.Tuni_1) + sum(high_conf2$Intensity.Tuni_2) + sum(high_conf2$Intensity.Tuni_3)
  hc2_ctrl <- sum(high_conf2$Intensity.Ctrl_4) + sum(high_conf2$Intensity.Ctrl_5) + sum(high_conf2$Intensity.Ctrl_6)
  cat("Ratio (high-conf): ", hc2_tuni/hc2_ctrl, "\n")
}

# Summary of what's happening
cat("\n=== EXPLANATION ===\n")
cat("Large differences occur when:\n")
cat("1. Low-prob PSMs have different intensity patterns than high-prob PSMs\n")
cat("2. Most intensity comes from low-prob PSMs that get filtered out\n")
cat("3. The remaining high-prob PSMs have very different Tuni/Ctrl ratios\n")
