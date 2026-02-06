# analyze_low_prob_psm.R
# Extract and summarize PSMs with confidence Level1 or Level1b and site probability < 0.75

library(tidyverse)

# =============================================================================
# File paths
# =============================================================================

site_files <- c(
  "HEK293T" = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/site/OGlcNAc_site_HEK293T.csv",
  "HepG2" = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/site/OGlcNAc_site_HepG2.csv",
  "Jurkat" = "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/site/OGlcNAc_site_Jurkat.csv"
)

# =============================================================================
# Function to parse site probability from string like "[5,I1V1,0.823]"
# =============================================================================

parse_site_probability <- function(prob_str) {
  prob <- str_extract(prob_str, "[0-9.]+(?=\\]$)")
  as.numeric(prob)
}

# =============================================================================
# Read and process data for each cell type
# =============================================================================

all_data <- list()

for (cell_type in names(site_files)) {
  cat("Processing", cell_type, "...\n")

  data <- read_csv(site_files[cell_type], show_col_types = FALSE)

  data <- data %>%
    mutate(
      Cell_Type = cell_type,
      Site_Prob = parse_site_probability(Site.Probabilities)
    )

  all_data[[cell_type]] <- data
  cat("  Total PSMs:", nrow(data), "\n")
}

combined_data <- bind_rows(all_data)

# =============================================================================
# Filter for Level1 or Level1b
# =============================================================================

level1_1b_data <- combined_data %>%
  filter(Confidence.Level %in% c("Level1", "Level1b"))

low_prob_data <- level1_1b_data %>%
  filter(Site_Prob < 0.75)

cat("\n")
cat(strrep("=", 70), "\n")
cat("PSMs with Confidence Level1/Level1b and Site Probability < 0.75\n")
cat(strrep("=", 70), "\n\n")

# =============================================================================
# Table 1: Proportion of Low Probability PSMs (< 0.75)
# =============================================================================

table1 <- level1_1b_data %>%
  group_by(Cell_Type, Confidence.Level) %>%
  summarise(
    Total_PSM = n(),
    Low_Prob_PSM = sum(Site_Prob < 0.75, na.rm = TRUE),
    Pct_Low_Prob = round(Low_Prob_PSM / Total_PSM * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(Cell_Type, Confidence.Level)

cat("Table 1: Proportion of Low Probability PSMs (< 0.75)\n")
cat(strrep("-", 70), "\n")
print(as.data.frame(table1))

# =============================================================================
# Table 2: Sites that would be lost if low prob PSMs are removed
# For each site with low prob PSMs, check if it has any high prob PSMs
# =============================================================================

# Get all sites from Level1/Level1b
sites_summary <- level1_1b_data %>%
  group_by(Cell_Type, site_index) %>%
  summarise(
    has_high_prob = any(Site_Prob >= 0.75, na.rm = TRUE),
    has_low_prob = any(Site_Prob < 0.75, na.rm = TRUE),
    .groups = "drop"
  )

# Sites that ONLY have low prob PSMs (would be lost if low prob removed)
sites_only_low_prob <- sites_summary %>%
  filter(has_low_prob & !has_high_prob)

# Sites that have BOTH low and high prob PSMs (would be retained)
sites_with_both <- sites_summary %>%
  filter(has_low_prob & has_high_prob)

# Calculate summary by cell type
table2 <- sites_summary %>%
  filter(has_low_prob) %>%  # Sites that have at least one low prob PSM
  group_by(Cell_Type) %>%
  summarise(
    Total_Sites_with_Low_Prob_PSM = n(),
    Sites_Only_Low_Prob = sum(!has_high_prob),
    Sites_Also_High_Prob = sum(has_high_prob),
    .groups = "drop"
  )

# Add total row
total_row <- sites_summary %>%
  filter(has_low_prob) %>%
  summarise(
    Cell_Type = "Total",
    Total_Sites_with_Low_Prob_PSM = n(),
    Sites_Only_Low_Prob = sum(!has_high_prob),
    Sites_Also_High_Prob = sum(has_high_prob)
  )

table2 <- bind_rows(table2, total_row)

cat("\n\nTable 2: Impact on Sites if Low Prob PSMs (< 0.75) are Removed\n")
cat(strrep("-", 70), "\n")
cat("Sites_Only_Low_Prob: Sites that would be LOST (no high prob PSM support)\n")
cat("Sites_Also_High_Prob: Sites that would be RETAINED (have high prob PSM)\n\n")
print(as.data.frame(table2))

# =============================================================================
# Table 3: Summary by Confidence Level (across all cell types)
# =============================================================================

table3 <- table1 %>%
  group_by(Confidence.Level) %>%
  summarise(
    Total_PSM = sum(Total_PSM),
    Low_Prob_PSM = sum(Low_Prob_PSM),
    Pct_Low_Prob = round(Low_Prob_PSM / Total_PSM * 100, 1),
    .groups = "drop"
  )

cat("\n\nTable 3: Summary by Confidence Level (All Cell Types)\n")
cat(strrep("-", 70), "\n")
print(as.data.frame(table3))

# =============================================================================
# Save extracted low probability PSMs
# =============================================================================

output_file <- "/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/site/OGlcNAc_low_prob_psm.csv"

low_prob_data %>%
  select(
    Cell_Type,
    Spectrum,
    Protein.ID,
    Gene,
    site_index,
    Peptide,
    Confidence.Level,
    Site.Probabilities,
    Site_Prob,
    Observed.M.Z,
    Charge,
    Hyperscore
  ) %>%
  write_csv(output_file)

cat("\n\nLow probability PSMs saved to:\n", output_file, "\n")
cat("Total rows:", nrow(low_prob_data), "\n")
