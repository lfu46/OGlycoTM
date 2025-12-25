# import packages
library(tidyverse)
library(jsonlite)

# import options.json
options <- fromJSON('/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/options.json')

# extract diff_mods and static_mods
diff_mods_df <- as.data.frame(options$diff_mods)
static_mods_df <- as.data.frame(options$static_mods)

# create mapping dataframes
diff_mod_mapping <- tibble(
  mass = diff_mods_df$mass,
  symbol = diff_mods_df$symbol
)

static_mod_mapping <- tibble(
  mass = static_mods_df$mass
)

# Function to parse Assigned.Modifications and extract position-mass mappings
# Format: "7C(555.2968), 12S(528.2859)" or "N-term(229.1629)" or "8K(229.1629), N-term(229.1629)"
parse_assigned_mods <- function(assigned_mods_string) {
  if (is.na(assigned_mods_string) || assigned_mods_string == "") {
    return(tibble(position = character(), mass = numeric()))
  }
  
  # Split by comma and extract position and mass
  mod_parts <- str_split(assigned_mods_string, ",\\s*")[[1]]
  
  positions <- character()
  masses <- numeric()
  
  for (mod_part in mod_parts) {
    mod_part <- str_trim(mod_part)  # Remove leading/trailing whitespace
    
    # Match pattern: position (number) + AA + (mass) OR N-term/C-term + (mass)
    
    if (str_detect(mod_part, "^N-term|^C-term")) {
      # Terminal modification: extract mass from within parentheses
      mass_value <- suppressWarnings(as.numeric(str_extract(mod_part, "(?<=\\()[0-9.]+(?=\\))")))
      if (!is.na(mass_value)) {
        if (str_detect(mod_part, "^N-term")) {
          positions <- c(positions, "N-term")
        } else {
          positions <- c(positions, "C-term")
        }
        masses <- c(masses, mass_value)
      }
    } else {
      # Position-specific modification: e.g., "6S(555.2968)"
      # Extract position as first one or more digits
      pos_value <- suppressWarnings(as.integer(str_extract(mod_part, "^[0-9]+")))
      # Extract mass from within parentheses
      mass_value <- suppressWarnings(as.numeric(str_extract(mod_part, "(?<=\\()[0-9.]+(?=\\))")))
      
      if (!is.na(pos_value) && !is.na(mass_value)) {
        positions <- c(positions, as.character(pos_value))
        masses <- c(masses, mass_value)
      }
    }
  }
  
  return(tibble(position = positions, mass = masses))
}

# Function to replace modifications in peptide string
replace_modifications <- function(pep_string, assigned_mods_string) {
  # Parse the assigned modifications to get mass values
  mods_present <- parse_assigned_mods(assigned_mods_string)
  
  result <- pep_string
  
  # Track which modifications we've already used
  mods_used <- rep(FALSE, nrow(mods_present))
  
  # Handle N-terminal modifications first
  if (str_detect(result, "^[nc]\\[")) {
    # Remove the entire n[XXX] or c[XXX] (whether it has a modification or not)
    result <- str_remove(result, "^[nc]\\[[0-9]+\\]")
    # Mark N-term as used if it exists
    n_term_idx <- which(mods_present$position == "N-term")
    if (length(n_term_idx) > 0) {
      mods_used[n_term_idx] <- TRUE
    }
  }
  
  # Now process all remaining AA[XXX] brackets
  # Process them in order of appearance (left to right)
  # and match them to modifications in order of appearance in Assigned.Modifications
  
  repeat {
    # Find first [XXX] pattern
    match <- str_locate(result, "\\[[0-9]+\\]")
    
    if (is.na(match[1, 1])) {
      break  # No more matches
    }
    
    start_pos <- match[1, "start"]
    end_pos <- match[1, "end"]
    
    # Get the amino acid before the bracket
    aa <- substr(result, start_pos - 1, start_pos - 1)
    
    replacement <- NA
    matched_idx <- NA
    
    # Find the first unused modification that hasn't been matched yet
    # This ensures brackets are matched in the order they appear
    for (k in 1:nrow(mods_present)) {
      if (!mods_used[k] && mods_present$position[k] != "N-term" && mods_present$position[k] != "C-term") {
        mod_mass <- mods_present$mass[k]
        
        # Check if this mass is a diff_mod
        for (j in 1:nrow(diff_mod_mapping)) {
          diff_mass <- diff_mod_mapping$mass[j]
          if (abs(diff_mass - mod_mass) < 0.01) {
            # Found matching diff_mod
            replacement <- paste0(aa, diff_mod_mapping$symbol[j])
            matched_idx <- k
            break
          }
        }
        
        # If not found, check static_mods
        if (is.na(replacement)) {
          for (j in 1:nrow(static_mod_mapping)) {
            static_mass <- static_mod_mapping$mass[j]
            if (abs(static_mass - mod_mass) < 0.01) {
              # Found matching static_mod
              replacement <- aa
              matched_idx <- k
              break
            }
          }
        }
        
        # If we found a match, use it and stop looking
        if (!is.na(replacement)) {
          break
        }
      }
    }
    
    # If still no match, remove bracket as safety fallback
    if (is.na(replacement)) {
      replacement <- aa
    }
    
    # Mark this modification as used
    if (!is.na(matched_idx)) {
      mods_used[matched_idx] <- TRUE
    }
    
    # Replace in string
    before <- substr(result, 1, start_pos - 2)
    after <- substr(result, end_pos + 1, nchar(result))
    result <- paste0(before, replacement, after)
  }
  
  return(result)
}

# Test the function on a few examples first
cat("\n=== Testing replacement function ===\n")

test_result_1 <- replace_modifications(
  "n[230]SENPTSQASQ",
  "N-term(229.1629)"
)
cat("Test 1: n[230]SENPTSQASQ with N-term(229.1629)\n")
cat("Expected: SENPTSQASQ\n")
cat("Result: ", test_result_1, "\n")
cat("Pass:", test_result_1 == "SENPTSQASQ", "\n\n")

test_result_2 <- replace_modifications(
  "SENPTS[642]QASQ",
  "6S(555.2968)"
)
cat("Test 2: SENPTS[642]QASQ with 6S(555.2968)\n")
cat("Expected: SENPTS~QASQ\n")
cat("Result: ", test_result_2, "\n")
cat("Pass:", test_result_2 == "SENPTS~QASQ", "\n\n")

test_result_3 <- replace_modifications(
  "n[230]SENPTM[147]SQASQ",
  "N-term(229.1629), 6M(15.9949)"
)
cat("Test 3: n[230]SENPTM[147]SQASQ with N-term(229.1629), 6M(15.9949)\n")
cat("Expected: SENPTM*SQASQ\n")
cat("Result: ", test_result_3, "\n")
cat("Pass:", test_result_3 == "SENPTM*SQASQ", "\n\n")

# Test 4: n[230]SENPTSQAS[615]Q
# After removing n[230]: SENPTSQAS[615]Q
# S E N P T S Q A S[615] Q
# 1 2 3 4 5 6 7 8 9      10
test_result_4 <- replace_modifications(
  "n[230]SENPTSQAS[615]Q",
  "N-term(229.1629), 9S(528.2859)"
)
cat("Test 4: n[230]SENPTSQAS[615]Q with N-term(229.1629), 9S(528.2859)\n")
cat("Expected: SENPTSQAS^Q\n")
cat("Result: ", test_result_4, "\n")
cat("Pass:", test_result_4 == "SENPTSQAS^Q", "\n\n")

# Test 5: n[230]S[615]ENPT[656]QAC[429]Q
# After removing n[230]: S[615]ENPT[656]QAC[429]Q
# S[615] E N P T[656] Q A C[429] Q
# 1      2 3 4 5      6 7 8      9
test_result_5 <- replace_modifications(
  "n[230]S[615]ENPT[656]QAC[429]Q",
  "N-term(229.1629), 1S(528.2859), 5T(555.2968), 8C(326.1339)"
)
cat("Test 5: n[230]S[615]ENPT[656]QAC[429]Q with N-term, S(528.2859), T(555.2968), C(326.1339)\n")
cat("Expected: S^ENPT~QAC#Q\n")
cat("Result: ", test_result_5, "\n")
cat("Pass:", test_result_5 == "S^ENPT~QAC#Q", "\n\n")

# Test 6: n[230]PEPTC[160]ASQC[658]
# After removing n[230]: PEPTC[160]ASQC[658]
# P E P T C[160] A S Q C[658]
# 1 2 3 4 5      6 7 8 9 10
test_result_6 <- replace_modifications(
  "n[230]PEPTC[160]ASQC[658]",
  "N-term(229.1629), 5C(57.0214), 10C(555.2968)"
)
cat("Test 6: n[230]PEPTC[160]ASQC[658] with N-term, C(57.0214), C(555.2968)\n")
cat("Expected: PEPTC$ASQC~\n")
cat("Result: ", test_result_6, "\n")
cat("Pass:", test_result_6 == "PEPTC$ASQC~", "\n\n")

# Test 7: n[230]M[147]ENPTS[386]SQASQ
# After removing n[230]: M[147]ENPTS[386]SQASQ
# M[147] E N P T S[386] S Q A S Q
# 1      2 3 4 5 6      7 8 9 10 11
test_result_7 <- replace_modifications(
  "n[230]M[147]ENPTS[386]SQASQ",
  "N-term(229.1629), 1M(15.9949), 6S(299.1230)"
)
cat("Test 7: n[230]M[147]ENPTS[386]SQASQ with N-term, M(15.9949), S(299.1230)\n")
cat("Expected: M*ENPTS@SQASQ\n")
cat("Result: ", test_result_7, "\n")
cat("Pass:", test_result_7 == "M*ENPTS@SQASQ", "\n\n")

# Test 8: n[230]S[386]ENPT[427]S[642]QASQ
# After removing n[230]: S[386]ENPT[427]S[642]QASQ
# S[386] E N P T[427] S[642] Q A S Q
# 1      2 3 4 5      6      7 8 9 10
test_result_8 <- replace_modifications(
  "n[230]S[386]ENPT[427]S[642]QASQ",
  "N-term(229.1629), 1S(299.1230), 5T(326.1339), 6S(555.2968)"
)
cat("Test 8: n[230]S[386]ENPT[427]S[642]QASQ with N-term, S(299.1230), T(326.1339), S(555.2968)\n")
cat("Expected: S@ENPT#S~QASQ\n")
cat("Result: ", test_result_8, "\n")
cat("Pass:", test_result_8 == "S@ENPT#S~QASQ", "\n\n")

# import Eclipse_LF_OGlycoTM_HEK293T_OG_1_10062025 results
Eclipse_LF_OGlycoTM_HEK293T_OG_1_10062025 <- read_tsv(
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/OGlycoTM_HEK293T/OGlyco/OGlyco_OPair_Search/OGlyco/psm.tsv',
  col_names = TRUE,
  name_repair = 'universal'
) |> 
  filter(str_detect(Spectrum, 'Eclipse_LF_OGlycoTM_HEK293T_OG_1_10062025')) |> 
  select(Spectrum, Modified.Peptide, Extended.Peptide, Observed.M.Z, Assigned.Modifications) |> 
  mutate(
    Scan.Number = as.integer(str_extract(Spectrum, '(?<=\\.)[0-9]+(?=\\.)')),
    .after = Spectrum
  )

# Apply the replacement function to full dataset
Eclipse_LF_OGlycoTM_HEK293T_OG_1_10062025 <- Eclipse_LF_OGlycoTM_HEK293T_OG_1_10062025 |>
  mutate(
    Modified.Peptide.Annotated = map2_chr(
      Modified.Peptide,
      Assigned.Modifications,
      replace_modifications
    ),
    # Then, create Extended Peptide format with one position before and after the dot
    # Format: "PP.SENPTSQASQ.-" -> extract P (before first dot) and - (after last dot)
    Prev.AA = str_extract(Extended.Peptide, "^."),  # Character before first dot
    Next.AA = str_extract(Extended.Peptide, ".$"),  # Character after last dot
    # Create the annotated Extended Peptide format
    Extended.Peptide.Annotated = paste0(
      Prev.AA, ".",
      Modified.Peptide.Annotated, ".",
      Next.AA
    ),
    .after = Modified.Peptide
  )

# Display results for verification
cat("\n=== Results from full dataset ===\n")
Eclipse_LF_OGlycoTM_HEK293T_OG_1_10062025 |>
  select(Spectrum, Scan.Number, Modified.Peptide, Modified.Peptide.Annotated, Assigned.Modifications) |>
  head(30)

write_tsv(
  Eclipse_LF_OGlycoTM_HEK293T_OG_1_10062025,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/data_source/Eclipse_LF_OGlycoTM_HEK293T_OG_1_10062025.tsv'
)

library(tidyverse)

Eclipse_LF_OGlycoTM_HEK293T_OG_1_10062025 <- read_tsv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/data_source/Eclipse_LF_OGlycoTM_HEK293T_OG_1_10062025.tsv'
)

Eclipse_LF_OGlycoTM_HEK293T_OG_1_10062025_528 <- Eclipse_LF_OGlycoTM_HEK293T_OG_1_10062025 |> 
  filter(str_detect(Modified.Peptide.Annotated, '\\^'))

write_tsv(
  Eclipse_LF_OGlycoTM_HEK293T_OG_1_10062025_528,
  file = '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_HEK293T/data_source/Eclipse_LF_OGlycoTM_HEK293T_OG_1_10062025_528.tsv'
)


