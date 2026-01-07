# import packages
library(tidyverse)

# source file path
source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'

# raw data
OGlyco_HEK293T_raw <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/OGlyco_HEK293T_raw.csv'
)
OGlyco_HepG2_raw <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/OGlyco_HepG2_raw.csv'
)
OGlyco_Jurkat_raw <- read_csv(
  '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/OGlyco_Jurkat_raw.csv'
)

