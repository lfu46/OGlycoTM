# OGalNAc_GSEA_exclusive.R
# GSEA for exclusive O-GalNAc proteins (not in O-GlcNAc dataset)

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'

for (cell in c("HEK293T", "HepG2", "Jurkat")) {
  cat("\n", strrep("=", 60), "\n")
  cat("GSEA (exclusive O-GalNAc):", cell, "\n")
  cat(strrep("=", 60), "\n")

  de_gal <- read_csv(paste0(source_file_path, "differential_analysis/OGalNAc_protein_DE_", cell, ".csv"), show_col_types = FALSE)
  de_glc <- read_csv(paste0(source_file_path, "differential_analysis/OGlcNAc_protein_DE_", cell, ".csv"), show_col_types = FALSE)

  # Exclusive O-GalNAc
  exclusive <- de_gal %>% filter(!Protein.ID %in% de_glc$Protein.ID)
  cat("  Exclusive O-GalNAc:", nrow(exclusive), "proteins\n")

  # Map to Entrez
  id_map <- bitr(exclusive$Protein.ID, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  mapped <- exclusive %>% inner_join(id_map, by = c("Protein.ID" = "UNIPROT"))

  gene_list <- mapped %>%
    group_by(ENTREZID) %>%
    slice_max(abs(logFC), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(logFC))

  ranked <- gene_list$logFC
  names(ranked) <- gene_list$ENTREZID
  cat("  Ranked list:", length(ranked), "genes\n")

  if (length(ranked) < 15) {
    cat("  Too few genes for GSEA, skipping\n")
    next
  }

  for (ont in c("BP", "CC", "MF")) {
    gsea <- gseGO(
      geneList = ranked, OrgDb = org.Hs.eg.db, ont = ont,
      minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.1, verbose = FALSE
    )
    sig <- gsea@result %>% filter(p.adjust < 0.05)
    cat("\n  ", ont, ":", nrow(sig), "significant terms\n")
    if (nrow(sig) > 0) {
      sig %>% arrange(p.adjust) %>% head(10) %>%
        mutate(dir = ifelse(NES > 0, "UP", "DOWN")) %>%
        dplyr::select(dir, NES, Description, p.adjust, setSize) %>%
        print()
    }
    write_csv(gsea@result, paste0(source_file_path, "enrichment/OGalNAc_exclusive_GSEA_", ont, "_", cell, ".csv"))
  }
}

cat("\nDone!\n")
