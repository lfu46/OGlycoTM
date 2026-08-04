# OGalNAc_GSEA.R
# GSEA analysis for O-GalNAc protein-level quantification
# Uses ranked logFC values — no cutoff needed, more powerful than ORA

library(tidyverse)
library(clusterProfiler)
library(org.Hs.eg.db)

source_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/data_source/'
figure_file_path <- '/Volumes/cos-lab-rwu60/Longping/OGlycoTM_Final_Version/Figures/'

cell_types <- c("HEK293T", "HepG2", "Jurkat")

for (cell in cell_types) {
  cat("\n", strrep("=", 60), "\n")
  cat("GSEA:", cell, "\n")
  cat(strrep("=", 60), "\n")

  de <- read_csv(
    paste0(source_file_path, "differential_analysis/OGalNAc_protein_DE_", cell, ".csv"),
    show_col_types = FALSE
  )

  # Convert UniProt to Entrez ID
  id_map <- bitr(de$Protein.ID, fromType = "UNIPROT", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
  de_mapped <- de %>%
    inner_join(id_map, by = c("Protein.ID" = "UNIPROT"))

  cat("  Mapped:", nrow(de_mapped), "/", nrow(de), "proteins to Entrez ID\n")

  # Create ranked gene list (logFC, named by Entrez ID)
  # If duplicates, keep the one with highest absolute logFC
  gene_list <- de_mapped %>%
    group_by(ENTREZID) %>%
    slice_max(abs(logFC), n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    arrange(desc(logFC))

  ranked <- gene_list$logFC
  names(ranked) <- gene_list$ENTREZID

  cat("  Ranked list:", length(ranked), "genes\n")
  cat("  logFC range: [", min(ranked), ",", max(ranked), "]\n")

  # GSEA - GO Biological Process
  gsea_bp <- gseGO(
    geneList = ranked,
    OrgDb = org.Hs.eg.db,
    ont = "BP",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.1,
    verbose = FALSE
  )

  # GSEA - GO Cellular Component
  gsea_cc <- gseGO(
    geneList = ranked,
    OrgDb = org.Hs.eg.db,
    ont = "CC",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.1,
    verbose = FALSE
  )

  # GSEA - GO Molecular Function
  gsea_mf <- gseGO(
    geneList = ranked,
    OrgDb = org.Hs.eg.db,
    ont = "MF",
    minGSSize = 10,
    maxGSSize = 500,
    pvalueCutoff = 0.1,
    verbose = FALSE
  )

  # Print results
  for (ont_name in c("BP", "CC", "MF")) {
    gsea_result <- get(paste0("gsea_", tolower(ont_name)))
    sig <- gsea_result@result %>% filter(p.adjust < 0.05)
    cat("\n  ", ont_name, ": ", nrow(sig), " significant terms (adj.P < 0.05)\n")
    if (nrow(sig) > 0) {
      sig %>%
        arrange(p.adjust) %>%
        head(10) %>%
        dplyr::select(Description, NES, p.adjust, setSize) %>%
        print()
    }
  }

  # Save results
  write_csv(gsea_bp@result, paste0(source_file_path, "enrichment/OGalNAc_GSEA_BP_", cell, ".csv"))
  write_csv(gsea_cc@result, paste0(source_file_path, "enrichment/OGalNAc_GSEA_CC_", cell, ".csv"))
  write_csv(gsea_mf@result, paste0(source_file_path, "enrichment/OGalNAc_GSEA_MF_", cell, ".csv"))

  cat("\n  Results saved to enrichment/OGalNAc_GSEA_{BP,CC,MF}_", cell, ".csv\n")
}

cat("\n", strrep("=", 60), "\n")
cat("O-GalNAc GSEA complete!\n")
cat(strrep("=", 60), "\n")
