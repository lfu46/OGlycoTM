library(openxlsx)
table_dir <- "/Users/longpingfu/Downloads/OGlycoTM/Manuscript/version_Mar24_2026/"

# All site examples mentioned in the manuscript:
# Figure 5: O-GlcNAc sites
#   PRDX6 Y89 (upregulated 2.83-fold)
#   EWSR1 S274 (fold change = 0.90)
# Figure 6: O-GalNAc sites
#   TFRC T104
#   CPD T44
#   PTPRC T139
#   PTPRC S146

cat("=== Checking manuscript site examples in supporting tables ===\n\n")

# O-GlcNAc sites: check S5 (site ID) and S6 (site DE)
oglcnac_sites <- list(
  list(gene="PRDX6", site="Y89", expected_fc=2.83),
  list(gene="EWSR1", site="S274", expected_fc=0.90)
)

cat("--- O-GlcNAc sites (S5 = site ID, S6 = site DE) ---\n\n")
cells <- c("HEK293T", "HepG2", "Jurkat")

for (ex in oglcnac_sites) {
  cat(sprintf("%s %s:\n", ex$gene, ex$site))

  # Check S5 (site identification)
  found_s5 <- FALSE
  for (i in 1:3) {
    d <- read.xlsx(paste0(table_dir, "supporting_table_S5.xlsx"), sheet=i, startRow=2, colNames=TRUE)
    # col 2 = Gene Symbol, col 3 = Site
    matches <- which(d[,2] == ex$gene & d[,3] == ex$site)
    if (length(matches) > 0) {
      found_s5 <- TRUE
      cat(sprintf("  S5%s (%s): FOUND (%d PSMs)\n", LETTERS[i], cells[i], length(matches)))
    }
  }
  if (!found_s5) cat("  S5: *** NOT FOUND IN ANY CELL TYPE ***\n")

  # Check S6 (site DE)
  found_s6 <- FALSE
  for (i in 1:3) {
    d <- read.xlsx(paste0(table_dir, "supporting_table_S6.xlsx"), sheet=i, startRow=2, colNames=TRUE)
    # col 2 = Gene Symbol, col 3 = Site
    matches <- which(d[,2] == ex$gene & d[,3] == ex$site)
    if (length(matches) > 0) {
      found_s6 <- TRUE
      logfc <- d[matches[1], 4]
      fc <- 2^logfc
      cat(sprintf("  S6%s (%s): FOUND, log2FC=%.2f, FC=%.2f (manuscript=%.2f)\n",
        LETTERS[i], cells[i], logfc, fc, ex$expected_fc))
    }
  }
  if (!found_s6) cat("  S6: *** NOT FOUND IN ANY CELL TYPE ***\n")
  cat("\n")
}

# O-GalNAc sites: check S9 (site ID) and S10 (site DE)
ogalnac_sites <- list(
  list(gene="TFRC", site="T104"),
  list(gene="CPD", site="T44"),
  list(gene="PTPRC", site="T139"),
  list(gene="PTPRC", site="S146")
)

cat("--- O-GalNAc sites (S9 = site ID, S10 = site DE) ---\n\n")

for (ex in ogalnac_sites) {
  cat(sprintf("%s %s:\n", ex$gene, ex$site))

  # Check S9 (site identification)
  found_s9 <- FALSE
  for (i in 1:3) {
    d <- read.xlsx(paste0(table_dir, "supporting_table_S9.xlsx"), sheet=i, startRow=2, colNames=TRUE)
    matches <- which(d[,2] == ex$gene & d[,3] == ex$site)
    if (length(matches) > 0) {
      found_s9 <- TRUE
      cat(sprintf("  S9%s (%s): FOUND (%d PSMs)\n", LETTERS[i], cells[i], length(matches)))
    }
  }
  if (!found_s9) cat("  S9: *** NOT FOUND IN ANY CELL TYPE ***\n")

  # Check S10 (site DE)
  found_s10 <- FALSE
  for (i in 1:3) {
    d <- read.xlsx(paste0(table_dir, "supporting_table_S10.xlsx"), sheet=i, startRow=2, colNames=TRUE)
    matches <- which(d[,2] == ex$gene & d[,3] == ex$site)
    if (length(matches) > 0) {
      found_s10 <- TRUE
      logfc <- d[matches[1], 4]
      cat(sprintf("  S10%s (%s): FOUND, log2FC=%.2f\n", LETTERS[i], cells[i], logfc))
    }
  }
  if (!found_s10) cat("  S10: *** NOT FOUND IN ANY CELL TYPE ***\n")
  cat("\n")
}

# Also check Figure 6D proteins are in S8 (O-GalNAc protein DE, Jurkat)
cat("--- Figure 6D proteins (S8C = O-GalNAc protein DE Jurkat) ---\n\n")
fig6d_proteins <- c("SEL1L", "EDEM1", "OS9", "CPD", "MAN1A2", "XXYLT1", "TFRC", "IGSF8", "PTPRC")

d8c <- read.xlsx(paste0(table_dir, "supporting_table_S8.xlsx"), sheet=3, startRow=2, colNames=TRUE)
for (gene in fig6d_proteins) {
  matches <- which(d8c[,2] == gene)
  if (length(matches) > 0) {
    logfc <- d8c[matches[1], 3]
    pval <- d8c[matches[1], 4]
    cat(sprintf("  %s: FOUND in S8C, log2FC=%.2f, adj.P=%.2e\n", gene, logfc, pval))
  } else {
    cat(sprintf("  %s: *** NOT FOUND in S8C ***\n", gene))
  }
}
