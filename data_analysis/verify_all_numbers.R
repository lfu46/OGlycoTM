library(openxlsx)
table_dir <- "/Users/longpingfu/Downloads/OGlycoTM/Manuscript/version_Mar24_2026/"
cells <- c("HEK293T", "HepG2", "Jurkat")

parse_percent_column <- function(x) {
  if (is.numeric(x)) return(as.numeric(x))
  suppressWarnings(as.numeric(gsub("%", "", x)))
}

check_protein_id_coverage <- function(table_path, table_label) {
  cat(sprintf("\n%s coverage sanity check:\n", table_label))
  for (i in 1:3) {
    d <- read.xlsx(table_path, sheet = i, startRow = 2, colNames = TRUE)
    cov <- suppressWarnings(as.numeric(d$Coverage))
    len <- suppressWarnings(as.numeric(d$Length))
    pct <- parse_percent_column(d$Coverage.Percentage)
    derived_pct <- ifelse(!is.na(len) & len > 0, 100 * cov / len, NA_real_)

    bad_pct_range <- sum(pct > 100 | pct < 0, na.rm = TRUE)
    bad_cov_len <- sum(cov > len, na.rm = TRUE)
    bad_consistency <- sum(abs(pct - derived_pct) > 1, na.rm = TRUE)

    cat(sprintf(
      "   %s: pct_out_of_range=%d, coverage_gt_length=%d, pct_vs_length_mismatch=%d\n",
      cells[i], bad_pct_range, bad_cov_len, bad_consistency
    ))
  }
}

cat("========== SUPPORTING TABLE INVENTORY ==========\n\n")
for (s in 1:11) {
  f <- paste0(table_dir, "supporting_table_S", s, ".xlsx")
  wb <- loadWorkbook(f)
  d1 <- read.xlsx(wb, sheet = 1, startRow = 2, colNames = TRUE)
  for (j in seq_along(names(wb))) {
    d <- read.xlsx(wb, sheet = j, colNames = FALSE)
    cat(sprintf("S%-2d %s: %4d rows | %s\n", s, LETTERS[j], nrow(d)-2,
      substr(as.character(d[1,1]), 1, 100)))
  }
}

cat("\n========== MANUSCRIPT NUMBER CHECKS ==========\n\n")

# 1. "1,109 O-GlcNAcylated proteins" (S1)
all_s1 <- c()
for (i in 1:3) {
  d <- read.xlsx(paste0(table_dir, "supporting_table_S1.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  all_s1 <- c(all_s1, d[,1])
}
n_unique_s1 <- length(unique(all_s1))
cat(sprintf("1. Total O-GlcNAc proteins: manuscript=1109, table=%d %s\n",
  n_unique_s1, ifelse(n_unique_s1==1109, "OK", "MISMATCH")))

# 2. Per-cell: 775/676/692 (S1)
for (i in 1:3) {
  d <- read.xlsx(paste0(table_dir, "supporting_table_S1.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  expected <- c(775, 676, 692)[i]
  cat(sprintf("2. S1 %s: manuscript=%d, table=%d %s\n",
    cells[i], expected, nrow(d), ifelse(nrow(d)==expected, "OK", "MISMATCH")))
}

# 3. "402 commonly identified" (need overlap of all 3 S1 sheets)
s1_prots <- list()
for (i in 1:3) {
  d <- read.xlsx(paste0(table_dir, "supporting_table_S1.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  s1_prots[[i]] <- d[,1]
}
common <- Reduce(intersect, s1_prots)
cat(sprintf("3. Common O-GlcNAc proteins: manuscript=402, table=%d %s\n",
  length(common), ifelse(length(common)==402, "OK", "MISMATCH")))

# 4. Up/down regulated: 163/113/62 up, 146/96/86 down (S2)
cat("\n4. Regulated O-GlcNAc proteins (S2):\n")
for (i in 1:3) {
  d <- read.xlsx(paste0(table_dir, "supporting_table_S2.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  logfc_col <- 3; pval_col <- 4
  n_up <- sum(d[,logfc_col] > 0.5 & d[,pval_col] < 0.05, na.rm=TRUE)
  n_down <- sum(d[,logfc_col] < -0.5 & d[,pval_col] < 0.05, na.rm=TRUE)
  exp_up <- c(163, 113, 62)[i]; exp_down <- c(146, 96, 86)[i]
  cat(sprintf("   %s: up=%d (exp %d %s), down=%d (exp %d %s)\n",
    cells[i], n_up, exp_up, ifelse(n_up==exp_up,"OK","MISMATCH"),
    n_down, exp_down, ifelse(n_down==exp_down,"OK","MISMATCH")))
}

# 5. "485 O-GalNAcylated proteins" (S7)
all_s7 <- c()
for (i in 1:3) {
  d <- read.xlsx(paste0(table_dir, "supporting_table_S7.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  all_s7 <- c(all_s7, d[,1])
}
n_unique_s7 <- length(unique(all_s7))
cat(sprintf("\n5. Total O-GalNAc proteins: manuscript=485, table=%d %s\n",
  n_unique_s7, ifelse(n_unique_s7==485, "OK", "MISMATCH")))

# 6. S-GlcNAc Cys sites S11 counts
cat("\n6. S-GlcNAc Cys sites (S11):\n")
for (i in 1:3) {
  d <- read.xlsx(paste0(table_dir, "supporting_table_S11.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  cat(sprintf("   %s: %d sites\n", cells[i], nrow(d)))
}

# 7. P2P Comment 3: proteins in correct tables
cat("\n7. P2P Comment 3 protein verification:\n")
checks <- list(
  list(id="Q9Y5M8", gene="SRPRB", t=c(1,2), s=c(1,1)),
  list(id="O00161", gene="SNAP23", t=c(1,2), s=c(1,1)),
  list(id="P11166", gene="SLC2A1", t=c(1,2), s=c(1,1)),
  list(id="Q92667", gene="AKAP1", t=c(1,2), s=c(2,2)),
  list(id="Q2M2I8", gene="AAK1", t=c(1,2), s=c(1,1))
)
for (chk in checks) {
  for (k in seq_along(chk$t)) {
    d <- read.xlsx(paste0(table_dir, "supporting_table_S", chk$t[k], ".xlsx"),
                   sheet=chk$s[k], startRow=2, colNames=TRUE)
    found <- chk$id %in% d[,1]
    cat(sprintf("   %s (%s) in S%d%s: %s\n",
      chk$gene, chk$id, chk$t[k], LETTERS[chk$s[k]], ifelse(found,"FOUND","MISSING")))
  }
}

# 8. Cross-table consistency: S8 subset of S7
cat("\n8. S8 proteins subset of S7:\n")
for (i in 1:3) {
  d7 <- read.xlsx(paste0(table_dir, "supporting_table_S7.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  d8 <- read.xlsx(paste0(table_dir, "supporting_table_S8.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  only_s8 <- sum(!d8[,1] %in% d7[,1])
  cat(sprintf("   %s: S7=%d, S8=%d, S8-only=%d %s\n",
    cells[i], nrow(d7), nrow(d8), only_s8, ifelse(only_s8==0,"OK","MISMATCH")))
}

# 9. S2 subset of S1
cat("\n9. S2 proteins subset of S1:\n")
for (i in 1:3) {
  d1 <- read.xlsx(paste0(table_dir, "supporting_table_S1.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  d2 <- read.xlsx(paste0(table_dir, "supporting_table_S2.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  only_s2 <- sum(!d2[,1] %in% d1[,1])
  cat(sprintf("   %s: S1=%d, S2=%d, S2-only=%d %s\n",
    cells[i], nrow(d1), nrow(d2), only_s2, ifelse(only_s2==0,"OK","MISMATCH")))
}

# 10. S4 matches S3 (same proteins)
cat("\n10. S4 matches S3 protein count:\n")
for (i in 1:3) {
  d3 <- read.xlsx(paste0(table_dir, "supporting_table_S3.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  d4 <- read.xlsx(paste0(table_dir, "supporting_table_S4.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  cat(sprintf("   %s: S3=%d, S4=%d %s\n",
    cells[i], nrow(d3), nrow(d4), ifelse(nrow(d3)==nrow(d4),"OK","MISMATCH")))
}

# 11. Numeric type check
cat("\n11. Numeric column type check:\n")
all_ok <- TRUE
for (s in 1:11) {
  f <- paste0(table_dir, "supporting_table_S", s, ".xlsx")
  d <- read.xlsx(f, sheet=1, startRow=2, colNames=TRUE)
  nc <- ncol(d)
  issues <- c()
  if (nc == 8) {
    for (j in 3:6) if (!is.numeric(d[1,j])) issues <- c(issues, colnames(d)[j])
  } else if (nc == 5) {
    for (j in 3:4) if (!is.numeric(d[1,j])) issues <- c(issues, colnames(d)[j])
  } else if (nc == 13) {
    for (j in c(4,7,8,9,10,12)) if (!is.numeric(d[1,j])) issues <- c(issues, colnames(d)[j])
  } else if (nc == 6) {
    for (j in 4:5) if (!is.numeric(d[1,j])) issues <- c(issues, colnames(d)[j])
  }
  status <- if (length(issues)==0) "OK" else paste("TEXT:", paste(issues, collapse=", "))
  cat(sprintf("   S%-2d: %s\n", s, status))
  if (length(issues) > 0) all_ok <- FALSE
}
if (all_ok) cat("   All numeric columns OK\n")

# 12. Column header consistency
cat("\n12. Column header consistency (paired tables):\n")
# S5 vs S9 (site ID)
d5 <- read.xlsx(paste0(table_dir, "supporting_table_S5.xlsx"), sheet=1, colNames=FALSE)
d9 <- read.xlsx(paste0(table_dir, "supporting_table_S9.xlsx"), sheet=1, colNames=FALSE)
d11 <- read.xlsx(paste0(table_dir, "supporting_table_S11.xlsx"), sheet=1, colNames=FALSE)
h5 <- paste(d5[2,], collapse="|"); h9 <- paste(d9[2,], collapse="|"); h11 <- paste(d11[2,], collapse="|")
cat(sprintf("   S5 vs S9:  %s\n", ifelse(h5==h9, "MATCH", "MISMATCH")))
cat(sprintf("   S5 vs S11: %s\n", ifelse(h5==h11, "MATCH", "MISMATCH")))

# S6 vs S10 (site DE)
d6 <- read.xlsx(paste0(table_dir, "supporting_table_S6.xlsx"), sheet=1, colNames=FALSE)
d10 <- read.xlsx(paste0(table_dir, "supporting_table_S10.xlsx"), sheet=1, colNames=FALSE)
h6 <- paste(d6[2,], collapse="|"); h10 <- paste(d10[2,], collapse="|")
cat(sprintf("   S6 vs S10: %s\n", ifelse(h6==h10, "MATCH", "MISMATCH")))

# S2 vs S8 (protein DE)
d2 <- read.xlsx(paste0(table_dir, "supporting_table_S2.xlsx"), sheet=1, colNames=FALSE)
d8 <- read.xlsx(paste0(table_dir, "supporting_table_S8.xlsx"), sheet=1, colNames=FALSE)
h2 <- paste(d2[2,], collapse="|"); h8 <- paste(d8[2,], collapse="|")
cat(sprintf("   S2 vs S8:  %s\n", ifelse(h2==h8, "MATCH", "MISMATCH")))

# 13. Coverage sanity checks for protein ID tables
cat("\n13. Coverage sanity checks:\n")
check_protein_id_coverage(paste0(table_dir, "supporting_table_S1.xlsx"), "   S1")
check_protein_id_coverage(paste0(table_dir, "supporting_table_S7.xlsx"), "   S7")
