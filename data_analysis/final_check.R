library(openxlsx)
table_dir <- "/Users/longpingfu/Downloads/OGlycoTM/Manuscript/version_Mar24_2026/"
cells <- c("HEK293T", "HepG2", "Jurkat")

ms <- readLines("/tmp/ms_final.txt")
p2p <- readLines("/tmp/p2p_final.txt")
si <- readLines("/tmp/si_final.txt")

errors <- c()
add_err <- function(msg) { errors <<- c(errors, msg) }

cat("================================================================\n")
cat("  COMPREHENSIVE PRE-SUBMISSION CHECK\n")
cat("================================================================\n\n")

# ======================================================================
# 1. TABLE INVENTORY & NUMERIC TYPES
# ======================================================================
cat("--- 1. Supporting Table Inventory ---\n\n")

for (s in 1:11) {
  f <- paste0(table_dir, "supporting_table_S", s, ".xlsx")
  wb <- loadWorkbook(f)
  for (j in seq_along(names(wb))) {
    d <- read.xlsx(wb, sheet = j, colNames = FALSE)
    n <- nrow(d) - 2
    title <- substr(as.character(d[1,1]), 1, 110)
    # Check title numbering
    expected_label <- paste0("Table S", s, LETTERS[j])
    if (!grepl(expected_label, title, fixed = TRUE)) {
      add_err(sprintf("S%d sheet %d: title says '%s' but expected '%s'",
        s, j, substr(title, 1, 15), expected_label))
    }
    cat(sprintf("  S%-2d%s: %4d rows  %s\n", s, LETTERS[j], n, title))
  }
  # Numeric type check
  d <- read.xlsx(f, sheet = 1, startRow = 2, colNames = TRUE)
  nc <- ncol(d)
  bad_cols <- c()
  check_cols <- if (nc == 8) 3:6 else if (nc == 5) 3:4 else if (nc == 13) c(4,7,8,9,10,12) else if (nc == 6) 4:5 else c()
  for (j in check_cols) {
    if (!is.numeric(d[1, j])) bad_cols <- c(bad_cols, colnames(d)[j])
  }
  if (length(bad_cols) > 0) add_err(sprintf("S%d: text in numeric columns: %s", s, paste(bad_cols, collapse=", ")))
}

# ======================================================================
# 2. MANUSCRIPT NUMBERS
# ======================================================================
cat("\n--- 2. Manuscript Number Verification ---\n\n")

# 2a. 1,109 O-GlcNAcylated proteins
all_s1 <- c(); for (i in 1:3) { d <- read.xlsx(paste0(table_dir,"supporting_table_S1.xlsx"), sheet=i, startRow=2, colNames=TRUE); all_s1 <- c(all_s1, d[,1]) }
n <- length(unique(all_s1))
cat(sprintf("  Total O-GlcNAc proteins: manuscript=1109, table=%d %s\n", n, ifelse(n==1109,"OK","*** MISMATCH ***")))
if (n != 1109) add_err(sprintf("O-GlcNAc total: manuscript=1109, table=%d", n))

# 2b. Per-cell 775/676/692
for (i in 1:3) {
  d <- read.xlsx(paste0(table_dir,"supporting_table_S1.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  exp <- c(775,676,692)[i]
  ok <- nrow(d)==exp
  cat(sprintf("  S1 %s: manuscript=%d, table=%d %s\n", cells[i], exp, nrow(d), ifelse(ok,"OK","*** MISMATCH ***")))
  if (!ok) add_err(sprintf("S1 %s: manuscript=%d, table=%d", cells[i], exp, nrow(d)))
}

# 2c. 402 common
s1p <- list(); for (i in 1:3) { d <- read.xlsx(paste0(table_dir,"supporting_table_S1.xlsx"), sheet=i, startRow=2, colNames=TRUE); s1p[[i]] <- d[,1] }
n_common <- length(Reduce(intersect, s1p))
cat(sprintf("  Common O-GlcNAc: manuscript=402, table=%d %s\n", n_common, ifelse(n_common==402,"OK","*** MISMATCH ***")))
if (n_common != 402) add_err(sprintf("Common O-GlcNAc: manuscript=402, table=%d", n_common))

# 2d. Up/down: 163/113/62 up, 146/96/86 down
cat("  Regulated (S2):\n")
exp_up <- c(163,113,62); exp_down <- c(146,96,86)
for (i in 1:3) {
  d <- read.xlsx(paste0(table_dir,"supporting_table_S2.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  n_up <- sum(d[,3] > 0.5 & d[,4] < 0.05, na.rm=TRUE)
  n_down <- sum(d[,3] < -0.5 & d[,4] < 0.05, na.rm=TRUE)
  ok_up <- n_up==exp_up[i]; ok_down <- n_down==exp_down[i]
  cat(sprintf("    %s: up=%d(%s) down=%d(%s)\n", cells[i],
    n_up, ifelse(ok_up,"OK","MISMATCH"), n_down, ifelse(ok_down,"OK","MISMATCH")))
  if (!ok_up) add_err(sprintf("S2 %s up: manuscript=%d, table=%d", cells[i], exp_up[i], n_up))
  if (!ok_down) add_err(sprintf("S2 %s down: manuscript=%d, table=%d", cells[i], exp_down[i], n_down))
}

# 2e. 485 O-GalNAc proteins
all_s7 <- c(); for (i in 1:3) { d <- read.xlsx(paste0(table_dir,"supporting_table_S7.xlsx"), sheet=i, startRow=2, colNames=TRUE); all_s7 <- c(all_s7, d[,1]) }
n <- length(unique(all_s7))
cat(sprintf("  Total O-GalNAc proteins: manuscript=485, table=%d %s\n", n, ifelse(n==485,"OK","*** MISMATCH ***")))
if (n != 485) add_err(sprintf("O-GalNAc total: manuscript=485, table=%d", n))

# ======================================================================
# 3. CROSS-TABLE CONSISTENCY
# ======================================================================
cat("\n--- 3. Cross-Table Consistency ---\n\n")

# S2 ⊆ S1
for (i in 1:3) {
  d1 <- read.xlsx(paste0(table_dir,"supporting_table_S1.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  d2 <- read.xlsx(paste0(table_dir,"supporting_table_S2.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  extra <- sum(!d2[,1] %in% d1[,1])
  cat(sprintf("  S2%s ⊆ S1%s: %s\n", LETTERS[i], LETTERS[i], ifelse(extra==0,"OK",paste(extra,"extra"))))
  if (extra > 0) add_err(sprintf("S2%s has %d proteins not in S1%s", LETTERS[i], extra, LETTERS[i]))
}

# S8 ⊆ S7
for (i in 1:3) {
  d7 <- read.xlsx(paste0(table_dir,"supporting_table_S7.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  d8 <- read.xlsx(paste0(table_dir,"supporting_table_S8.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  extra <- sum(!d8[,1] %in% d7[,1])
  cat(sprintf("  S8%s ⊆ S7%s: %s\n", LETTERS[i], LETTERS[i], ifelse(extra==0,"OK",paste(extra,"extra"))))
  if (extra > 0) add_err(sprintf("S8%s has %d proteins not in S7%s", LETTERS[i], extra, LETTERS[i]))
}

# S4 = S3 count
for (i in 1:3) {
  d3 <- read.xlsx(paste0(table_dir,"supporting_table_S3.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  d4 <- read.xlsx(paste0(table_dir,"supporting_table_S4.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  ok <- nrow(d3)==nrow(d4)
  cat(sprintf("  S3%s=S4%s count: %d vs %d %s\n", LETTERS[i], LETTERS[i], nrow(d3), nrow(d4), ifelse(ok,"OK","MISMATCH")))
  if (!ok) add_err(sprintf("S3%s=%d vs S4%s=%d", LETTERS[i], nrow(d3), LETTERS[i], nrow(d4)))
}

# S6 sites ⊆ S5 sites (by UniProt+Site)
for (i in 1:3) {
  d5 <- read.xlsx(paste0(table_dir,"supporting_table_S5.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  d6 <- read.xlsx(paste0(table_dir,"supporting_table_S6.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  s5_keys <- paste(d5[,1], d5[,3])
  s6_keys <- paste(d6[,1], d6[,3])
  extra <- sum(!s6_keys %in% s5_keys)
  cat(sprintf("  S6%s ⊆ S5%s: %s\n", LETTERS[i], LETTERS[i], ifelse(extra==0,"OK",paste(extra,"extra"))))
  if (extra > 0) add_err(sprintf("S6%s has %d sites not in S5%s", LETTERS[i], extra, LETTERS[i]))
}

# S10 sites ⊆ S9 sites
for (i in 1:3) {
  d9 <- read.xlsx(paste0(table_dir,"supporting_table_S9.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  d10 <- read.xlsx(paste0(table_dir,"supporting_table_S10.xlsx"), sheet=i, startRow=2, colNames=TRUE)
  s9_keys <- paste(d9[,1], d9[,3])
  s10_keys <- paste(d10[,1], d10[,3])
  extra <- sum(!s10_keys %in% s9_keys)
  cat(sprintf("  S10%s ⊆ S9%s: %s\n", LETTERS[i], LETTERS[i], ifelse(extra==0,"OK",paste(extra,"extra"))))
  if (extra > 0) add_err(sprintf("S10%s has %d sites not in S9%s", LETTERS[i], extra, LETTERS[i]))
}

# ======================================================================
# 4. COLUMN HEADER CONSISTENCY
# ======================================================================
cat("\n--- 4. Column Header Consistency ---\n\n")

get_headers <- function(s) {
  d <- read.xlsx(paste0(table_dir,"supporting_table_S",s,".xlsx"), sheet=1, colNames=FALSE)
  paste(d[2,], collapse="|")
}

pairs <- list(c(1,7,"Protein ID"), c(2,8,"Protein DE"), c(5,9,"Site ID"), c(5,11,"Site ID (Cys)"), c(6,10,"Site DE"))
for (p in pairs) {
  h1 <- get_headers(p[[1]]); h2 <- get_headers(p[[2]])
  ok <- h1 == h2
  cat(sprintf("  S%s vs S%s (%s): %s\n", p[[1]], p[[2]], p[[3]], ifelse(ok,"MATCH","*** MISMATCH ***")))
  if (!ok) add_err(sprintf("Headers S%s vs S%s mismatch", p[[1]], p[[2]]))
}

# ======================================================================
# 5. P2P RESPONSE CHECKS
# ======================================================================
cat("\n--- 5. P2P Response Verification ---\n\n")

# Comment 3 proteins
checks <- list(
  list(id="Q9Y5M8", gene="SRPRB", t=c(1,2), s=c(1,1)),
  list(id="O00161", gene="SNAP23", t=c(1,2), s=c(1,1)),
  list(id="P11166", gene="SLC2A1", t=c(1,2), s=c(1,1)),
  list(id="Q92667", gene="AKAP1", t=c(1,2), s=c(2,2)),
  list(id="Q2M2I8", gene="AAK1", t=c(1,2), s=c(1,1))
)
for (chk in checks) {
  for (k in seq_along(chk$t)) {
    d <- read.xlsx(paste0(table_dir,"supporting_table_S",chk$t[k],".xlsx"), sheet=chk$s[k], startRow=2, colNames=TRUE)
    found <- chk$id %in% d[,1]
    cat(sprintf("  %s in S%d%s: %s\n", chk$gene, chk$t[k], LETTERS[chk$s[k]], ifelse(found,"FOUND","*** MISSING ***")))
    if (!found) add_err(sprintf("%s (%s) missing from S%d%s", chk$gene, chk$id, chk$t[k], LETTERS[chk$s[k]]))
  }
}

# P2P says "Table S7 and Table S8 ... Table S9 and Table S10" for O-GalNAc
p2p_text <- paste(p2p, collapse = " ")
if (grepl("Table S7 and Table S8", p2p_text)) cat("  P2P S7/S8 reference: OK\n") else add_err("P2P missing S7/S8 reference")
if (grepl("Table S9 and Table S10", p2p_text)) cat("  P2P S9/S10 reference: OK\n") else add_err("P2P missing S9/S10 reference")
if (grepl("Table S11", p2p_text)) cat("  P2P S11 reference: OK\n") else add_err("P2P missing S11 reference")

# ======================================================================
# 6. SI DOCUMENT TABLE LIST
# ======================================================================
cat("\n--- 6. SI Document Table List ---\n\n")

si_text <- paste(si, collapse = " ")
expected_si <- list(
  "S1"  = "Identification of O-GlcNAcylated proteins",
  "S2"  = "Abundance changes of O-GlcNAcylated proteins",
  "S3"  = "Identification of total proteins",
  "S4"  = "Abundance changes of total proteins",
  "S5"  = "Identification of O-GlcNAcylation sites",
  "S6"  = "Abundance changes of O-GlcNAcylation sites",
  "S7"  = "Identification of O-GalNAcylated proteins",
  "S8"  = "Abundance changes of O-GalNAcylated proteins",
  "S9"  = "Identification of O-GalNAcylation sites",
  "S10" = "Identification of S-GlcNAcylation",  # S11 in tables
  "S11" = "Abundance changes of O-GalNAcylation sites"  # S10 in tables
)

for (s in paste0("Table S", 1:11)) {
  found <- grepl(s, si_text)
  cat(sprintf("  %s in SI: %s\n", s, ifelse(found, "FOUND", "*** MISSING ***")))
  if (!found) add_err(sprintf("%s not listed in SI document", s))
}

# ======================================================================
# 7. MANUSCRIPT ↔ SI TABLE REFERENCE CONSISTENCY
# ======================================================================
cat("\n--- 7. Manuscript Table References ---\n\n")

ms_text <- paste(ms, collapse = " ")
ms_tables <- unique(regmatches(ms_text, gregexpr("Table S\\d+", ms_text))[[1]])
cat("  Tables referenced in manuscript:", paste(sort(ms_tables), collapse=", "), "\n")

expected_ms_tables <- paste0("Table S", c(1,2,5,6,7,8,9,10,11))
for (t in expected_ms_tables) {
  found <- t %in% ms_tables
  if (!found) {
    cat(sprintf("  WARNING: %s not referenced in manuscript\n", t))
    # Not an error for S3/S4 (WP) which may not be cited in text
  }
}

# Check S3/S4 are at least in SI
cat("  S3 in SI listing:", ifelse(grepl("Table S3", si_text), "YES", "NO"), "\n")
cat("  S4 in SI listing:", ifelse(grepl("Table S4", si_text), "YES", "NO"), "\n")

# ======================================================================
# 8. FIGURE NUMBER CONSISTENCY
# ======================================================================
cat("\n--- 8. Figure References ---\n\n")

ms_figs <- unique(regmatches(ms_text, gregexpr("Figure [0-9S]+[A-F]?", ms_text))[[1]])
cat("  Figures referenced:", paste(sort(unique(gsub("[A-F]$","",ms_figs))), collapse=", "), "\n")

# Check no duplicate Figure captions
fig_captions <- grep("^Figure [0-9]\\.", ms, value = TRUE)
cat("  Figure captions found:", length(fig_captions), "\n")
for (fc in fig_captions) cat("    ", substr(fc, 1, 80), "\n")

# ======================================================================
# FINAL SUMMARY
# ======================================================================
cat("\n================================================================\n")
if (length(errors) == 0) {
  cat("  ALL CHECKS PASSED - Ready for advisor review\n")
} else {
  cat(sprintf("  %d ISSUE(S) FOUND:\n", length(errors)))
  for (e in errors) cat("  - ", e, "\n")
}
cat("================================================================\n")
