library(readxl)
library(clusterProfiler)
library(org.Hs.eg.db)

args <- commandArgs(trailingOnly = TRUE)

input_xlsx <- if (length(args) >= 1) {
  args[[1]]
} else {
  "/Volumes/cos-lab-rwu60/Longping/Weekly_Reading/Mar26_2026/ac1c05438_si_004.xlsx"
}

output_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(getwd(), "..", "Manuscript", "ac1c05438_si_004_GO_analysis")
}

site_mode <- if (length(args) >= 3) {
  args[[3]]
} else {
  "all"
}

if (!site_mode %in% c("all", "st_only")) {
  stop("site_mode must be one of: all, st_only")
}

if (length(args) < 2 && identical(site_mode, "st_only")) {
  output_dir <- paste0(output_dir, "_ST_only")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

sheet_names <- excel_sheets(input_xlsx)

extract_cell_type <- function(sheet_name) {
  if (startsWith(sheet_name, "Jurkat")) {
    return("Jurkat")
  }
  if (startsWith(sheet_name, "MCF7")) {
    return("MCF7")
  }
  if (startsWith(sheet_name, "K562")) {
    return("K562")
  }
  NA_character_
}

extract_gene_symbol <- function(protein_name) {
  if (is.na(protein_name) || !nzchar(protein_name)) {
    return(NA_character_)
  }
  match <- regexpr("GN=([^ ]+)", protein_name, perl = TRUE)
  if (match[1] == -1) {
    return(NA_character_)
  }
  sub("GN=([^ ]+)", "\\1", regmatches(protein_name, match))
}

extract_accession <- function(protein_name) {
  if (is.na(protein_name) || !nzchar(protein_name)) {
    return(NA_character_)
  }
  parts <- strsplit(protein_name, "\\|")[[1]]
  if (length(parts) >= 3) {
    return(parts[2])
  }
  protein_name
}

extract_site_support <- function(mod_string) {
  if (is.na(mod_string) || !nzchar(mod_string)) {
    return(c(has_st = FALSE, has_c = FALSE))
  }
  has_st <- grepl("S[0-9]+\\(ogly\\)|T[0-9]+\\(ogly\\)", mod_string, perl = TRUE)
  has_c <- grepl("C[0-9]+\\(ogly\\)", mod_string, perl = TRUE)
  c(has_st = has_st, has_c = has_c)
}

cell_gene_lists <- list()
cell_accession_lists <- list()
count_summary <- data.frame()

for (sheet_name in sheet_names) {
  cell_type <- extract_cell_type(sheet_name)
  if (is.na(cell_type)) {
    next
  }

  sheet_df <- read_excel(input_xlsx, sheet = sheet_name, skip = 1)
  protein_col <- grep("Protein.*Name", names(sheet_df), value = TRUE)
  mods_col <- grep("Mods", names(sheet_df), value = TRUE)
  if (length(protein_col) != 1) {
    stop(sprintf("Could not uniquely identify protein column in sheet: %s", sheet_name))
  }
  if (length(mods_col) != 1) {
    stop(sprintf("Could not uniquely identify mods column in sheet: %s", sheet_name))
  }

  protein_names <- as.character(sheet_df[[protein_col]])
  mod_strings <- as.character(sheet_df[[mods_col]])
  is_contaminant <- grepl("Common contaminant protein", protein_names, fixed = TRUE)
  accessions <- vapply(protein_names, extract_accession, character(1))
  gene_symbols <- vapply(protein_names, extract_gene_symbol, character(1))
  site_support <- t(vapply(mod_strings, extract_site_support, logical(2)))

  keep_rows <- !is_contaminant
  if (identical(site_mode, "st_only")) {
    keep_rows <- keep_rows & site_support[, "has_st"]
  }

  clean_gene_symbols <- unique(gene_symbols[!is.na(gene_symbols) & keep_rows])
  clean_accessions <- unique(accessions[!is.na(accessions) & keep_rows])
  cell_gene_lists[[cell_type]] <- unique(c(cell_gene_lists[[cell_type]], clean_gene_symbols))
  cell_accession_lists[[cell_type]] <- unique(c(cell_accession_lists[[cell_type]], clean_accessions))

  count_summary <- rbind(
    count_summary,
    data.frame(
      site_mode = site_mode,
      cell_type = cell_type,
      sheet = sheet_name,
      unique_proteins_including_contaminants = length(unique(accessions[!is.na(accessions)])),
      unique_proteins_excluding_contaminants = length(unique(accessions[!is.na(accessions) & keep_rows])),
      unique_genes_excluding_contaminants = length(clean_gene_symbols),
      stringsAsFactors = FALSE
    )
  )
}

combined_summary <- do.call(
  rbind,
  lapply(names(cell_gene_lists), function(cell_type) {
    data.frame(
      site_mode = site_mode,
      cell_type = cell_type,
      sheet = "combined",
      unique_proteins_including_contaminants = NA_integer_,
      unique_proteins_excluding_contaminants = length(unique(cell_accession_lists[[cell_type]])),
      unique_genes_excluding_contaminants = length(unique(cell_gene_lists[[cell_type]])),
      stringsAsFactors = FALSE
    )
  })
)

count_summary <- rbind(count_summary, combined_summary)
write.csv(count_summary, file.path(output_dir, "protein_and_gene_counts.csv"), row.names = FALSE)

gene_table <- do.call(
  rbind,
  lapply(names(cell_gene_lists), function(cell_type) {
    data.frame(
      site_mode = site_mode,
      cell_type = cell_type,
      gene_symbol = sort(unique(cell_gene_lists[[cell_type]])),
      stringsAsFactors = FALSE
    )
  })
)
write.csv(gene_table, file.path(output_dir, "celltype_gene_lists.csv"), row.names = FALSE)

all_entrez <- keys(org.Hs.eg.db, keytype = "ENTREZID")
curated_map <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = all_entrez,
  columns = c("SYMBOL", "UNIPROT"),
  keytype = "ENTREZID"
)
curated_map <- unique(curated_map[!is.na(curated_map$UNIPROT), c("ENTREZID", "SYMBOL", "UNIPROT")])
curated_universe <- unique(curated_map$ENTREZID)

write.csv(
  data.frame(
    site_mode = site_mode,
    universe_name = "Human genes with UniProt protein mapping from org.Hs.eg.db",
    universe_entrez_count = length(curated_universe),
    universe_symbol_count = length(unique(curated_map$SYMBOL[!is.na(curated_map$SYMBOL)])),
    stringsAsFactors = FALSE
  ),
  file.path(output_dir, "background_universe_summary.csv"),
  row.names = FALSE
)

for (cell_type in names(cell_gene_lists)) {
  mapped <- bitr(
    cell_gene_lists[[cell_type]],
    fromType = "SYMBOL",
    toType = "ENTREZID",
    OrgDb = org.Hs.eg.db
  )

  mapping_summary <- data.frame(
    site_mode = site_mode,
    cell_type = cell_type,
    input_gene_count = length(cell_gene_lists[[cell_type]]),
    mapped_gene_count = nrow(mapped),
    unmapped_gene_count = length(setdiff(cell_gene_lists[[cell_type]], mapped$SYMBOL)),
    unmapped_genes = paste(sort(setdiff(cell_gene_lists[[cell_type]], mapped$SYMBOL)), collapse = ";"),
    stringsAsFactors = FALSE
  )
  write.csv(
    mapping_summary,
    file.path(output_dir, sprintf("%s_mapping_summary.csv", cell_type)),
    row.names = FALSE
  )

  for (ont in c("BP", "CC", "MF")) {
    ego <- suppressMessages(
      enrichGO(
        gene = unique(mapped$ENTREZID),
        universe = curated_universe,
        OrgDb = org.Hs.eg.db,
        keyType = "ENTREZID",
        ont = ont,
        pAdjustMethod = "BH",
        pvalueCutoff = 0.05,
        qvalueCutoff = 0.2,
        readable = TRUE
      )
    )

    out_path <- file.path(output_dir, sprintf("%s_GO_%s.csv", cell_type, ont))
    if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
      write.csv(
        data.frame(
          Description = character(),
          GeneRatio = character(),
          BgRatio = character(),
          pvalue = numeric(),
          p.adjust = numeric(),
          qvalue = numeric(),
          geneID = character(),
          Count = integer()
        ),
        out_path,
        row.names = FALSE
      )
    } else {
      write.csv(as.data.frame(ego), out_path, row.names = FALSE)
    }
  }
}

message("GO analysis written to: ", normalizePath(output_dir))
