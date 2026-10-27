library(shiny)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)
library(plotly)
library(RColorBrewer)
library(data.table)
library(DBI)
library(RSQLite)
library(clusterProfiler)
library(org.Hs.eg.db)
library(uwot)
library(bslib)
library(circlize)
select <- dplyr::select
filter <- dplyr::filter
rename <- dplyr::rename
slice  <- dplyr::slice
TROPICAL    <- c("#ffffd9", "#c7e9b4", "#41b6c4", "#225ea8", "#081d58")
BERRY_SIG   <- "#c51b8a"
BERRY_LIGHT <- "#fcc5c0"
GREY_NONSIG <- "#bdbdbd"

tissue_colors <- function(tissues) {
  base <- c(
    "#e6194b",
    "#3cb44b",
    "#6A5ACD",
    "#f58231",
    "#911eb4",
    "#42d4f4",
    "#f032e6",
    "#bfef45",
    "#fabed4",
    "#469990",
    "#dcbeff",
    "#1E88E5",
    "#E53935",
    "#aaffc3",
    "#C0CA33",
    "#5C6BC0",
    "#FFC107",
    "#C71585",
    "#ffe119",
    "#ffd8b1",
    "#26A69A",
    "#FF4500",
    "#ff6db6",
    "#1abc9c",
    "#AB47BC",
    "#ff1493",
    "#6A5ACD",
    "#556b2f",
    "#2f4f4f",
    "#ff1493"
  )
  n   <- length(tissues)
  pal <- if (n <= length(base))
    base[seq_len(n)]
  else
    colorRampPalette(base)(n)
  setNames(pal, sort(tissues))
}

collapse_isoforms <- function(df_long, n = 5) {
  top <- df_long %>%
    group_by(display_label) %>%
    summarise(total = sum(TPM, na.rm = TRUE), .groups = "drop") %>%
    arrange(desc(total)) %>% head(n) %>% pull(display_label)
  df_long %>% mutate(isoform_group = ifelse(display_label %in% top, display_label, "Other"))
}

get_colors <- function(n) {
  base8 <- c(
    "#66C2A5",
    "#FC8D62",
    "#FFD92F",
    "#E78AC3",
    "#A6D854",
    "#8DA0CB",
    "#E5C494",
    "#B3B3B3"
  )
  if (n <= 8)
    base8[1:n]
  else
    colorRampPalette(base8)(n)
}

shannon_entropy <- function(x) {
  x <- x[x > 0]
  if (length(x) == 0)
    return(0)
  n <- length(x)
  if (n == 1)
    return(0)
  p <- x / sum(x)
  raw <- -sum(p * log2(p))
  raw / log2(n)
}

symbols_to_entrez <- function(symbols) {
  mapped <- AnnotationDbi::mapIds(
    org.Hs.eg.db,
    keys = symbols,
    column = "ENTREZID",
    keytype = "SYMBOL",
    multiVals = "first"
  )
  mapped[!is.na(mapped)]
}

DATA_ROOT <- "shiny_data"

read_dataset_meta <- function(id) {
  p <- file.path(DATA_ROOT, id, "quants.db")
  if (!file.exists(p)) return(NULL)
  con <- dbConnect(RSQLite::SQLite(), p)
  on.exit(dbDisconnect(con))
  if (!"dataset_meta" %in% dbListTables(con)) return(NULL)
  m <- dbReadTable(con, "dataset_meta")
  m$db_path <- p
  m$deseq_path <- file.path(DATA_ROOT, id, paste0("deseq2_", id, ".db"))
  m$has_deseq <- file.exists(m$deseq_path)
  m
}

DATASETS <- do.call(rbind, Filter(Negate(is.null),
  lapply(list.dirs(DATA_ROOT, recursive = FALSE, full.names = FALSE), read_dataset_meta)))

if (is.null(DATASETS) || nrow(DATASETS) == 0)
  stop("No datasets found under ", DATA_ROOT)

REF_ID <- DATASETS$dataset_id[DATASETS$is_reference == 1][1]
if (is.na(REF_ID)) stop("No reference dataset")
QUERY_IDS <- DATASETS$dataset_id[DATASETS$is_reference == 0]

ref_row <- DATASETS[DATASETS$dataset_id == REF_ID, ]
cat("Reference:", REF_ID, "|", ref_row$n_samples, "samples\n")
cat("Query datasets:", if (length(QUERY_IDS)) paste(QUERY_IDS, collapse = ", ") else "none", "\n")

db <- dbConnect(RSQLite::SQLite(), ref_row$db_path)

deseq_db_path <- ref_row$deseq_path
HAS_DESEQ2 <- ref_row$has_deseq
deseq_db <- if (HAS_DESEQ2) {
  cat("Connecting to DESeq2 database...\n")
  dbConnect(RSQLite::SQLite(), deseq_db_path)
} else {
  cat("DESeq2 database not found. Wilcoxon-only mode.\n")
  NULL
}

attached_query <- NULL

attach_query <- function(id) {
  if (!is.null(attached_query)) {
    try(dbExecute(db, "DETACH DATABASE query"), silent = TRUE)
    attached_query <<- NULL
  }
  if (is.null(id) || id == REF_ID) return(invisible(NULL))
  row <- DATASETS[DATASETS$dataset_id == id, ]
  if (nrow(row) != 1) return(invisible(NULL))
  dbExecute(db, sprintf("ATTACH DATABASE '%s' AS query", row$db_path))
  attached_query <<- id
  invisible(id)
}

ensure_db <- function() {
  if (!dbIsValid(db)) {
    db <<- dbConnect(RSQLite::SQLite(), ref_row$db_path)
    attached_query <<- NULL
  }
  if (HAS_DESEQ2 && !is.null(deseq_db) && !dbIsValid(deseq_db)) {
    deseq_db <<- dbConnect(RSQLite::SQLite(), deseq_db_path)
  }
}

cat("Loading precomputed tables...\n")
sample_meta_df       <- dbReadTable(db, "samples")
sample_meta_df$Tissue   <- sample_meta_df$group_label
sample_meta_df$CellLine <- sample_meta_df$subject_id
tissue_avg_df        <- as.data.frame(fread(file.path(DATA_ROOT, REF_ID, "group_avg.csv")))

tx_annot_path <- "shiny_data/transcript_annotations.csv"
protein_annot_path <- "shiny_data/protein_annotations_collapsed.csv"

HAS_TX_ANNOT <- file.exists(tx_annot_path)
HAS_PROTEIN_ANNOT <- file.exists(protein_annot_path)

tx_annot_df <- if (HAS_TX_ANNOT) {
  cat("Loading transcript annotations...\n")
  as.data.frame(fread(tx_annot_path))
} else {
  NULL
}

protein_annot_df <- if (HAS_PROTEIN_ANNOT) {
  cat("Loading protein annotations...\n")
  as.data.frame(fread(protein_annot_path))
} else {
  NULL
}

tx_exons_path <- "shiny_data/transcript_exons.csv.gz"
HAS_TX_EXONS <- file.exists(tx_exons_path)
tx_exons_df <- if (HAS_TX_EXONS) {
  cat("Loading transcript exon structure...\n")
  as.data.frame(fread(tx_exons_path))
} else {
  NULL
}

join_annotations <- function(df) {
  if (!"transcript_id" %in% names(df))
    return(df)
  df$tx_id_clean <- df$transcript_id
  if (!is.null(tx_annot_df)) {
    df <- df %>% left_join(
      tx_annot_df %>%
        select(transcript_id, biotype, n_exons, tx_length) %>%
        rename(
          Biotype = biotype,
          Exons = n_exons,
          `Transcript Length (bp)` = tx_length
        ),
      by = c("tx_id_clean" = "transcript_id")
    )
  }
  if (HAS_REFSEQ && !is.null(refseq_map_df)) {
    df <- df %>% left_join(
      refseq_map_df %>%
        rename(RefSeq = refseq),
      by = c("tx_id_clean" = "ensembl_transcript_id")
    )
    df$RefSeq[is.na(df$RefSeq)] <- ""
  }
  if (!is.null(protein_annot_df)) {
    df <- df %>% left_join(
      protein_annot_df %>%
        select(
          transcript_id,
          cds_length,
          uniprot,
          pfam_domains,
          interpro_domains
        ) %>%
        rename(
          `CDS (bp)` = cds_length,
          UniProt = uniprot,
          Pfam = pfam_domains,
          InterPro = interpro_domains
        ),
      by = c("tx_id_clean" = "transcript_id")
    )
    if ("Biotype" %in% names(df)) {
      non_coding <- !is.na(df$Biotype) & df$Biotype != "protein_coding"
      if ("UniProt" %in% names(df))
        df$UniProt[non_coding &
                     is.na(df$UniProt)] <- "N/A (non-coding)"
      if ("Pfam" %in% names(df))
        df$Pfam[non_coding & is.na(df$Pfam)] <- "N/A (non-coding)"
      if ("InterPro" %in% names(df))
        df$InterPro[non_coding &
                      is.na(df$InterPro)] <- "N/A (non-coding)"
      if ("CDS (bp)" %in% names(df))
        df$`CDS (bp)`[non_coding & is.na(df$`CDS (bp)`)] <- NA
    }
  }
  df$tx_id_clean <- NULL
  if ("RefSeq" %in% names(df) && "transcript_id" %in% names(df)) {
    tx_pos <- which(names(df) == "transcript_id")
    other_cols <- setdiff(names(df), "RefSeq")
    new_order <- append(other_cols, "RefSeq", after = tx_pos)
    df <- df[, new_order, drop = FALSE]
  }
  df
}

tissue_avg_df$display_label <- paste0(tissue_avg_df$gene, " | ", tissue_avg_df$transcript_id)

get_active_avg <- function(id) {
  if (is.null(id) || id == REF_ID) return(tissue_avg_df)
  if (!is.null(AVG_CACHE[[id]])) return(AVG_CACHE[[id]])
  p <- file.path(DATA_ROOT, id, "group_avg.csv")
  if (!file.exists(p)) return(tissue_avg_df)
  q <- fread(p)
  gcols <- setdiff(names(q), c("gene", "transcript_id", "display_label"))
  setnames(q, gcols, pretty_group(gcols, id))
  q <- q[, setdiff(names(q), c("display_label", names(tissue_avg_df)[!names(tissue_avg_df) %in% c("gene", "transcript_id")])), with = FALSE]
  m <- merge(as.data.table(tissue_avg_df)[, display_label := NULL], q,
             by = c("gene", "transcript_id"), all = TRUE)
  num <- names(m)[vapply(m, is.numeric, logical(1))]
  for (col in num) set(m, which(is.na(m[[col]])), col, 0)
  m <- as.data.frame(m)
  m$display_label <- paste0(m$gene, " | ", m$transcript_id)
  AVG_CACHE[[id]] <- m
  m
}
diversity_range_df   <- as.data.frame(fread(file.path(DATA_ROOT, REF_ID, "diversity_range.csv")))
if ("most_diverse_group" %in% names(diversity_range_df)) {
  diversity_range_df$most_diverse_tissue  <- diversity_range_df$most_diverse_group
  diversity_range_df$least_diverse_tissue <- diversity_range_df$least_diverse_group
}
dominant_isoforms_df <- as.data.frame(fread(file.path(DATA_ROOT, REF_ID, "dominant_isoforms.csv")))
if ("group_label" %in% names(dominant_isoforms_df)) {
  dominant_isoforms_df$tissue <- dominant_isoforms_df$group_label
}
gdsc_response_df     <- as.data.frame(fread("shiny_data/gdsc2_matched.csv"))
gdsc1_response_df <- if (file.exists("shiny_data/gdsc1_matched.csv")) {
  as.data.frame(fread("shiny_data/gdsc1_matched.csv"))
} else {
  NULL
}
transcript_coords_df <- as.data.frame(fread("shiny_data/transcript_coords.csv"))
if ("gene" %in% names(transcript_coords_df) && !"gene_name" %in% names(transcript_coords_df))
  names(transcript_coords_df)[names(transcript_coords_df) == "gene"] <- "gene_name"

gene_tpm_path <- file.path(DATA_ROOT, REF_ID, "salmon_gene_tpm.csv.gz")
HAS_GENE_TPM <- file.exists(gene_tpm_path)
gene_tpm_df <- if (HAS_GENE_TPM) {
  cat("Loading gene-level TPM matrix...\n")
  g <- as.data.frame(fread(gene_tpm_path))
  id2label <- setNames(sample_meta_df$sample_label, sample_meta_df$sample_id)
  hit <- names(g) %in% names(id2label)
  names(g)[hit] <- id2label[names(g)[hit]]
  g
} else {
  NULL
}

if (!is.null(protein_annot_df) && !is.null(tx_annot_df) &&
    !any(grepl("\\.", head(protein_annot_df$transcript_id, 200)))) {
  ver_map <- setNames(tx_annot_df$transcript_id, sub("\\..*", "", tx_annot_df$transcript_id))
  ver_map <- ver_map[!duplicated(names(ver_map))]
  protein_annot_df$transcript_id <- unname(ver_map[protein_annot_df$transcript_id])
  protein_annot_df <- protein_annot_df[!is.na(protein_annot_df$transcript_id), ]
}
cat("Versioned IDs: annotations", any(grepl("\\.", head(tx_annot_df$transcript_id))),
    "| exons", any(grepl("\\.", head(tx_exons_df$transcript_id))),
    "| coords", any(grepl("\\.", head(transcript_coords_df$transcript_id))), "\n")

QUERY_CON <- new.env()
DESEQ_CON <- new.env()
AVG_CACHE <- new.env()

query_con <- function(id) {
  if (is.null(id) || id == REF_ID) return(NULL)
  con <- QUERY_CON[[id]]
  if (is.null(con) || !dbIsValid(con)) {
    row <- DATASETS[DATASETS$dataset_id == id, ]
    if (nrow(row) != 1) return(NULL)
    con <- dbConnect(RSQLite::SQLite(), row$db_path)
    QUERY_CON[[id]] <- con
  }
  con
}

deseq_con <- function(id) {
  if (id == REF_ID) return(deseq_db)
  con <- DESEQ_CON[[id]]
  if (is.null(con) || !dbIsValid(con)) {
    row <- DATASETS[DATASETS$dataset_id == id, ]
    if (nrow(row) != 1 || !row$has_deseq) return(NULL)
    con <- dbConnect(RSQLite::SQLite(), row$deseq_path)
    DESEQ_CON[[id]] <- con
  }
  con
}

tissue_rename <- c(
  "LARGE" = "LARGE INTESTINE",
  "SMALL" = "SMALL INTESTINE",
  "SOFT" = "SOFT TISSUE",
  "HAEMATOPOIETIC_AND_LYMPHOID_TISSUE" = "HAEMATOPOIETIC"
)
tissue_unrename <- setNames(names(tissue_rename), tissue_rename)
rename_tissue <- function(x) {
  r <- tissue_rename[x]
  ifelse(is.na(r), x, r)
}
unrename_tissue <- function(x) {
  r <- tissue_unrename[x]
  ifelse(is.na(r), x, r)
}

sample_meta_df$Tissue <- rename_tissue(sample_meta_df$Tissue)
tissue_avg_df <- tissue_avg_df %>% rename_with(~ rename_tissue(.x), any_of(names(tissue_rename)))
diversity_range_df$most_diverse_tissue <- rename_tissue(diversity_range_df$most_diverse_tissue)
diversity_range_df$least_diverse_tissue <- rename_tissue(diversity_range_df$least_diverse_tissue)
dominant_isoforms_df$tissue <- rename_tissue(dominant_isoforms_df$tissue)

gene_list      <- sort(unique(tissue_avg_df$gene))
tissue_list    <- sort(unique(sample_meta_df$Tissue))
sample_cols    <- sample_meta_df$sample_label
cell_line_list <- sort(unique(sample_meta_df$CellLine))
drug_list      <- sort(unique(gdsc_response_df$DRUG_NAME))
drug_list_1    <- if (!is.null(gdsc1_response_df)) {
  sort(unique(gdsc1_response_df$DRUG_NAME))
} else {
  character(0)
}
TISSUE_COLORS  <- tissue_colors(tissue_list)

GROUP_LABELS <- c(
  "Adjacent" = "Adjacent",
  "Tumor" = "Tumor",
  "BmRC" = "Primary (bone-met)",
  "NBmRC" = "Primary (non-met)",
  "Bone_met" = "Bone met",
  "LN" = "LN met",
  "PRE" = "Pre-G207",
  "POST" = "Post-G207"
)
DATASET_PREFIX <- c("hcc" = "HCC", "pca" = "PCa", "gbm" = "GBM")

pretty_group <- function(g, ds) {
  ds <- rep_len(ds, length(g))
  lbl <- ifelse(g %in% names(GROUP_LABELS), GROUP_LABELS[g], g)
  pre <- DATASET_PREFIX[ds]
  ifelse(is.na(pre), lbl, paste(pre, lbl))
}

QUERY_GROUPS <- do.call(rbind, lapply(QUERY_IDS, function(id) {
  con <- dbConnect(RSQLite::SQLite(), DATASETS$db_path[DATASETS$dataset_id == id])
  on.exit(dbDisconnect(con))
  g <- dbGetQuery(con, "SELECT group_label, COUNT(*) AS n FROM samples GROUP BY group_label")
  data.frame(dataset = id, raw = g$group_label, n = g$n,
             pretty = pretty_group(g$group_label, id), stringsAsFactors = FALSE)
}))

GROUP_FULL <- c(
  "Adjacent" = "Adjacent normal",
  "Tumor" = "Tumor",
  "BmRC" = "Primary, bone-metastatic",
  "NBmRC" = "Primary, non-metastatic",
  "Bone_met" = "Bone metastasis",
  "LN" = "Lymph node metastasis",
  "PRE" = "Pre-G207",
  "POST" = "Post-G207"
)
DATASET_FULL <- c("hcc" = "Liver", "pca" = "Prostate", "gbm" = "Glioblastoma")
if (!is.null(QUERY_GROUPS) && nrow(QUERY_GROUPS) > 0) {
  full_lbl <- ifelse(QUERY_GROUPS$raw %in% names(GROUP_FULL), GROUP_FULL[QUERY_GROUPS$raw], QUERY_GROUPS$raw)
  QUERY_GROUPS$full <- paste0(DATASET_FULL[QUERY_GROUPS$dataset], ": ", full_lbl)
}

legend_label <- function(x) {
  if (is.null(QUERY_GROUPS) || nrow(QUERY_GROUPS) == 0) return(x)
  r <- QUERY_GROUPS$full[match(x, QUERY_GROUPS$pretty)]
  ifelse(is.na(r), x, r)
}

raw_group <- function(x) {
  r <- QUERY_GROUPS$raw[match(x, QUERY_GROUPS$pretty)]
  ifelse(is.na(r), x, r)
}

if (!is.null(QUERY_GROUPS) && nrow(QUERY_GROUPS) > 0) {
  q_pal <- c("#1b9e77", "#d95f02", "#7570b3", "#e7298a", "#66a61e", "#e6ab02",
             "#a6761d", "#1f78b4", "#b2182b", "#4d4d4d")
  q_pal <- if (nrow(QUERY_GROUPS) <= length(q_pal)) q_pal[seq_len(nrow(QUERY_GROUPS))]
           else colorRampPalette(q_pal)(nrow(QUERY_GROUPS))
  TISSUE_COLORS <- c(TISSUE_COLORS, setNames(q_pal, QUERY_GROUPS$pretty))
}

expr_ok <- function(a, b, thr, rule) {
  if (identical(rule, "either")) pmax(a, b) >= thr else pmin(a, b) >= thr
}

CCLE_SAMPLES <- sample_meta_df$subject_id
GBM_SAMPLES  <- character(0)

biotype_list <- if (HAS_TX_ANNOT) {
  sort(unique(tx_annot_df$biotype))
} else {
  character(0)
}
tissue_counts <- as.integer(table(sample_meta_df$Tissue)[tissue_list])
tissue_order <- order(tissue_counts, decreasing = TRUE)
tissue_list <- tissue_list[tissue_order]
tissue_counts <- tissue_counts[tissue_order]
tissue_choices <- setNames(tissue_list, paste0(tissue_list, " (n=", tissue_counts, ")"))
CCLE_TISSUES <- tissue_list
GBM_TISSUES  <- character(0)

default_gene <- if ("FGFR1" %in% gene_list) {
  "FGFR1"
} else if ("EGFR" %in% gene_list) {
  "EGFR"
} else if ("KRAS" %in% gene_list) {
  "KRAS"
} else {
  gene_list[which(!grepl("^[0-9]", gene_list))[1]]
}

refseq_map_path <- "shiny_data/enst_to_refseq.csv"
HAS_REFSEQ <- file.exists(refseq_map_path)
refseq_map_df <- if (HAS_REFSEQ) {
  cat("Loading RefSeq mapping...\n")
  as.data.frame(fread(refseq_map_path))
} else {
  NULL
}

cat("Startup complete.\n")

isoscope_logo_light <- function(size = 40) {
  tags$svg(
    xmlns = "http://www.w3.org/2000/svg",
    width = paste0(size, "px"),
    height = paste0(size, "px"),
    viewBox = "0 0 144 144",
    style = "display: inline-block; vertical-align: middle;",
    tags$circle(
      cx = "60",
      cy = "60",
      r = "55",
      fill = "none",
      stroke = "white",
      `stroke-width` = "3"
    ),
    tags$circle(
      cx = "60",
      cy = "60",
      r = "46",
      fill = "none",
      stroke = "#41b6c4",
      `stroke-width` = "1.5",
      `stroke-dasharray` = "4,3"
    ),
    tags$rect(
      x = "22",
      y = "38",
      width = "76",
      height = "7",
      rx = "3.5",
      fill = "white",
      opacity = "0.95"
    ),
    tags$rect(
      x = "30",
      y = "50",
      width = "60",
      height = "7",
      rx = "3.5",
      fill = "#41b6c4",
      opacity = "0.9"
    ),
    tags$rect(
      x = "22",
      y = "62",
      width = "45",
      height = "7",
      rx = "3.5",
      fill = "#c51b8a",
      opacity = "0.85"
    ),
    tags$rect(
      x = "35",
      y = "74",
      width = "55",
      height = "7",
      rx = "3.5",
      fill = "white",
      opacity = "0.8"
    ),
    tags$line(
      x1 = "99",
      y1 = "99",
      x2 = "130",
      y2 = "130",
      stroke = "white",
      `stroke-width` = "5",
      `stroke-linecap` = "round"
    )
  )
}

isoscope_logo <- function(size = 120) {
  tags$svg(
    xmlns = "http://www.w3.org/2000/svg",
    width = paste0(size * 1.2, "px"),
    height = paste0(size * 1.2, "px"),
    viewBox = "0 0 144 144",
    style = "display: inline-block; vertical-align: middle;",
    tags$circle(
      cx = "60",
      cy = "60",
      r = "55",
      fill = "none",
      stroke = "#225ea8",
      `stroke-width` = "3"
    ),
    tags$circle(
      cx = "60",
      cy = "60",
      r = "46",
      fill = "none",
      stroke = "#41b6c4",
      `stroke-width` = "1.5",
      `stroke-dasharray` = "4,3"
    ),
    tags$rect(
      x = "22",
      y = "38",
      width = "76",
      height = "7",
      rx = "3.5",
      fill = "#225ea8",
      opacity = "0.9"
    ),
    tags$rect(
      x = "30",
      y = "50",
      width = "60",
      height = "7",
      rx = "3.5",
      fill = "#41b6c4",
      opacity = "0.85"
    ),
    tags$rect(
      x = "22",
      y = "62",
      width = "45",
      height = "7",
      rx = "3.5",
      fill = "#c51b8a",
      opacity = "0.8"
    ),
    tags$rect(
      x = "35",
      y = "74",
      width = "55",
      height = "7",
      rx = "3.5",
      fill = "#225ea8",
      opacity = "0.7"
    ),
    tags$line(
      x1 = "99",
      y1 = "99",
      x2 = "130",
      y2 = "130",
      stroke = "#225ea8",
      `stroke-width` = "5",
      `stroke-linecap` = "round"
    )
  )
}

# --- UI ---

tab_gene_viewer <- tabPanel("Gene Viewer", sidebarLayout(
  sidebarPanel(
    width = 3,
    selectizeInput(
      "gene",
      "Search Gene:",
      choices = NULL,
      options = list(placeholder = "Type a gene name...")
    ),
    selectInput(
      "tissue_filter",
      "View Mode:",
      choices = c("All tissues (averaged)" = "All", tissue_choices),
      selected = "All"
    ),
    selectizeInput(
      "gene_tissue_multi",
      "Show Tissues:",
      choices = tissue_choices,
      selected = NULL,
      multiple = TRUE,
      options = list(placeholder = "All tissues (leave empty for all)")
    ),
    helpText(
      tags$b("View Mode:"),
      "'All tissues' shows one bar per tissue or patient group (averaged across samples).",
      "Selecting a specific tissue or group shows its individual samples.",
      br(),
      tags$b("Show Tissues:"),
      "In 'All tissues' mode, pick a subset of tissues or groups to display. Leave empty for all."
    ),
    sliderInput(
      "top_n",
      "Show top N isoforms:",
      min = 3,
      max = 20,
      value = 5,
      step = 1
    ),
    radioButtons(
      "plot_type",
      "Plot Type:",
      choices = c(
        "Stacked Bar (proportions)" = "stacked",
        "Grouped Bar (absolute TPM)" = "grouped",
        "Heatmap" = "heatmap"
      ),
      selected = "stacked"
    ),
    checkboxInput("log_scale", "Log2(TPM + 1) scale", value = FALSE),
    if (HAS_TX_ANNOT)
      selectInput(
        "gene_biotype_filter",
        "Filter by Biotype:",
        choices = c("All", biotype_list),
        selected = "All"
      ),
    actionButton("run_gene", "Update Plot", class = "btn-primary"),
    hr(),
    downloadButton("dl_plot_gene", "Download Plot (PDF)"),
    downloadButton("dl_table_gene", "Download Table (CSV)"),
    hr(),
    bookmarkButton(label = "Bookmark this view"),
    hr(),
    helpText(
      "'All tissues' averages expression per tissue or group. Select one to see its individual samples."
    )
  ),
  mainPanel(
    width = 9,
    plotlyOutput("isoform_plot", height = "550px"),
    hr(),
    h4("Isoform Expression Table"),
    DTOutput("isoform_table")
  )
))

tab_switching <- tabPanel("Isoform Switching", sidebarLayout(
  sidebarPanel(
    width = 3,
    radioButtons(
      "switch_mode",
      "Comparison Mode:",
      choices = c(
        "Group vs Group" = "tissue",
        "Sample vs Sample" = "cellline"
      ),
      selected = "tissue"
    ),
    conditionalPanel(
      condition = "input.switch_mode == 'tissue'",
      selectInput(
        "switch_tissue1",
        "Group 1:",
        choices = tissue_choices,
        selected = if ("LUNG" %in% tissue_list)
          "LUNG"
        else
          tissue_list[1]
      ),
      selectInput(
        "switch_tissue2",
        "Group 2:",
        choices = tissue_choices,
        selected = if ("BREAST" %in% tissue_list)
          "BREAST"
        else
          tissue_list[min(2, length(tissue_list))]
      )
    ),
    conditionalPanel(
      condition = "input.switch_mode == 'cellline'",
      selectizeInput(
        "switch_cl1",
        "Sample 1:",
        choices = NULL,
        options = list(placeholder = "Type a cell line or patient...")
      ),
      selectizeInput(
        "switch_cl2",
        "Sample 2:",
        choices = NULL,
        options = list(placeholder = "Type a cell line or patient...")
      ),
      uiOutput("switch_cl_tissues")
    ),
    radioButtons(
      "switch_expr_rule",
      "TPM threshold must be met by:",
      choices = c("Both groups" = "both", "Either group" = "either"),
      selected = "both",
      inline = TRUE
    ),
    sliderInput(
      "min_expr",
      "Min avg TPM:",
      min = 0,
      max = 50,
      value = 1,
      step = 0.5
    ),
    sliderInput(
      "min_switch_score",
      "Min switch score:",
      min = 0,
      max = 1,
      value = 0.2,
      step = 0.05
    ),
    conditionalPanel(
      condition = "input.switch_mode == 'tissue'",
      selectInput(
        "switch_confidence",
        "Min confidence:",
        choices = c(
          "All" = "all",
          "High only" = "high",
          "High + Medium" = "high_med"
        ),
        selected = "high_med"
      )
    ),
    if (HAS_TX_ANNOT)
      selectInput(
        "switch_biotype_filter",
        "Filter by Biotype:",
        choices = c("All", biotype_list),
        selected = "All"
      ),
    actionButton("run_switch", "Detect Switches", class = "btn-primary"),
    hr(),
    downloadButton("dl_switch_table", "Download Results (CSV)"),
    downloadButton("dl_switch_plot_pdf", "Download Plot (PDF)"),
    hr(),
    conditionalPanel(
      condition = "input.switch_mode == 'tissue'",
      helpText(
        tags$b("What is isoform switching?"),
        br(),
        "A gene can produce multiple transcript isoforms. Isoform switching occurs when",
        "different tissues, patient groups or samples use a different dominant isoform for the same gene.",
        br(),
        br(),
        tags$b("How it works:"),
        "For each gene, the tool identifies the most highly expressed",
        "isoform in each group. If the dominant isoform differs between the two groups, that is a switch.",
        br(),
        br(),
        tags$b("Reading the table:"),
        br(),
        tags$b("Prop"),
        "= proportion of total gene TPM from the dominant isoform (0-1).",
        br(),
        tags$b("n"),
        "= number of samples in that group.",
        br(),
        tags$b("Switch score"),
        "= |p1 - p2|, where p1 and p2 are the proportions of each group's own dominant isoform.",
        "Because the two dominant isoforms are different transcripts, the score reflects how concentrated",
        "expression is in each group rather than the change in a single isoform; a balanced shift",
        "(for example 60/40 to 40/60) scores close to 0.",
        br(),
        tags$b("Confidence"),
        "= based on bootstrap resampling (1000 iterations).",
        br(),
        "High = in each group, the lower 95% CI bound of that group's dominant isoform exceeds that group's proportion of the other group's dominant isoform (both groups n >= 3).",
        br(),
        "Medium = the dominant isoform differs, but this holds in only one group or neither.",
        br(),
        "Low = fewer than 3 samples in a group.",
        br(),
        br(),
        "Click a row to see the isoform proportion breakdown below."
      )
    ),
    conditionalPanel(
      condition = "input.switch_mode == 'cellline'",
      helpText(
        tags$b("Sample comparison:"),
        "Compares dominant isoforms between two individual samples (cell lines or patient samples).",
        br(),
        br(),
        tags$b("Prop"),
        "= proportion of total gene TPM from the dominant isoform.",
        br(),
        tags$b("Switch score"),
        "= |p1 - p2|, the difference between each sample's own dominant-isoform proportion.",
        br(),
        "No confidence intervals (single sample per group).",
        br(),
        "Click a row to see the isoform proportion breakdown below."
      )
    )
  ),
  mainPanel(
    width = 9,
    h4("Isoform Switching Events"),
    DTOutput("switch_table"),
    hr(),
    h4("Selected Gene Detail"),
    plotlyOutput("switch_plot", height = "400px"),
    hr(),
    h4("Exon Structure Comparison"),
    helpText(
      "The two dominant isoforms shown side by side. Blue tiles = exon present in both isoforms (shared). ",
      "Pink tiles = exon present in only one of the two (the structural difference driving the switch). ",
      "Each tile is labeled with that exon's internal number; x-axis shows the gene-wide exon slot 5' to 3'."
    ),
    plotlyOutput("switch_exon_plot", height = "360px"),
    hr(),
    h4("Single Transcript Detail"),
    helpText(
      "Pick either of the two dominant isoforms from the dropdown (or click a row in the grid above)."
    ),
    selectizeInput(
      "switch_zoom_tx",
      "Transcript:",
      choices = NULL,
      options = list(placeholder = "Pick a transcript or click a row in the grid above"),
      width = "100%"
    ),
    plotlyOutput("switch_zoom_plot", height = "220px"),
    downloadButton("dl_switch_exon_zoom", "Download Transcript Detail (PDF)")
  )
))

tab_qc <- tabPanel("QC & Distribution", sidebarLayout(
  sidebarPanel(
    width = 3,
    selectizeInput(
      "qc_tissue",
      "Filter by Tissue(s):",
      choices = tissue_choices,
      selected = NULL,
      multiple = TRUE,
      options = list(placeholder = "All tissues (click to filter)")
    ),
    radioButtons(
      "qc_metric",
      "Metric:",
      choices = c(
        "Unique Isoforms Detected" = "n_isoforms",
        "Unique Genes Detected" = "n_genes",
        "Isoforms per Gene (avg)" = "iso_per_gene"
      ),
      selected = "n_isoforms"
    ),
    actionButton("run_qc", "Run", class = "btn-primary"),
    hr(),
    downloadButton("dl_qc_table", "Download QC Table (CSV)"),
    downloadButton("dl_qc_plot", "Download Plot (PDF)"),
    hr(),
    helpText(
      "Quality control overview of the isoform quantification data.",
      "Each point is one sample. Boxplots show the distribution within each tissue or patient group."
    )
  ),
  mainPanel(
    width = 9,
    plotlyOutput("qc_plot", height = "500px"),
    hr(),
    h4("Per-Sample Summary"),
    DTOutput("qc_table")
  )
))

tab_diversity <- tabPanel("Isoform Diversity", sidebarLayout(
  sidebarPanel(
    width = 3,
    selectizeInput(
      "div_gene",
      "Search Gene:",
      choices = NULL,
      options = list(placeholder = "Type a gene name...")
    ),
    radioButtons(
      "div_mode",
      "View Mode:",
      choices = c("By Tissue (averaged)" = "tissue", "By Sample" = "cellline"),
      selected = "tissue"
    ),
    conditionalPanel(
      condition = "input.div_mode == 'cellline'",
      selectizeInput(
        "div_tissue_filter",
        "Filter by Tissue(s):",
        choices = tissue_choices,
        selected = NULL,
        multiple = TRUE,
        options = list(placeholder = "All tissues (click to filter)")
      )
    ),
    radioButtons(
      "div_plot_type",
      "Plot Type:",
      choices = c("Bar Chart" = "bar", "Bubble Plot" = "bubble"),
      selected = "bubble"
    ),
    sliderInput(
      "div_min_tpm",
      "Min total TPM (filter noise):",
      min = 0,
      max = 10,
      value = 1,
      step = 0.5
    ),
    actionButton("run_diversity", "Run", class = "btn-primary"),
    hr(),
    helpText(
      "Normalized Shannon entropy measures how evenly isoforms are expressed.",
      "1 = all isoforms equally expressed. 0 = one isoform dominates.",
      "Normalized by isoform count so genes with different numbers of isoforms are comparable.",
      br(),
      br(),
      tags$b("Interpreting the scale:"),
      br(),
      "0.0 = one isoform dominates entirely",
      br(),
      "0.5 = moderate diversity",
      br(),
      "1.0 = perfect evenness"
    ),
    hr(),
    downloadButton("dl_diversity", "Download Per-Gene Diversity (CSV)"),
    downloadButton("dl_diversity_ranking", "Download Global Ranking (CSV)"),
    downloadButton("dl_diversity_plot", "Download Plot (PDF)")
  ),
  mainPanel(
    width = 9,
    plotlyOutput("diversity_plot", height = "450px"),
    hr(),
    h4("Isoform Exon Structure"),
    helpText(
      "Each row is an isoform of this gene; columns are unique exons across all isoforms (ordered 5' to 3'). ",
      "A filled tile means the isoform contains that exon. ",
      "Click a row in the per-tissue table below to highlight that tissue's dominant isoform in pink."
    ),
    plotlyOutput("diversity_exon_plot", height = "500px"),
    downloadButton("dl_diversity_exon_grid", "Download Exon Grid (PDF)"),
    hr(),
    h4("Single Transcript Detail"),
    helpText(
      "Click any row in the composite grid above to zoom into that transcript, or pick one from the dropdown. ",
      "The detail view shows true genomic proportions: exons as solid blocks sized by their real length in base pairs, ",
      "and introns as thin connecting lines. Each exon is labeled with its internal number for that isoform."
    ),
    selectizeInput(
      "div_zoom_tx",
      "Transcript:",
      choices = NULL,
      options = list(placeholder = "Pick a transcript or click a row in the grid above"),
      width = "100%"
    ),
    plotlyOutput("diversity_zoom_plot", height = "220px"),
    downloadButton("dl_diversity_exon_zoom", "Download Transcript Detail (PDF)"),
    hr(),
    h4("Per-Tissue Normalized Entropy"),
    DTOutput("diversity_gene_table"),
    hr(),
    h4("Genes Ranked by Diversity Range"),
    helpText("All genes ranked by normalized entropy range across tissues."),
    DTOutput("diversity_range_table")
  )
))

tab_tissue_specific <- tabPanel("Tissue-Specific", sidebarLayout(
  sidebarPanel(
    width = 3,
    selectizeInput(
      "spec_gene",
      "Search a specific gene (optional):",
      choices = NULL,
      options = list(placeholder = "Leave blank to scan all genes")
    ),
    helpText(
      style = "color: #666; font-size: 12px; margin-top: -8px; margin-bottom: 10px;",
      "Leave empty to find all tissue-specific isoforms in the selected tissue. ",
      "Type a gene name to check whether that one gene is tissue-specific."
    ),
    selectInput(
      "spec_tissue",
      "Tissue:",
      choices = tissue_choices,
      selected = tissue_list[1]
    ),
    sliderInput(
      "spec_threshold",
      "Tissue fraction threshold:",
      min = 0.3,
      max = 1.0,
      value = 0.5,
      step = 0.05
    ),
    sliderInput(
      "spec_min_tpm",
      "Min avg TPM in tissue:",
      min = 0,
      max = 20,
      value = 0.5,
      step = 0.5
    ),
    if (HAS_TX_ANNOT)
      selectInput(
        "spec_biotype_filter",
        "Filter by Biotype:",
        choices = c("All", biotype_list),
        selected = "All"
      ),
    actionButton("run_specific", "Find Tissue-Specific Isoforms", class = "btn-primary"),
    hr(),
    radioButtons(
      "spec_detail_mode",
      "Detail Plot:",
      choices = c(
        "Tissue averages" = "tissue",
        "Per-cell-line distribution" = "cellline"
      ),
      selected = "tissue"
    ),
    downloadButton("dl_specific", "Download Results (CSV)"),
    downloadButton("dl_specific_plot", "Download Plot (PDF)"),
    hr(),
    helpText(
      tags$b("What is tissue specificity?"),
      br(),
      "An isoform is called tissue-specific if its mean TPM in the selected tissue accounts for ",
      "a large fraction of its summed mean TPM across all tissues in the active dataset.",
      br(),
      br(),
      tags$b("Avg TPM:"),
      "Mean TPM of this isoform across all samples in the selected tissue or group.",
      br(),
      tags$b("Tissue Fraction:"),
      "target group mean / sum of group means (0 to 1). With 26 or more groups, ",
      "a fraction of 0.5 already implies strong enrichment.",
      br(),
      br(),
      tags$b("Caveat:"),
      "CCLE samples are cancer cell lines, not normal tissue. Differentiation markers ",
      "(e.g. A4GNT in stomach, normal hepatocyte markers in liver) are often lost in ",
      "transformed cells and may show very low TPM even in their tissue of origin.",
      br(),
      br(),
      tags$b("If a single-gene query returns no hits,"),
      "the table below will show that gene's per-isoform TPM and fraction so you can ",
      "see why it failed the cutoffs.",
      br(),
      br(),
      "Click a row to see the expression detail plot below."
    )
  ),
  mainPanel(
    width = 9,
    h4("Tissue-Specific Isoforms"),
    DTOutput("specific_table"),
    hr(),
    h4("Selected Isoform Detail"),
    plotlyOutput("specific_plot", height = "400px"),
    hr(),
    h4("Exon Structure"),
    helpText(
      "The selected tissue-specific isoform (highlighted in pink) shown alongside all other isoforms of the gene. ",
      "Use this to see which exon(s) make this isoform structurally distinct from the rest of the gene's transcripts."
    ),
    plotlyOutput("specific_exon_plot", height = "450px"),
    hr(),
    h4("Single Transcript Detail"),
    helpText(
      "Click any row in the grid above to zoom into that transcript, or pick one from the dropdown."
    ),
    selectizeInput(
      "spec_zoom_tx",
      "Transcript:",
      choices = NULL,
      options = list(placeholder = "Pick a transcript or click a row in the grid above"),
      width = "100%"
    ),
    plotlyOutput("spec_zoom_plot", height = "220px"),
    downloadButton("dl_specific_exon_zoom", "Download Transcript Detail (PDF)"),
    hr(),
    h4(textOutput("specific_detail_header")),
    DTOutput("specific_detail_table")
  )
))

tab_clustering <- tabPanel("Clustering", sidebarLayout(
  sidebarPanel(
    width = 3,
    selectizeInput(
      "clust_gene",
      "Gene (or All):",
      choices = NULL,
      options = list(placeholder = "Type gene or leave blank for all")
    ),
    radioButtons(
      "clust_reduction",
      "Method:",
      choices = c("PCA" = "pca", "UMAP" = "umap"),
      selected = "pca"
    ),
    if (HAS_GENE_TPM)
      radioButtons(
        "clust_level",
        "Expression Level:",
        choices = c("Isoform (Salmon TPM)" = "isoform", "Gene (tximport TPM)" = "gene"),
        selected = "isoform"
      ),
    checkboxInput(
      "clust_tissue_avg",
      "Cluster by group averages instead of individual samples",
      value = FALSE
    ),
    checkboxInput("clust_3d", "3D plot (use 3 components)", value = FALSE),
    selectizeInput(
      "clust_tissue_filter",
      "Filter by Tissue(s):",
      choices = tissue_choices,
      selected = NULL,
      multiple = TRUE,
      options = list(placeholder = "All tissues (click to filter)")
    ),
    helpText(
      "Leave empty for all tissues. Select specific tissues to focus the clustering.",
      "At least 3 samples are needed for PCA and 5 for UMAP."
    ),
    if (HAS_TX_ANNOT)
      selectInput(
        "clust_biotype_filter",
        "Filter by Biotype:",
        choices = c("All", biotype_list),
        selected = "All"
      ),
    actionButton("run_cluster", "Run Clustering", class = "btn-primary"),
    hr(),
    downloadButton("dl_cluster_plot", "Download Plot (PDF)"),
    conditionalPanel(
      condition = "input.clust_reduction == 'pca'",
      selectInput(
        "clust_method",
        "Distance Metric:",
        choices = c("Correlation" = "cor", "Euclidean" = "euclidean"),
        selected = "cor"
      )
    ),
    conditionalPanel(
      condition = "input.clust_reduction == 'umap'",
      sliderInput(
        "umap_neighbors",
        "N neighbors:",
        min = 5,
        max = 100,
        value = 15,
        step = 5
      ),
      sliderInput(
        "umap_min_dist",
        "Min distance:",
        min = 0.01,
        max = 1.0,
        value = 0.1,
        step = 0.05
      )
    ),
    sliderInput(
      "clust_top_var",
      "Use top N variable features:",
      min = 50,
      max = 5000,
      value = 500,
      step = 50
    ),
    sliderInput(
      "clust_dot_size",
      "Dot size:",
      min = 1,
      max = 10,
      value = 3,
      step = 0.5
    ),
    hr(),
    helpText(
      tags$b("What this does:"),
      "Samples are clustered by expression patterns (log2(TPM+1) transformed).",
      "When a query dataset is selected, its samples are projected onto the reference embedding (Isoform level).",
      br(),
      br(),
      tags$b("Isoform (Salmon TPM):"),
      "Captures splicing differences between samples.",
      br(),
      tags$b("Gene (tximport TPM):"),
      "Gene-level expression summarised from Salmon with tximport. Reference cell lines only.",
      br(),
      br(),
      tags$b("PCA:"),
      "Linear projection.",
      br(),
      tags$b("UMAP:"),
      "Non-linear embedding.",
      br(),
      br(),
      tags$b("Group averages:"),
      "One point per tissue or group instead of one per sample."
    )
  ),
  mainPanel(
    width = 9,
    h4("Clustering by Expression"),
    helpText(
      style = "color: #888; font-size: 12px;",
      "For 3D plots, use the camera icon in the plot toolbar to save as PNG."
    ),
    plotlyOutput("cluster_plot", height = "600px"),
    hr(),
    h4("Top Features Driving Clustering"),
    helpText(
      "For PCA: features ranked by combined loading magnitude on the first principal components ",
      "(higher = more influence on sample positioning). For UMAP: features ranked by variance across samples ",
      "(UMAP does not produce loadings, so variance is used as a proxy for contribution)."
    ),
    DTOutput("cluster_top_features_table"),
    downloadButton("dl_cluster_top_features", "Download Table (CSV)"),
    br(),
    br(),
    hr(),
    conditionalPanel(
      condition = "input.clust_reduction == 'pca' && !input.clust_3d",
      h4("PC Variance Explained"),
      plotlyOutput("scree_plot", height = "250px")
    ),
    hr(),
    h4("Top Features Driving Clustering"),
    DTOutput("cluster_features_table")
  )
))

tab_corr_heatmap <- tabPanel("Tissue Correlation",
                             sidebarLayout(
                               sidebarPanel(
                                 width = 3,
                                 selectizeInput(
                                   "corr_gene",
                                   "Gene (or All):",
                                   choices = NULL,
                                   options = list(placeholder = "Type gene or leave blank for all")
                                 ),
                                 radioButtons(
                                   "corr_mode",
                                   "Compare:",
                                   choices = c(
                                     "Tissue vs Tissue" = "tissue",
                                     "Sample vs Sample" = "cellline"
                                   ),
                                   selected = "tissue"
                                 ),
                                 selectizeInput(
                                   "corr_tissue_filter",
                                   "Filter by Tissue(s):",
                                   choices = tissue_choices,
                                   selected = NULL,
                                   multiple = TRUE,
                                   options = list(placeholder = "All tissues (click to filter)")
                                 ),
                                 helpText("Leave empty for all tissues. Select specific tissues to filter."),
                                 radioButtons(
                                   "corr_method",
                                   "Correlation Method:",
                                   choices = c("Pearson" = "pearson", "Spearman" = "spearman"),
                                   selected = "pearson"
                                 ),
                                 sliderInput(
                                   "corr_top_var",
                                   "Use top N variable isoforms:",
                                   min = 50,
                                   max = 5000,
                                   value = 500,
                                   step = 50
                                 ),
                                 if (HAS_TX_ANNOT)
                                   selectInput(
                                     "corr_biotype_filter",
                                     "Filter by Biotype:",
                                     choices = c("All", biotype_list),
                                     selected = "All"
                                   ),
                                 actionButton("run_corr", "Compute Correlation", class = "btn-primary"),
                                 hr(),
                                 downloadButton("dl_corr_plot", "Download Plot (PDF)"),
                                 downloadButton("dl_corr_table", "Download Correlation Table (CSV)"),
                                 hr(),
                                 helpText("Uses log2(TPM+1) transformed expression.")
                               ),
                               mainPanel(
                                 width = 9,
                                 h4("Correlation Heatmap"),
                                 plotlyOutput("corr_heatmap", height = "600px"),
                                 hr(),
                                 h4("Correlation Table"),
                                 DTOutput("corr_table")
                               )
                             ))

de_method_choices <- c("Wilcoxon Transcript (TPM)" = "wilcoxon")
if (any(DATASETS$has_deseq))
  de_method_choices <- c(
    de_method_choices,
    "DESeq2 Gene-Level" = "deseq2_gene",
    "DESeq2 Transcript-Level (Isoforms)" = "deseq2_tx"
  )

tab_de <- tabPanel("Differential Expression",
                   sidebarLayout(
                     sidebarPanel(
                       width = 3,
                       radioButtons(
                         "de_method",
                         "Method:",
                         choices = de_method_choices,
                         selected = "wilcoxon"
                       ),
                       radioButtons(
                         "de_mode",
                         "Comparison Mode:",
                         choices = c(
                           "Group vs Group" = "tissue",
                           "Sample vs Sample" = "cellline"
                         ),
                         selected = "tissue"
                       ),
                       conditionalPanel(
                         condition = "input.de_mode == 'tissue'",
                         selectInput(
                           "de_tissue1",
                           "Group 1:",
                           choices = tissue_choices,
                           selected = if ("LUNG" %in% tissue_list)
                             "LUNG"
                           else
                             tissue_list[1]
                         ),
                         selectInput(
                           "de_tissue2",
                           "Group 2:",
                           choices = tissue_choices,
                           selected = if ("BREAST" %in% tissue_list)
                             "BREAST"
                           else
                             tissue_list[min(2, length(tissue_list))]
                         ),
                         uiOutput("de_availability")
                       ),
                       conditionalPanel(
                         condition = "input.de_mode == 'cellline'",
                         selectizeInput(
                           "de_cl1",
                           "Sample 1:",
                           choices = NULL,
                           options = list(placeholder = "Type a cell line or patient...")
                         ),
                         selectizeInput(
                           "de_cl2",
                           "Sample 2:",
                           choices = NULL,
                           options = list(placeholder = "Type a cell line or patient...")
                         ),
                         uiOutput("de_cl_tissues")
                       ),
                       conditionalPanel(
                         condition = "input.de_method == 'wilcoxon'",
                         radioButtons(
                           "de_expr_rule",
                           "TPM threshold must be met by:",
                           choices = c("Both groups" = "both", "Either group" = "either"),
                           selected = "both",
                           inline = TRUE
                         ),
                         sliderInput(
                           "de_min_tpm",
                           "Min mean TPM:",
                           min = 0,
                           max = 20,
                           value = 1,
                           step = 0.5
                         )
                       ),
                       if (HAS_TX_ANNOT)
                         selectInput(
                           "de_biotype_filter",
                           "Filter by Biotype:",
                           choices = c("All", biotype_list),
                           selected = "All"
                         ),
                       actionButton("run_de", "Run Analysis", class = "btn-primary"),
                       hr(),
                       sliderInput(
                         "de_volcano_fc",
                         "Volcano |log2FC| cutoff line:",
                         min = 0.5,
                         max = 5,
                         value = 1,
                         step = 0.25
                       ),
                       sliderInput(
                         "de_volcano_pval",
                         "Volcano adj. p-value cutoff line:",
                         min = 0.001,
                         max = 0.1,
                         value = 0.05,
                         step = 0.005
                       ),
                       hr(),
                       sliderInput(
                         "de_legend_size",
                         "Legend dot size:",
                         min = 3,
                         max = 12,
                         value = 6,
                         step = 1
                       ),
                       downloadButton("dl_de_table", "Download Results (CSV)"),
                       downloadButton("dl_de_sig", "Download Significant (CSV)"),
                       downloadButton("dl_de_plot", "Download Plot (PDF)"),
                       hr(),
                       helpText(
                         tags$b("Differential Expression Analysis"),
                         br(),
                         br(),
                         tags$b("Group vs Group:"),
                         br(),
                         "Compares two groups of samples: cell line tissues from the reference,",
                         "or patient groups from the selected query dataset.",
                         br(),
                         tags$b("Wilcoxon:"),
                         "Rank-sum test on Salmon transcript TPM. Works for any two groups,",
                         "including a cell line tissue against a patient group.",
                         br(),
                         tags$b("DESeq2 Gene:"),
                         "Negative binomial model on gene counts summarised from Salmon with tximport.",
                         br(),
                         tags$b("DESeq2 Transcript:"),
                         "Same model on Salmon transcript counts.",
                         br(),
                         "DESeq2 results are precomputed within each dataset, using that dataset's design,",
                         "so both groups must come from the same dataset.",
                         br(),
                         br(),
                         tags$b("DESeq2 designs:"),
                         br(),
                         tags$b("CCLE:"), "~ tissue; 319 of 325 tissue pairs (the 6 pairs among BILIARY, CERVIX,",
                         "SALIVARY GLAND and SMALL INTESTINE have too few samples).",
                         br(),
                         tags$b("Liver:"), "Tumor vs Adjacent, paired by patient (6 patients); one flagged",
                         "outlier sample excluded.",
                         br(),
                         tags$b("Prostate:"), "~ group; all 6 pairs of the 4 groups.",
                         br(),
                         tags$b("Glioblastoma:"), "Post vs Pre G207, paired by patient (3 patients). Pre and post",
                         "samples differ in library preparation and sequencing platform, so differences",
                         "are confounded with these technical effects; treat results as exploratory only.",
                         br(),
                         br(),
                         tags$b("Sample vs Sample:"),
                         br(),
                         "Compares two individual samples. Statistical testing is not possible with n=1,",
                         "so only log2 fold change is computed and results are ranked by its magnitude (MA plot)."
                       )
                     ),
                     mainPanel(
                       width = 9,
                       h4("Results"),
                       uiOutput("gbm_banner_de"),
                       uiOutput("de_method_info"),
                       plotlyOutput("de_plot", height = "500px"),
                       hr(),
                       h4("Transcript Table"),
                       DTOutput("de_table")
                     )
                   ))

tab_kegg <- tabPanel("Pathway Enrichment",
                     sidebarLayout(
                       sidebarPanel(
                         width = 3,
                         radioButtons(
                           "kegg_source",
                           "Gene Source:",
                           choices = c(
                             "DE results (current)" = "de",
                             "Switching results (current)" = "switch"
                           ),
                           selected = "de"
                         ),
                         conditionalPanel(
                           condition = "input.kegg_source == 'de'",
                           sliderInput(
                             "kegg_fc_cutoff",
                             "|log2FC| cutoff:",
                             min = 0,
                             max = 5,
                             value = 1,
                             step = 0.25
                           ),
                           sliderInput(
                             "kegg_pval_cutoff",
                             "Adj. p-value cutoff:",
                             min = 0.001,
                             max = 0.1,
                             value = 0.05,
                             step = 0.005
                           )
                         ),
                         conditionalPanel(
                           condition = "input.kegg_source == 'switch'",
                           helpText(
                             "Uses all genes with detected isoform switching events. Run 'Detect Switches' first."
                           )
                         ),
                         if (HAS_TX_ANNOT)
                           selectInput(
                             "kegg_biotype_filter",
                             "Filter by Biotype:",
                             choices = c("All", biotype_list),
                             selected = "All"
                           ),
                         actionButton("run_kegg", "Run Pathway Enrichment", class = "btn-primary"),
                         hr(),
                         sliderInput(
                           "kegg_top_n",
                           "Top N pathways to display:",
                           min = 5,
                           max = 50,
                           value = 20,
                           step = 5
                         ),
                         downloadButton("dl_kegg_table", "Download Pathway Table (CSV)"),
                         downloadButton("dl_kegg_plot", "Download Plots (PDF)"),
                         hr(),
                         helpText(
                           tags$b("Two modes:"),
                           br(),
                           br(),
                           tags$b("From DE results:"),
                           "Takes significantly differentially expressed genes (above the fold-change and p-value cutoffs)",
                           "and splits them by direction: genes upregulated in Group 1 and genes upregulated in Group 2.",
                           "Runs separate KEGG enrichments for each direction, so you can see which pathways are",
                           "specifically activated in each group. Requires running Group vs Group DE first.",
                           br(),
                           br(),
                           tags$b("From Switching results:"),
                           "Takes all genes where the dominant isoform differs between the two groups.",
                           "Unlike DE, switching is not directional (it is not up or down, just a different isoform),",
                           "so all switching genes are pooled into a single enrichment.",
                           "This reveals pathways where alternative isoform usage is concentrated.",
                           br(),
                           br(),
                           tags$b("Gene universe:"),
                           "Both modes use a proper background: all tested genes from the source analysis,",
                           "not the whole genome. This prevents inflated significance.",
                           br(),
                           br(),
                           tags$b("Reading the dot plot:"),
                           "X-axis = gene ratio (fraction of pathway genes in your list).",
                           "Dot size = number of genes. Dot color = adjusted p-value (darker = more significant).",
                           "Click a row in the table below to see which of your genes belong to that pathway."
                         )
                       ),
                       mainPanel(
                         width = 9,
                         uiOutput("gbm_banner_kegg"),
                         uiOutput("kegg_source_info"),
                         hr(),
                         h4(textOutput("kegg_up_label")),
                         plotlyOutput("kegg_dot_up", height = "400px"),
                         hr(),
                         h4(textOutput("kegg_down_label")),
                         plotlyOutput("kegg_dot_down", height = "400px"),
                         hr(),
                         h4("Pathway Table (combined)"),
                         DTOutput("kegg_table"),
                         hr(),
                         h4("Genes in Selected Pathway"),
                         DTOutput("kegg_gene_table")
                       )
                     ))

tab_drug <- tabPanel("Drug Response", sidebarLayout(
  sidebarPanel(
    width = 3,
    radioButtons(
      "drug_dataset",
      "GDSC Dataset:",
      choices = c("GDSC1" = "gdsc1", "GDSC2" = "gdsc2"),
      selected = "gdsc2"
    ),
    selectizeInput(
      "drug_select",
      "Select Drug:",
      choices = NULL,
      options = list(placeholder = "Type a drug name...")
    ),
    radioButtons(
      "drug_metric",
      "Coloring:",
      choices = c(
        "AUC (continuous)" = "auc",
        "Binary (per-tissue median)" = "binary",
        "Binary (global median)" = "binary_global"
      ),
      selected = "auc"
    ),
    radioButtons(
      "drug_reduction",
      "Embedding:",
      choices = c("UMAP" = "umap", "PCA" = "pca"),
      selected = "umap"
    ),
    if (HAS_GENE_TPM)
      radioButtons(
        "drug_level",
        "Expression Level:",
        choices = c("Isoform (Salmon TPM)" = "isoform", "Gene (tximport TPM)" = "gene"),
        selected = "isoform"
      ),
    checkboxInput("drug_3d", "3D plot", value = FALSE),
    checkboxInput("drug_split", "Split by response (side by side)", value = FALSE),
    sliderInput(
      "drug_top_var",
      "Top N variable isoforms:",
      min = 50,
      max = 5000,
      value = 500,
      step = 50
    ),
    selectizeInput(
      "drug_tissue_filter",
      "Filter by Tissue(s):",
      choices = tissue_choices,
      selected = NULL,
      multiple = TRUE,
      options = list(placeholder = "All tissues (click to filter)")
    ),
    helpText("Leave empty for all. Select tissues to restrict the embedding. At least 3 cell lines are needed for PCA and 5 for UMAP."),
    actionButton("run_drug", "Run", class = "btn-primary"),
    hr(),
    sliderInput(
      "drug_legend_size",
      "Legend dot size:",
      min = 3,
      max = 12,
      value = 6,
      step = 1
    ),
    downloadButton("dl_drug_plot", "Download Plot (PDF)"),
    downloadButton("dl_drug_table", "Download Table (CSV)"),
    hr(),
    helpText(
      "Cell lines embedded by isoform expression, colored by GDSC drug response.",
      "Query samples have no drug data; they are projected onto the same embedding to show which cell lines they resemble.",
      "Overlaying drug response on an embedding is hypothesis generating and does not by itself show isoform-driven response.",
      br(),
      br(),
      tags$b("GDSC1 vs GDSC2:"),
      "GDSC1 covers more compounds; GDSC2 uses updated screening protocols",
      "and different maximum drug concentrations for some drugs.",
      "AUC values between the two datasets are not directly comparable",
      "and any cross-dataset comparison should be interpreted with caution.",
      br(),
      br(),
      tags$b("Concentration range:"),
      "AUC and IC50 depend on the range of concentrations screened, and the maximum tested",
      "concentration (MAX_CONC) can differ between drugs and between GDSC1 and GDSC2.",
      "An LN(IC50) above ln(MAX_CONC) is extrapolated beyond the highest dose tested;",
      "the table flags these cases. Compare cell lines within one drug rather than across drugs.",
      br(),
      br(),
      tags$b("AUC:"),
      "Low = sensitive, high = resistant.",
      br(),
      br(),
      tags$b("Per-tissue binary:"),
      "Splits at median AUC within each tissue.",
      br(),
      br(),
      tags$b("Global binary:"),
      "Splits at overall median AUC."
    )
  ),
  mainPanel(
    width = 9,
    h4("Drug Response on Isoform Embedding"),
    plotlyOutput("drug_plot", height = "600px"),
    hr(),
    h4("Cell Line Response Table"),
    DTOutput("drug_table")
  )
))

tab_circos <- tabPanel("Circos", sidebarLayout(
  sidebarPanel(
    width = 3,
    selectInput(
      "circos_tissue",
      "Expression Tissue:",
      choices = tissue_choices,
      selected = tissue_list[1]
    ),
    selectInput(
      "circos_chr_detail",
      "Chromosome Detail:",
      choices = c("None", paste0("chr", 1:22), "chrX", "chrY"),
      selected = "None"
    ),
    checkboxInput("circos_show_switches", "Overlay switching events", value = FALSE),
    checkboxInput(
      "circos_tissue_specific",
      "Show tissue-specific isoforms only",
      value = FALSE
    ),
    checkboxInput(
      "circos_show_arcs",
      "Show arcs between switching gene pairs",
      value = FALSE
    ),
    helpText(
      style = "color: #888; font-size: 12px;",
      "Switching overlays require running 'Detect Switches' first."
    ),
    selectInput(
      "circos_palette",
      "Expression Color Scheme:",
      choices = c(
        "Default (blue-teal-pink)" = "default",
        "Viridis" = "viridis",
        "Heat (yellow-orange-red)" = "heat",
        "Cool (white-blue-navy)" = "cool"
      ),
      selected = "default"
    ),
    sliderInput(
      "circos_window_kb",
      "Window size (kb):",
      min = 100,
      max = 5000,
      value = 500,
      step = 100
    ),
    sliderInput(
      "circos_low_thresh",
      "Low/moderate boundary (% of max):",
      min = 10,
      max = 60,
      value = 40,
      step = 5
    ),
    sliderInput(
      "circos_high_thresh",
      "Moderate/high boundary (% of max):",
      min = 50,
      max = 90,
      value = 70,
      step = 5
    ),
    sliderInput(
      "circos_top_switch",
      "Max switching genes to plot:",
      min = 10,
      max = 200,
      value = 50,
      step = 10
    ),
    actionButton("run_circos", "Generate Circos Plot", class = "btn-primary"),
    hr(),
    downloadButton("dl_circos", "Download Plot (PDF)"),
    hr(),
    helpText(
      tags$b("Circos Plot"),
      br(),
      "Genome-wide view of isoform expression for the selected tissue or group.",
      br(),
      br(),
      tags$b("Outer track (bars):"),
      "Average isoform expression per genomic window.",
      br(),
      tags$b("Inner track (dots):"),
      "Isoform switching events (run Detect Switches first).",
      br(),
      tags$b("Arcs:"),
      "Connect the top switching genes in rank order by switch score.",
      "Thicker arcs = stronger switching between the two groups.",
      "These show genomic spread of switching events, not functional relationships.",
      br(),
      tags$b("Chromosome Detail:"),
      "Select a chromosome to see an interactive bar chart."
    )
  ),
  mainPanel(
    width = 9,
    h4("Genome-Wide Isoform Landscape"),
    plotOutput("circos_plot", height = "750px", width = "750px"),
    hr(),
    h4("Chromosome Detail (hover for isoform info)"),
    plotlyOutput("circos_chr_plot", height = "350px"),
    hr(),
    h4("Top Switching Genes (Plotted as Arcs)"),
    helpText(
      tags$b("What this shows: "),
      "The genes drawn as arcs in the Circos plot above, ranked by switch score. ",
      "Each row is one gene with its chromosome, genomic position (Mb), the switch score, ",
      "and the dominant isoform in each of the two compared groups.",
      br(),
      br(),
      tags$b("How to use: "),
      "Set your tissue pair in the Isoform Switching tab and click 'Detect Switches' first. ",
      "Adjust 'Max switching genes to plot' in the sidebar to change how many genes appear here and in the Circos arcs."
    ),
    DTOutput("circos_tab10a"),
    downloadButton("dl_circos_tab10a", "Download Switching Genes Table (CSV)"),
    hr(),
    h4("Per-Chromosome Summary"),
    helpText(
      tags$b("What this shows: "),
      "One row per chromosome (chr1-22, X, Y) summarizing the data behind the Circos outer track. ",
      "For the Expression Tissue selected in the sidebar: the average TPM across all transcripts on that chromosome, ",
      "the total transcript count, and the number of switching events detected on that chromosome.",
      br(),
      br(),
      tags$b("How to use: "),
      "Change the Expression Tissue dropdown in the sidebar to update this table for a different tissue. ",
      "Switching event counts require 'Detect Switches' to have been run in the Isoform Switching tab."
    ),
    DTOutput("circos_tab10b"),
    downloadButton("dl_circos_tab10b", "Download Chromosome Summary (CSV)")
  )
))

tab_about <- tabPanel("About", fluidRow(
  column(
    8,
    offset = 2,
    div(
      style = "text-align: center; margin-top: 30px; margin-bottom: 20px;",
      isoscope_logo(120),
      h2("IsoScope", style = "margin-top: 0; color: #225ea8; font-weight: bold; letter-spacing: 1px;"),
      p(
        style = "color: #41b6c4; font-size: 14px; font-weight: 600; margin-bottom: 0;",
        paste0(ref_row$n_samples, " cell line samples \u2022 ", ref_row$n_groups, " tissues \u2022 ",
               length(QUERY_IDS), " patient cohorts \u2022 Transcript-level resolution")
      ),
      h4(
        "Interactive visualization of transcript isoform expression across cancer cell lines and patient cohorts",
        style = "color: #666; font-weight: normal; margin-top: 10px;"
      )
    ),
    hr(),
    h4("Why Isoform-Level Analysis?"),
    p(
      "Gene-level expression summaries collapse all transcript variants into a single measurement,",
      "masking biologically important variation. Different isoforms of the same gene can have distinct",
      "or opposing functions, and cancer cells frequently exploit alternative isoforms to evade apoptosis,",
      "promote proliferation, or resist treatment."
    ),
    h4("Data"),
    p(
      "RNA-seq data is from the Cancer Cell Line Encyclopedia (CCLE), BioProject",
      tags$a(
        href = "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA523380",
        target = "_blank",
        "PRJNA523380",
        style = "color: #225ea8;"
      ),
      ". Raw reads were aligned to GRCh38/hg38 using STAR 2.7.11b with GENCODE v49.",
      "Transcript-level expression was quantified with Salmon 1.10.0 in alignment-based mode against",
      "GENCODE v49 basic transcripts. CCLE comprises 668 samples from 667 cell lines (the TT cell line",
      "is represented by separate oesophagus and thyroid samples).",
      "Patient cohorts (glioblastoma, hepatocellular carcinoma, prostate adenocarcinoma) were processed",
      "from raw reads with STAR alignment and Salmon quantification against the same reference.",
      "Gene-level values were summarised with tximport (lengthScaledTPM), and differential expression",
      "was computed with DESeq2 1.52.0 separately within each dataset. These results were precomputed",
      "in R 4.6; the hosted application itself runs the R version shown at the bottom of this page.",
      "Drug response data (AUC, IC50) from GDSC1 and GDSC2 via",
      tags$a(
        href = "https://cellmodelpassports.sanger.ac.uk/downloads",
        target = "_blank",
        "Cell Model Passports",
        style = "color: #225ea8;"
      ),
      "."
    ),
    p(
      tags$b("Patient cohort notes."),
      "Hepatocellular carcinoma: DESeq2 uses the 6 patients with paired tumor and adjacent samples;",
      "one adjacent sample flagged as an outlier in the source metadata is excluded.",
      "Glioblastoma: pre- and post-G207 samples differ in library preparation, read length and",
      "sequencing platform, so differences between them are confounded with these technical effects",
      "and should be treated as exploratory only."
    ),
    p(
      "Transcript annotations (biotype, exon count, transcript length) from GENCODE v49 GTF.",
      "Protein domain annotations (Pfam, InterPro) and UniProt mappings from Ensembl BioMart release 115."
    ),
    
    tags$ol(
      style = "margin-top: 8px; color: #666; font-size: 13px; padding-left: 20px;",
      tags$li(
        "Ghandi, M. et al.",
        tags$i(
          "Next-generation characterization of the Cancer Cell Line Encyclopedia."
        ),
        "Nature, 569, 503\u2013508, 2019."
      ),
      tags$li(
        "Garnett, M.J. et al.",
        tags$i(
          "Systematic identification of genomic markers of drug sensitivity in cancer cells."
        ),
        "Nature, 483, 570\u2013575, 2012."
      ),
      tags$li(
        "Yang, W. et al.",
        tags$i("Genomics of Drug Sensitivity in Cancer (GDSC): a resource for therapeutic biomarker discovery in cancer cells."),
        "Nucleic Acids Research, 41, D955\u2013D961, 2013."
      ),
      tags$li(
        "Iorio, F. et al.",
        tags$i("A landscape of pharmacogenomic interactions in cancer."),
        "Cell, 166, 740\u2013754, 2016."
      ),
      tags$li(
        "van der Meer, D. et al.",
        tags$i("Cell Model Passports: a hub for clinical, genetic and functional datasets of preclinical cancer models."),
        "Nucleic Acids Research, 47, D923\u2013D929, 2019."
      ),
      tags$li(
        "Dobin, A. et al.",
        tags$i("STAR: ultrafast universal RNA-seq aligner."),
        "Bioinformatics, 29, 15\u201321, 2013."
      ),
      tags$li(
        "Frankish, A. et al.",
        tags$i("GENCODE 2021."),
        "Nucleic Acids Research, 49, D916\u2013D923, 2021."
      ),
      tags$li(
        "McInnes, L., Healy, J., Melville, J.",
        tags$i("UMAP: Uniform Manifold Approximation and Projection for dimension reduction."),
        "arXiv:1802.03426, 2018."
      ),
      tags$li(
        "Kanehisa, M., Goto, S.",
        tags$i("KEGG: Kyoto Encyclopedia of Genes and Genomes."),
        "Nucleic Acids Research, 28, 27\u201330, 2000."
      ),
      tags$li(
        "Yu, G., Wang, L.-G., Han, Y., He, Q.-Y.",
        tags$i("clusterProfiler: an R package for comparing biological themes among gene clusters."),
        "OMICS, 16, 284\u2013287, 2012."
      ),
      tags$li(
        "Patro, R., Duggal, G., Love, M.I., Irizarry, R.A., Kingsford, C.",
        tags$i("Salmon provides fast and bias-aware quantification of transcript expression."),
        "Nature Methods, 14, 417\u2013419, 2017."
      ),
      tags$li(
        "Soneson, C., Love, M.I., Robinson, M.D.",
        tags$i("Differential analyses for RNA-seq: transcript-level estimates improve gene-level inferences."),
        "F1000Research, 4, 1521, 2015."
      ),
      tags$li(
        "Love, M.I., Huber, W., Anders, S.",
        tags$i(
          "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2."
        ),
        "Genome Biology, 15, 550, 2014."
      ),
      tags$li(
        "Martin, F.J. et al.",
        tags$i("Ensembl 2023."),
        "Nucleic Acids Research, 51, D933\u2013D941, 2023."
      ),
      tags$li(
        "The UniProt Consortium.",
        tags$i("UniProt: the Universal Protein Knowledgebase in 2023."),
        "Nucleic Acids Research, 51, D523\u2013D531, 2023."
      ),
      tags$li(
        "Mistry, J. et al.",
        tags$i("Pfam: The protein families database in 2021."),
        "Nucleic Acids Research, 49, D412\u2013D419, 2021."
      ),
      tags$li(
        "Paysan-Lafosse, T. et al.",
        tags$i("InterPro in 2022."),
        "Nucleic Acids Research, 51, D418\u2013D427, 2023."
      ),
      tags$li(
        "Markert, J.M. et al.",
        tags$i(
          "Phase Ib trial of mutant herpes simplex virus G207 inoculated pre- and post-tumor resection for recurrent GBM."
        ),
        "Molecular Therapy, 17, 199\u2013207, 2009.",
        "RNA-seq data:",
        tags$a(
          href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162643",
          target = "_blank",
          "GSE162643",
          style = "color: #225ea8;"
        )
      ),
      tags$li(
        "Patient cohorts: hepatocellular carcinoma, SRA BioProject",
        tags$a(href = "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1427671", target = "_blank",
               "PRJNA1427671", style = "color: #225ea8;", .noWS = "after"),
        "; prostate adenocarcinoma, SRA BioProject",
        tags$a(href = "https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1497887", target = "_blank",
               "PRJNA1497887", style = "color: #225ea8;", .noWS = "after"),
        "."
      )
    ),
    h4("Citation"),
    tags$blockquote(
      style = "border-left: 3px solid #225ea8; padding-left: 12px; color: #444;",
      "Sajjad A, Martinez S, Johnson L, Ballesteros Prieto V, Arestakesyan H, Dias J, Alhossiny M, and Horvath A. IsoScope: Interactive Pan-Cancer Isoform-Level Analysis. 2026."
    ),
    h4("Horvath Lab"),
    p(
      "This application was developed in the Horvath Lab at The George Washington University, Department of Biochemistry and Molecular Medicine.",
    ),
    h4("Data Availability"),
    p(
      "IsoScope is publicly available at ",
      tags$a(
        href = "https://huggingface.co/spaces/HorvathLab/IsoScope",
        target = "_blank",
        "https://huggingface.co/spaces/HorvathLab/IsoScope",
        style = "color: #225ea8;"
      ),
      ". Source code and preprocessing scripts are available at ",
      tags$a(
        href = "https://github.com/HorvathLab/IsoScope",
        target = "_blank",
        "https://github.com/HorvathLab/IsoScope",
        style = "color: #225ea8;"
      ),
      "."
    ),
    p(
      style = "color: #999; margin-top: 30px; text-align: center;",
      paste0(
        "IsoScope v2.0 | Built with R ",
        R.version$major,
        ".",
        R.version$minor,
        " + Shiny"
      )
    )
  )
))

ui <- function(request) {
  navbarPage(
    title = tags$span(
      isoscope_logo_light(48),
      tags$span(style = "font-weight: 800; font-size: 22px; margin-left: 8px; letter-spacing: 0.5px; color: white;", "IsoScope")
    ),
    id = "main_nav",
    theme = bs_theme(
      bg = "#ffffff",
      fg = "#333333",
      primary = "#0d3b6e",
      "navbar-bg" = "#0d3b6e",
      "nav-link-font-size" = "1rem",
      "nav-link-padding-y" = "0.75rem"
    ),
    tabPanel("Dataset", icon = icon("database"), fluidRow(column(
      10, offset = 1, div(
        style = "margin-top: 30px;",
        div(
          class = "ds-panel",
          style = "position: relative; overflow: hidden;",
          style = "background: #E6F1FB; border-radius: 16px; padding: 32px;",
          tags$style(HTML("
            .ds-panel .shiny-options-group { display: grid !important; grid-template-columns: repeat(auto-fit, minmax(230px, 1fr)); gap: 18px; width: 100%; }
            .ds-panel .radio-inline, .ds-panel .form-check { display: block; width: 100%; margin: 0 !important; padding: 0 !important; }
            .ds-panel .radio-inline label, .ds-panel .form-check label, .ds-panel label.radio-inline { display: block; width: 100%; height: 100%; margin: 0; }
            .ds-panel input[type=radio] { position: absolute; opacity: 0; pointer-events: none; }
            .ds-panel .radio-inline > span, .ds-panel .form-check label > span { display: block; height: 100%; }
            .ds-panel .ds-card { height: 100%; box-sizing: border-box; }
            .ds-panel .ds-card { cursor: pointer; }
            .ds-panel .radio-inline:hover .ds-card { box-shadow: 0 4px 12px rgba(4,44,83,0.14); }
            .ds-panel label:has(input[type=radio]:checked) .ds-card { outline: 3px solid #378ADD; outline-offset: -2px; box-shadow: 0 4px 14px rgba(55,138,221,0.25); }
            .ds-panel label:has(input[type=radio]:checked) .ds-card::after { content: 'SELECTED'; display: block; margin-top: 12px; font-size: 13px; letter-spacing: 0.07em; color: #378ADD; font-weight: 600; }
          ")),
          HTML('<svg viewBox="0 0 260 220" preserveAspectRatio="xMidYMid meet" style="position:absolute;top:28px;right:28px;width:240px;height:auto;opacity:0.9;pointer-events:none;" aria-hidden="true"><g fill="#E8B84B"><rect x="8" y="26" width="240" height="2" opacity="0.4"/><rect x="8" y="19" width="34" height="16" rx="3"/><rect x="76" y="19" width="52" height="16" rx="3"/><rect x="168" y="19" width="40" height="16" rx="3"/><rect x="8" y="82" width="240" height="2" opacity="0.4"/><rect x="8" y="75" width="34" height="16" rx="3" opacity="0.85"/><rect x="168" y="75" width="40" height="16" rx="3" opacity="0.85"/><rect x="8" y="138" width="240" height="2" opacity="0.4"/><rect x="76" y="131" width="52" height="16" rx="3" opacity="0.7"/><rect x="168" y="131" width="80" height="16" rx="3" opacity="0.7"/><rect x="8" y="194" width="240" height="2" opacity="0.4"/><rect x="8" y="187" width="34" height="16" rx="3" opacity="0.55"/><rect x="76" y="187" width="52" height="16" rx="3" opacity="0.55"/><rect x="168" y="187" width="40" height="16" rx="3" opacity="0.55"/></g></svg>'),
          div(
            style = "font-size: 16px; font-weight: 600; color: #185FA5; letter-spacing: 0.1em; margin-bottom: 10px;",
            "ISOSCOPE"
          ),
          h3(style = "font-size: 40px; font-weight: 500; color: #042C53; margin: 0 0 10px;", "Pan-cancer isoform explorer"),
          p(
            style = "font-size: 19px; color: #185FA5; margin: 0 0 28px; max-width: 46em;",
            "Transcript-level expression across cancer cell lines and patient cohorts."
          ),
          div(
            style = "display: flex; gap: 40px; padding-bottom: 24px; border-bottom: 1px solid #B5D4F4; margin-bottom: 24px;",
            lapply(list(
              list(sum(DATASETS$n_samples), "samples"),
              list(format(ref_row$n_transcripts, big.mark = ","), "transcripts in reference"),
              list(nrow(DATASETS), "datasets")
            ), function(s) div(
              div(style = "font-size: 36px; font-weight: 500; color: #042C53;", s[[1]]),
              div(style = "font-size: 16px; color: #185FA5; margin-top: 6px;", s[[2]])
            ))
          ),
          radioButtons(
            "dataset_mode", NULL,
            choiceNames = lapply(seq_len(nrow(DATASETS)), function(i) {
              d <- DATASETS[i, ]
              is_ref <- d$is_reference == 1
              accent <- if (is_ref) "#378ADD" else c("#1D9E75", "#7F77DD", "#D85A30")[((i - 1) %% 3) + 1]
              caveat <- if (!is.na(d$mapping_rate) && nzchar(as.character(d$mapping_rate)))
                paste0("Remapped from ", d$remapped_from, ", ", round(as.numeric(d$mapping_rate) * 100, 1), "% retained")
              else NULL
              div(
                class = "ds-card",
                style = paste0("background: #ffffff !important; border-radius: 14px; padding: 22px; border-top: 5px solid ", accent, "; min-height: 230px; box-shadow: 0 1px 3px rgba(4,44,83,0.08); transition: box-shadow 0.15s;"),
                div(style = paste0("font-size: 13px; font-weight: 600; letter-spacing: 0.08em; color: ", accent, "; margin-bottom: 12px;"),
                    if (is_ref) "REFERENCE" else "QUERY"),
                div(style = "font-size: 19px; font-weight: 500; margin-bottom: 16px;", d$display_name),
                div(style = "font-size: 44px; font-weight: 500; line-height: 1;", d$n_samples),
                div(style = "font-size: 16px; color: #666; margin: 6px 0 14px;", "samples"),
                div(style = "font-size: 16px; color: #555; padding-top: 12px; border-top: 1px solid #eee; line-height: 1.6;",
                    if (is_ref) {
                      tagList(paste0(d$n_subjects, " cell lines"), tags$br(),
                              paste0(d$n_groups, " tissues"))
                    } else {
                      g <- QUERY_GROUPS[QUERY_GROUPS$dataset == d$dataset_id, ]
                      g <- g[order(-g$n), ]
                      g_lbl <- ifelse(g$raw %in% names(GROUP_FULL), GROUP_FULL[g$raw], g$raw)
                      tagList(
                        tags$div(style = "font-weight: 500; color: #333; margin-bottom: 4px;",
                                 paste0(d$n_subjects, " patients")),
                        lapply(seq_len(nrow(g)), function(k)
                          tags$div(tags$span(style = "font-weight: 600; color: #042C53; display: inline-block; min-width: 2.2em;", g$n[k]),
                                   g_lbl[k]))
                      )
                    },
                    if (!is.null(caveat)) tagList(tags$br(), tags$span(style = "color: #993C1D;", caveat))
                )
              )
            }),
            choiceValues = DATASETS$dataset_id,
            selected = REF_ID,
            inline = TRUE
          )
        ),
        div(
          style = "margin-top: 24px;",
          h4(style = "font-size: 28px; font-weight: 500; margin-bottom: 16px;", "How datasets are used"),
          p(style = "color: #444; font-size: 18px; line-height: 1.75;",
            "Selecting a query dataset loads it alongside the reference. Each dataset keeps its own database,",
            "so nothing is merged and the reference is never altered."
          ),
          p(style = "color: #444; font-size: 18px; line-height: 1.75;",
            "In embedding views, query samples are projected onto the reference coordinate space, showing where",
            "individual tumours sit relative to established cell lines without distorting the reference landscape.",
            "Drug response data is available for cell lines only, but the spatial position of patient samples",
            "relative to sensitive and resistant lines still carries translational information."
          ),
          p(style = "color: #444; font-size: 18px; line-height: 1.75;",
            "Differential expression is run within each dataset using its own design, so effect sizes are not",
            "directly comparable between datasets. New datasets can be added by building a database to the",
            "published schema and placing it in the data folder."
          )
        )
      )
    ))),
    tab_gene_viewer,
    tab_switching,
    tab_diversity,
    tab_tissue_specific,
    tab_clustering,
    tab_corr_heatmap,
    tab_drug,
    tab_de,
    tab_kegg,
    tab_circos,
    tab_qc,
    tab_about
  )
}

# --- Server ---

server <- function(input, output, session) {
  onStop(function() {
    if (DBI::dbIsValid(db))
      DBI::dbDisconnect(db)
    if (!is.null(deseq_db) &&
        DBI::dbIsValid(deseq_db))
      DBI::dbDisconnect(deseq_db)
  })
  
  updateSelectizeInput(
    session,
    "gene",
    choices = gene_list,
    selected = default_gene,
    server = TRUE
  )
  updateSelectizeInput(
    session,
    "div_gene",
    choices = gene_list,
    selected = default_gene,
    server = TRUE
  )
  updateSelectizeInput(
    session,
    "spec_gene",
    choices = c("All genes" = "", gene_list),
    selected = "",
    server = TRUE
  )
  updateSelectizeInput(
    session,
    "clust_gene",
    choices = c("All genes" = "__ALL__", gene_list),
    selected = "__ALL__",
    server = TRUE
  )
  updateSelectizeInput(
    session,
    "corr_gene",
    choices = c("All genes" = "__ALL__", gene_list),
    selected = "__ALL__",
    server = TRUE
  )
  updateSelectizeInput(
    session,
    "drug_select",
    choices = drug_list,
    selected = drug_list[1],
    server = TRUE
  )
  log_timing <- function(label) {
    t0 <- Sys.time()
    function() {
      mem <- tryCatch(sum(gc()[, 2]), error = function(e) NA)
      cat(sprintf("[timing] %s | dataset %s | %.1f s | R heap %.0f MB\n",
                  label, isolate(active_id()), as.numeric(Sys.time() - t0, units = "secs"), mem))
    }
  }

  merge_quants <- function(ref, q) {
    if (is.null(q) || nrow(q) == 0) return(ref)
    if (nrow(ref) == 0) return(q)
    m <- merge(as.data.table(ref), as.data.table(q), by = c("gene", "transcript_id"), all = TRUE)
    num <- names(m)[vapply(m, is.numeric, logical(1))]
    for (col in num) set(m, which(is.na(m[[col]])), col, 0)
    as.data.frame(m)
  }

  fetch_quants <- function(sql, params = NULL) {
    ref <- dbGetQuery(db, sql, params = params)
    qc <- query_con(active_id())
    if (is.null(qc)) return(ref)
    q <- tryCatch(dbGetQuery(qc, sql, params = params), error = function(e) NULL)
    merge_quants(ref, q)
  }

  quant_fields <- function() {
    qc <- query_con(active_id())
    union(dbListFields(db, 'quants'), if (!is.null(qc)) dbListFields(qc, 'quants'))
  }

  fetch_quant_cols <- function(cols, where = NULL, params = NULL) {
    one <- function(con) {
      keep <- intersect(cols, dbListFields(con, "quants"))
      sel <- paste(paste0('"', c("gene", "transcript_id", keep), '"'), collapse = ", ")
      sql <- paste0("SELECT ", sel, " FROM quants", if (!is.null(where)) paste0(" WHERE ", where))
      if (is.null(params)) dbGetQuery(con, sql) else dbGetQuery(con, sql, params = params)
    }
    out <- one(db)
    qc <- query_con(active_id())
    if (!is.null(qc)) out <- merge_quants(out, one(qc))
    for (col in setdiff(cols, names(out))) out[[col]] <- 0
    out[, c("gene", "transcript_id", cols), drop = FALSE]
  }

  fetch_union <- function(sql) {
    ref <- dbGetQuery(db, sql)
    qc <- query_con(active_id())
    if (is.null(qc)) return(ref)
    q <- tryCatch(dbGetQuery(qc, sql), error = function(e) NULL)
    if (is.null(q)) return(ref)
    as.data.frame(rbindlist(list(ref, q), fill = TRUE))
  }

  active_id <- reactive({
    if (is.null(input$dataset_mode)) REF_ID else input$dataset_mode
  })

  gbm_banner <- function() renderUI({
    if (!identical(active_id(), "gbm")) return(NULL)
    div(style = "background: #fdf1e7; border-left: 4px solid #D85A30; padding: 10px 14px; margin-bottom: 12px; font-size: 14px;",
        tags$b("Glioblastoma: "),
        "pre- and post-G207 samples differ in library preparation (total RNA vs poly(A)), read length and",
        "sequencing platform. Differences between these groups are confounded with these technical effects",
        "and should be treated as exploratory only.")
  })
  output$gbm_banner_de <- gbm_banner()
  output$gbm_banner_kegg <- gbm_banner()

  query_samples <- reactive({
    id <- active_id()
    if (id == REF_ID) return(NULL)
    row <- DATASETS[DATASETS$dataset_id == id, ]
    if (nrow(row) != 1) return(NULL)
    con <- dbConnect(RSQLite::SQLite(), row$db_path)
    on.exit(dbDisconnect(con))
    s <- dbReadTable(con, "samples")
    s$Tissue   <- pretty_group(s$group_label, id)
    s$CellLine <- s$subject_id
    s$dataset  <- id
    s
  })

  active_meta <- reactive({
    q <- query_samples()
    if (is.null(q)) {
      m <- sample_meta_df
      m$dataset <- REF_ID
      m
    } else {
      m <- sample_meta_df
      m$dataset <- REF_ID
      cols <- intersect(names(m), names(q))
      rbind(m[, cols, drop = FALSE], q[, cols, drop = FALSE])
    }
  })

  active_tissues <- reactive({
    sort(unique(active_meta()$Tissue))
  })

  query_tissues <- reactive({
    q <- query_samples()
    if (is.null(q)) character(0) else sort(unique(q$Tissue))
  })

  active_avg <- reactive({
    get_active_avg(active_id())
  })

  active_cell_lines <- reactive({
    sort(unique(active_meta()$CellLine))
  })

  
  observeEvent(input$dataset_mode, {
    at <- active_tissues()
    am <- active_meta()
    tc <- as.integer(table(am$Tissue)[at])
    ds <- am$dataset[match(at, am$Tissue)]
    is_q <- !(at %in% CCLE_TISSUES)
    o <- order(!is_q, -tc)
    at <- at[o]; tc <- tc[o]; ds <- ds[o]
    new_choices <- setNames(at, paste0(at, " (n=", tc, ")"))
    for (id in c(
      "tissue_filter",
      "switch_tissue1",
      "switch_tissue2",
      "de_tissue1",
      "de_tissue2",
      "spec_tissue",
      "circos_tissue"
    )) {
      updateSelectInput(session, id, choices = c(if (id == "tissue_filter")
        c("All tissues (averaged)" = "All"), new_choices),
        selected = if (id %in% c("switch_tissue2", "de_tissue2") && length(new_choices) > 1)
          new_choices[[2]] else NULL)
    }
    for (id in c(
      "gene_tissue_multi",
      "clust_tissue_filter",
      "corr_tissue_filter",
      "drug_tissue_filter",
      "qc_tissue",
      "div_tissue_filter"
    )) {
      updateSelectizeInput(session, id, choices = new_choices, selected = NULL)
    }
    acl <- active_cell_lines()
    for (id in c("switch_cl1", "de_cl1")) {
      updateSelectizeInput(
        session,
        id,
        choices = acl,
        selected = acl[1],
        server = TRUE
      )
    }
    for (id in c("switch_cl2", "de_cl2")) {
      updateSelectizeInput(
        session,
        id,
        choices = acl,
        selected = acl[min(2, length(acl))],
        server = TRUE
      )
    }
  }, ignoreNULL = FALSE)
  
  observeEvent(input$drug_dataset, {
    if (input$drug_dataset == "gdsc1" && !is.null(gdsc1_response_df))
      updateSelectizeInput(
        session,
        "drug_select",
        choices = drug_list_1,
        selected = drug_list_1[1],
        server = TRUE
      )
    else
      updateSelectizeInput(
        session,
        "drug_select",
        choices = drug_list,
        selected = drug_list[1],
        server = TRUE
      )
  })
  updateSelectizeInput(
    session,
    "switch_cl1",
    choices = cell_line_list,
    selected = cell_line_list[1],
    server = TRUE
  )
  updateSelectizeInput(
    session,
    "switch_cl2",
    choices = cell_line_list,
    selected = cell_line_list[min(2, length(cell_line_list))],
    server = TRUE
  )
  updateSelectizeInput(
    session,
    "de_cl1",
    choices = cell_line_list,
    selected = cell_line_list[1],
    server = TRUE
  )
  updateSelectizeInput(
    session,
    "de_cl2",
    choices = cell_line_list,
    selected = cell_line_list[min(2, length(cell_line_list))],
    server = TRUE
  )
  
  observeEvent(input$de_mode, {
    if (input$de_mode == "cellline" &&
        input$de_method %in% c("deseq2_gene", "deseq2_tx")) {
      updateRadioButtons(session, "de_method", selected = "wilcoxon")
      showNotification(
        "DESeq2 requires biological replicates. Switched to Wilcoxon.",
        type = "warning",
        duration = 5
      )
    }
  })
  
  DE_NOTES <- c(
    ccle = "CCLE, ~ tissue",
    hcc = "Liver, paired by patient",
    pca = "Prostate, ~ group",
    gbm = "Glioblastoma, paired by patient; confounded with library prep"
  )

  de_pair_status <- reactive({
    g <- c(input$de_tissue1, input$de_tissue2)
    if (any(!nzchar(g)) || length(g) < 2) return(list(kind = "none"))
    if (g[1] == g[2]) return(list(kind = "same"))
    in_q <- g %in% query_tissues()
    if (any(in_q) && !all(in_q)) return(list(kind = "cross"))
    id <- if (all(in_q)) active_id() else REF_ID
    list(kind = if (is.null(deseq_con(id))) "nodeseq" else "within", id = id)
  })

  output$de_availability <- renderUI({
    st <- de_pair_status()
    box <- function(bg, fg, txt) tags$div(
      style = paste0("background: ", bg, "; color: ", fg, "; padding: 7px 10px; border-radius: 5px; font-size: 13px; margin: -4px 0 10px;"),
      txt)
    switch(st$kind,
      same = box("#fdecea", "#8a1f11", "Pick two different groups."),
      cross = box("#fff4e5", "#7a4b00", "Groups from different datasets: Wilcoxon only."),
      nodeseq = box("#fff4e5", "#7a4b00", "No DESeq2 results for this dataset: Wilcoxon only."),
      within = box("#e8f4ec", "#1e5631", paste0("DESeq2 gene and transcript available (",
                                              DE_NOTES[st$id], ")")),
      NULL)
  })

  observeEvent(list(input$de_tissue1, input$de_tissue2, input$de_method), {
    if (identical(input$de_mode, "tissue") &&
        input$de_method %in% c("deseq2_gene", "deseq2_tx") &&
        de_pair_status()$kind %in% c("cross", "nodeseq")) {
      updateRadioButtons(session, "de_method", selected = "wilcoxon")
      showNotification("DESeq2 is computed within each dataset. Switched to Wilcoxon.",
                       type = "warning", duration = 6)
    }
  }, ignoreInit = TRUE)

  output$de_cl_tissues <- renderUI({
    req(input$de_mode == "cellline")
    t1 <- active_meta()$Tissue[match(input$de_cl1, active_meta()$CellLine)]
    t2 <- active_meta()$Tissue[match(input$de_cl2, active_meta()$CellLine)]
    tags$div(
      style = "background: #f0f4ff; padding: 6px 10px; border-radius: 5px; font-size: 12px; margin-bottom: 8px;",
      tags$span(style = "font-weight: 600;", input$de_cl1),
      tags$span(style = "color: #666;", paste0(" (", t1, ")")),
      tags$span(style = "margin: 0 4px;", "vs"),
      tags$span(style = "font-weight: 600;", input$de_cl2),
      tags$span(style = "color: #666;", paste0(" (", t2, ")"))
    )
  })
  
  # === Gene Viewer ===
  
  gene_data <- reactive({
    input$run_gene
    ensure_db()
    done <- log_timing("Gene Viewer")
    on.exit(done())
    req(input$gene)
    if (input$tissue_filter == "All") {
      gdata <- active_avg() %>% filter(gene == input$gene)
      if (nrow(gdata) == 0)
        return(NULL)
      tiss_cols <- intersect(names(gdata), active_tissues())
      if (length(input$gene_tissue_multi) > 0)
        tiss_cols <- intersect(tiss_cols, input$gene_tissue_multi)
      bind_rows(lapply(tiss_cols, function(tissue) {
        tibble(
          transcript_id = gdata$transcript_id,
          display_label = gdata$display_label,
          sample = tissue,
          TPM = gdata[[tissue]],
          Tissue = tissue,
          view_type = "tissue"
        )
      }))
    } else {
      df <- fetch_quants(
                       "SELECT * FROM quants WHERE gene = ?",
                       params = list(input$gene))
      if (nrow(df) == 0)
        return(NULL)
      df$display_label <- paste0(df$gene, " | ", df$transcript_id)
      cols <- intersect(
        active_meta() %>% filter(Tissue == input$tissue_filter) %>% pull(sample_label),
        names(df)
      )
      if (length(cols) == 0)
        return(NULL)
      df %>% select(transcript_id, display_label, all_of(cols)) %>%
        pivot_longer(-c(transcript_id, display_label),
                     names_to = "sample",
                     values_to = "TPM") %>%
        mutate(Tissue = input$tissue_filter, view_type = "cell_line")
    }
  })
  
  gene_plot_obj <- reactive({
    df <- gene_data()
    if (is.null(df))
      return(NULL)
    df <- collapse_isoforms(df, input$top_n)
    groups <- unique(df$isoform_group)
    non_other <- groups[groups != "Other"]
    other_present <- "Other" %in% groups
    pal_base <- get_colors(length(non_other))
    pal <- if (other_present) {
      setNames(c(pal_base, "#79aaff"), c(non_other, "Other"))
    } else {
      setNames(pal_base, non_other)
    }
    if (input$log_scale &&
        input$plot_type != "stacked") {
      df$TPM <- log2(df$TPM + 1)
      y_lab <- "log2(TPM + 1)"
    } else {
      y_lab <- "TPM"
    }
    if (input$plot_type == "stacked") {
      df_agg <- df %>% group_by(sample, isoform_group) %>%
        summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")
      totals <- df_agg %>% group_by(sample) %>%
        summarise(total = sum(TPM), .groups = "drop")
      df_agg <- df_agg %>% left_join(totals, by = "sample") %>%
        mutate(
          proportion = ifelse(total > 0, TPM / total, 0),
          sample = factor(sample, levels = sort(unique(
            as.character(sample)
          )))
        )
      ggplot(df_agg,
             aes(
               x = sample,
               y = proportion,
               fill = isoform_group,
               text = paste0(
                 "Sample: ",
                 sample,
                 "\nIsoform: ",
                 isoform_group,
                 "\nProportion: ",
                 round(proportion, 3),
                 "\nTPM: ",
                 round(TPM, 2)
               )
             )) +
        geom_bar(stat = "identity", position = "stack") +
        labs(
          title = paste("Isoform Proportions:", input$gene),
          x = NULL,
          y = "Proportion",
          fill = "Isoform"
        ) +
        theme_minimal() + scale_fill_manual(values = pal) +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            size = max(10, min(14, round(
              300 / length(unique(df_agg$sample))
            )))
          ),
          legend.position = "bottom",
          plot.title = element_text(size = 18, face = "bold"),
          axis.text.y = element_text(size = 13),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 15)
        )
    } else if (input$plot_type == "grouped") {
      df_agg <- df %>% group_by(sample, isoform_group) %>%
        summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")
      ggplot(df_agg,
             aes(
               x = sample,
               y = TPM,
               fill = isoform_group,
               text = paste0(
                 "Sample: ",
                 sample,
                 "\nIsoform: ",
                 isoform_group,
                 "\nTPM: ",
                 round(TPM, 2)
               )
             )) +
        geom_bar(stat = "identity", position = "dodge") +
        labs(
          title = paste("Isoform Expression:", input$gene),
          x = NULL,
          y = y_lab,
          fill = "Isoform"
        ) +
        theme_minimal() + scale_fill_manual(values = pal) +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            size = max(10, min(14, round(
              300 / length(unique(df_agg$sample))
            )))
          ),
          legend.position = "bottom",
          plot.title = element_text(size = 18, face = "bold"),
          axis.text.y = element_text(size = 13),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 15)
        )
    } else {
      df_agg <- df %>% group_by(sample, isoform_group) %>%
        summarise(TPM = sum(TPM, na.rm = TRUE), .groups = "drop")
      ggplot(df_agg,
             aes(
               x = sample,
               y = isoform_group,
               fill = TPM,
               text = paste0(
                 "Sample: ",
                 sample,
                 "\nIsoform: ",
                 isoform_group,
                 "\nTPM: ",
                 round(TPM, 2)
               )
             )) +
        geom_tile() + scale_fill_gradientn(colours = TROPICAL) +
        labs(
          title = paste("Isoform Heatmap:", input$gene),
          x = NULL,
          y = "Isoform",
          fill = y_lab
        ) +
        theme_minimal(base_size = 13) +
        theme(
          axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            size = max(10, min(14, round(
              300 / length(unique(df_agg$sample))
            )))
          ),
          plot.title = element_text(size = 18, face = "bold"),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 14),
          panel.grid = element_blank()
        )
    }
  })
  
  output$isoform_plot <- renderPlotly({
    p <- gene_plot_obj()
    if (is.null(p))
      return(NULL)
    n_shown <- length(unique(gene_data()$sample))
    legend_size <- max(10, min(13, round(280 / n_shown)))
    tick_size <- max(10, min(13, round(300 / n_shown)))
    ggplotly(p, tooltip = "text") %>%
      layout(
        legend = list(
          orientation = "h",
          y = -0.45,
          font = list(size = legend_size),
          title = list(
            text = "Isoform",
            side = "top",
            font = list(size = 12)
          )
        ),
        xaxis = list(tickfont = list(size = tick_size)),
        title = list(font = list(size = 20)),
        margin = list(b = 180)
      )
  })
  output$dl_plot_gene <- downloadHandler(
    filename = function()
      paste0("isoform_", input$gene, "_", input$tissue_filter, ".pdf"),
    content = function(file) {
      p <- gene_plot_obj()
      if (!is.null(p))
        ggsave(file, p, width = 16, height = 7)
    }
  )
  
  gene_table_data <- reactive({
    req(input$gene)
    if (input$tissue_filter == "All") {
      tcols <- intersect(names(active_avg()), active_tissues())
      if (length(input$gene_tissue_multi) > 0)
        tcols <- intersect(tcols, input$gene_tissue_multi)
      result <- active_avg() %>% filter(gene == input$gene) %>%
        select(gene, transcript_id, all_of(tcols))
    } else {
      df <- fetch_quants(
                       "SELECT * FROM quants WHERE gene = ?",
                       params = list(input$gene))
      cols <- intersect(
        active_meta() %>% filter(Tissue == input$tissue_filter) %>% pull(sample_label),
        names(df)
      )
      result <- df %>% select(gene, transcript_id, all_of(cols))
    }
    result <- join_annotations(result)
    if (HAS_TX_ANNOT &&
        !is.null(input$gene_biotype_filter) &&
        input$gene_biotype_filter != "All")
      result <- result %>% filter(Biotype == input$gene_biotype_filter)
    result
  })
  
  output$isoform_table <- renderDT({
    result <- gene_table_data()
    skip_cols <- c("gene",
                   "transcript_id",
                   "RefSeq",
                   "Biotype",
                   "UniProt",
                   "Pfam",
                   "InterPro")
    int_cols <- intersect(c("Exons", "Transcript Length (bp)", "CDS (bp)"),
                          names(result))
    round_cols <- setdiff(names(result), c(skip_cols, int_cols))
    datatable(result,
              options = list(scrollX = TRUE, pageLength = 10),
              rownames = FALSE) %>%
      formatRound(columns = round_cols, digits = 2) %>%
      formatRound(columns = int_cols, digits = 0)
  })
  output$dl_table_gene <- downloadHandler(
    filename = function()
      paste0(
        "isoform_table_",
        input$gene,
        "_",
        input$tissue_filter,
        if (HAS_TX_ANNOT &&
            !is.null(input$gene_biotype_filter) &&
            input$gene_biotype_filter != "All")
          paste0("_", input$gene_biotype_filter)
        else
          "",
        ".csv"
      ),
    content = function(file)
      write.csv(gene_table_data(), file, row.names = FALSE)
  )
  
  # === Isoform Switching ===
  
  switch_label1 <- reactive({
    if (input$switch_mode == "tissue")
      input$switch_tissue1
    else
      input$switch_cl1
  })
  switch_label2 <- reactive({
    if (input$switch_mode == "tissue")
      input$switch_tissue2
    else
      input$switch_cl2
  })
  
  output$switch_cl_tissues <- renderUI({
    req(input$switch_mode == "cellline")
    t1 <- active_meta()$Tissue[match(input$switch_cl1, active_meta()$CellLine)]
    t2 <- active_meta()$Tissue[match(input$switch_cl2, active_meta()$CellLine)]
    tags$div(
      style = "background: #f0f4ff; padding: 6px 10px; border-radius: 5px; font-size: 12px; margin-bottom: 8px;",
      tags$span(style = "font-weight: 600;", input$switch_cl1),
      tags$span(style = "color: #666;", paste0(" (", t1, ")")),
      tags$span(style = "margin: 0 4px;", "vs"),
      tags$span(style = "font-weight: 600;", input$switch_cl2),
      tags$span(style = "color: #666;", paste0(" (", t2, ")"))
    )
  })
  
  switch_results <- eventReactive(input$run_switch, {
    ensure_db()
    dataTableProxy("switch_table") %>% selectRows(NULL)
    if (input$switch_mode == "tissue") {
      t1 <- input$switch_tissue1
      t2 <- input$switch_tissue2
      d1 <- dominant_isoforms_df %>% filter(tissue == t1) %>%
        select(gene,
               dominant_tx,
               dominant_prop,
               total_tpm,
               ci_lower,
               ci_upper,
               n_samples) %>%
        rename(
          dominant_tx1 = dominant_tx,
          prop1 = dominant_prop,
          tpm1 = total_tpm,
          ci_lo1 = ci_lower,
          ci_hi1 = ci_upper,
          n1 = n_samples
        )
      d2 <- dominant_isoforms_df %>% filter(tissue == t2) %>%
        select(gene,
               dominant_tx,
               dominant_prop,
               total_tpm,
               ci_lower,
               ci_upper,
               n_samples) %>%
        rename(
          dominant_tx2 = dominant_tx,
          prop2 = dominant_prop,
          tpm2 = total_tpm,
          ci_lo2 = ci_lower,
          ci_hi2 = ci_upper,
          n2 = n_samples
        )
      n_small <- min(c(d1$n1, d2$n2, Inf), na.rm = TRUE)
      if (is.finite(n_small) && n_small < 3)
        showNotification(
          "One of the selected groups has fewer than 3 samples, so every switch in this comparison is Low confidence.",
          type = "warning", duration = 8)
      av <- tissue_avg_df[, c("transcript_id", t1, t2)]
      names(av) <- c("transcript_id", "avg_g1", "avg_g2")
      result <- inner_join(d1, d2, by = "gene") %>%
        filter(dominant_tx1 != dominant_tx2,
               expr_ok(tpm1, tpm2, input$min_expr, input$switch_expr_rule)) %>%
        left_join(av %>% select(transcript_id, tx2_avg_g1 = avg_g1),
                  by = c("dominant_tx2" = "transcript_id")) %>%
        left_join(av %>% select(transcript_id, tx1_avg_g2 = avg_g2),
                  by = c("dominant_tx1" = "transcript_id")) %>%
        mutate(
          switch_score = abs(prop1 - prop2),
          other_prop1 = ifelse(tpm1 > 0, coalesce(tx2_avg_g1, 0) / tpm1, 0),
          other_prop2 = ifelse(tpm2 > 0, coalesce(tx1_avg_g2, 0) / tpm2, 0),
          confidence = case_when(
            pmin(n1, n2) < 3 ~ "Low",
            is.na(ci_lo1) | is.na(ci_lo2) ~ "Low",
            ci_lo1 > other_prop1 & ci_lo2 > other_prop2 ~ "High",
            TRUE ~ "Medium"
          ),
          dominant_iso_group1 = paste0(gene, " | ", dominant_tx1),
          dominant_iso_group2 = paste0(gene, " | ", dominant_tx2)
        ) %>%
        filter(switch_score >= input$min_switch_score) %>%
        arrange(desc(switch_score))
      if (input$switch_confidence == "high")
        result <- result %>% filter(confidence == "High")
      else if (input$switch_confidence == "high_med")
        result <- result %>% filter(confidence %in% c("High", "Medium"))
      result
    } else {
      cl1_label <- active_meta() %>% filter(CellLine == input$switch_cl1) %>%
        pull(sample_label) %>% head(1)
      cl2_label <- active_meta() %>% filter(CellLine == input$switch_cl2) %>%
        pull(sample_label) %>% head(1)
      if (length(cl1_label) == 0 || length(cl2_label) == 0)
        return(NULL)
      df <- fetch_quant_cols(c(cl1_label, cl2_label))
      names(df)[3:4] <- c("tpm_1", "tpm_2")
      df %>% group_by(gene) %>%
        summarise(
          dominant_tx1 = transcript_id[which.max(tpm_1)],
          prop1 = ifelse(sum(tpm_1) > 0, max(tpm_1) / sum(tpm_1), 0),
          tpm1 = sum(tpm_1),
          dominant_tx2 = transcript_id[which.max(tpm_2)],
          prop2 = ifelse(sum(tpm_2) > 0, max(tpm_2) / sum(tpm_2), 0),
          tpm2 = sum(tpm_2),
          .groups = "drop"
        ) %>%
        filter(dominant_tx1 != dominant_tx2,
               expr_ok(tpm1, tpm2, input$min_expr, input$switch_expr_rule)) %>%
        mutate(
          switch_score = abs(prop1 - prop2),
          confidence = "N/A",
          n1 = 1,
          n2 = 1,
          dominant_iso_group1 = paste0(gene, " | ", dominant_tx1),
          dominant_iso_group2 = paste0(gene, " | ", dominant_tx2)
        ) %>%
        filter(switch_score >= input$min_switch_score) %>%
        arrange(desc(switch_score))
    }
  })
  
  # Shared filtered dataframe - single source of truth for table and downstream plots
  switch_results_filtered <- reactive({
    df <- switch_results()
    if (is.null(df) || nrow(df) == 0)
      return(df)
    if (!is.null(tx_annot_df)) {
      df <- df %>%
        mutate(
          tx1_clean = identity(dominant_tx1),
          tx2_clean = identity(dominant_tx2)
        ) %>%
        left_join(
          tx_annot_df %>% select(transcript_id, biotype) %>% rename(biotype1 = biotype),
          by = c("tx1_clean" = "transcript_id")
        ) %>%
        left_join(
          tx_annot_df %>% select(transcript_id, biotype) %>% rename(biotype2 = biotype),
          by = c("tx2_clean" = "transcript_id")
        ) %>%
        select(-tx1_clean, -tx2_clean)
    }
    if (HAS_TX_ANNOT &&
        !is.null(input$switch_biotype_filter) &&
        input$switch_biotype_filter != "All") {
      df <- df %>% filter(
        biotype1 == input$switch_biotype_filter |
          biotype2 == input$switch_biotype_filter
      )
    }
    df
  })
  
  output$switch_table <- renderDT({
    df <- switch_results_filtered()
    if (is.null(df) || nrow(df) == 0)
      return(NULL)
    lbl1 <- switch_label1()
    lbl2 <- switch_label2()
    if (input$switch_mode == "tissue") {
      if (!is.null(tx_annot_df)) {
        display_df <- df %>% select(
          gene,
          dominant_iso_group1,
          biotype1,
          prop1,
          n1,
          dominant_iso_group2,
          biotype2,
          prop2,
          n2,
          switch_score,
          confidence
        )
        datatable(
          display_df,
          selection = "single",
          options = list(pageLength = 20, scrollX = TRUE),
          rownames = FALSE,
          colnames = c(
            "Gene",
            paste0("Top (", lbl1, ")"),
            "Biotype 1",
            "Prop 1",
            "n1",
            paste0("Top (", lbl2, ")"),
            "Biotype 2",
            "Prop 2",
            "n2",
            "Score",
            "Confidence"
          )
        ) %>%
          formatRound(columns = c("prop1", "prop2", "switch_score"),
                      digits = 3) %>%
          formatStyle(
            "confidence",
            backgroundColor = styleEqual(
              c("High", "Medium", "Low"),
              c("#dce8f5", "#f5e6ce", "#f0e0e0")
            ),
            fontWeight = styleEqual(
              c("High", "Medium", "Low"),
              c("bold", "normal", "normal")
            )
          )
      } else {
        display_df <- df %>% select(
          gene,
          dominant_iso_group1,
          prop1,
          n1,
          dominant_iso_group2,
          prop2,
          n2,
          switch_score,
          confidence
        )
        datatable(
          display_df,
          selection = "single",
          options = list(pageLength = 20, scrollX = TRUE),
          rownames = FALSE,
          colnames = c(
            "Gene",
            paste0("Top (", lbl1, ")"),
            "Prop 1",
            "n1",
            paste0("Top (", lbl2, ")"),
            "Prop 2",
            "n2",
            "Score",
            "Confidence"
          )
        ) %>%
          formatRound(columns = c("prop1", "prop2", "switch_score"),
                      digits = 3) %>%
          formatStyle(
            "confidence",
            backgroundColor = styleEqual(
              c("High", "Medium", "Low"),
              c("#dce8f5", "#f5e6ce", "#f0e0e0")
            ),
            fontWeight = styleEqual(
              c("High", "Medium", "Low"),
              c("bold", "normal", "normal")
            )
          )
      }
    } else {
      if (!is.null(tx_annot_df)) {
        display_df <- df %>% select(
          gene,
          dominant_iso_group1,
          biotype1,
          prop1,
          dominant_iso_group2,
          biotype2,
          prop2,
          switch_score
        )
        datatable(
          display_df,
          selection = "single",
          options = list(pageLength = 20, scrollX = TRUE),
          rownames = FALSE,
          colnames = c(
            "Gene",
            paste0("Top (", lbl1, ")"),
            "Biotype 1",
            "Prop 1",
            paste0("Top (", lbl2, ")"),
            "Biotype 2",
            "Prop 2",
            "Score"
          )
        ) %>%
          formatRound(columns = c("prop1", "prop2", "switch_score"),
                      digits = 3)
      } else {
        display_df <- df %>% select(gene,
                                    dominant_iso_group1,
                                    prop1,
                                    dominant_iso_group2,
                                    prop2,
                                    switch_score)
        datatable(
          display_df,
          selection = "single",
          options = list(pageLength = 20, scrollX = TRUE),
          rownames = FALSE,
          colnames = c(
            "Gene",
            paste0("Top (", lbl1, ")"),
            "Prop 1",
            paste0("Top (", lbl2, ")"),
            "Prop 2",
            "Score"
          )
        ) %>%
          formatRound(columns = c("prop1", "prop2", "switch_score"),
                      digits = 3)
      }
    }
  })
  
  output$switch_plot <- renderPlotly({
    df <- switch_results_filtered()
    sel <- input$switch_table_rows_selected
    if (is.null(df) || is.null(sel) || length(sel) == 0)
      return(NULL)
    gene_name <- df$gene[sel]
    lbl1 <- switch_label1()
    lbl2 <- switch_label2()
    if (input$switch_mode == "tissue") {
      gdata <- fetch_quants("SELECT * FROM quants WHERE gene = ?", params = list(gene_name))
      if (nrow(gdata) == 0)
        return(NULL)
      gdata$display_label <- paste0(gdata$gene, " | ", gdata$transcript_id)
      cols1 <- intersect(
        active_meta() %>% filter(Tissue == input$switch_tissue1) %>% pull(sample_label),
        names(gdata)
      )
      cols2 <- intersect(
        active_meta() %>% filter(Tissue == input$switch_tissue2) %>% pull(sample_label),
        names(gdata)
      )
      if (length(cols1) == 0 || length(cols2) == 0)
        return(NULL)
      plot_df <- tibble(
        transcript = rep(gdata$display_label, 2),
        group = rep(c(lbl1, lbl2), each = nrow(gdata)),
        avg_tpm = c(
          rowMeans(gdata[, cols1, drop = FALSE], na.rm = TRUE),
          rowMeans(gdata[, cols2, drop = FALSE], na.rm = TRUE)
        )
      )
    } else {
      cl1_label <- active_meta() %>% filter(CellLine == input$switch_cl1) %>%
        pull(sample_label) %>% head(1)
      cl2_label <- active_meta() %>% filter(CellLine == input$switch_cl2) %>%
        pull(sample_label) %>% head(1)
      gdata <- fetch_quant_cols(c(cl1_label, cl2_label), "gene = ?", list(gene_name))
      if (nrow(gdata) == 0)
        return(NULL)
      names(gdata)[3:4] <- c("tpm_1", "tpm_2")
      gdata$display_label <- paste0(gdata$gene, " | ", gdata$transcript_id)
      plot_df <- tibble(
        transcript = rep(gdata$display_label, 2),
        group = rep(c(lbl1, lbl2), each = nrow(gdata)),
        avg_tpm = c(gdata$tpm_1, gdata$tpm_2)
      )
    }
    top5 <- plot_df %>% group_by(transcript) %>%
      summarise(total = sum(avg_tpm), .groups = "drop") %>%
      arrange(desc(total)) %>% head(5) %>% pull(transcript)
    plot_df <- plot_df %>%
      mutate(isoform_group = ifelse(transcript %in% top5, transcript, "Other"))
    groups <- unique(plot_df$isoform_group)
    non_other <- groups[groups != "Other"]
    pal_base <- get_colors(length(non_other))
    pal <- if ("Other" %in% groups)
      setNames(c(pal_base, "#79aaff"), c(non_other, "Other"))
    else
      setNames(pal_base, non_other)
    p <- ggplot(plot_df,
                aes(
                  x = group,
                  y = avg_tpm,
                  fill = isoform_group,
                  text = paste0(
                    "Group: ",
                    group,
                    "\nIsoform: ",
                    isoform_group,
                    "\nTPM: ",
                    round(avg_tpm, 2)
                  )
                )) +
      geom_bar(stat = "identity", position = "fill") +
      labs(
        title = paste("Isoform Switch:", gene_name),
        x = NULL,
        y = "Proportion",
        fill = "Isoform"
      ) +
      theme_minimal() + scale_fill_manual(values = pal)
    ggplotly(p, tooltip = "text") %>%
      layout(legend = list(
        orientation = "h",
        y = -0.2,
        font = list(size = 13)
      ))
  })
  
  switch_exon_data <- reactive({
    if (!HAS_TX_EXONS || is.null(tx_exons_df))
      return(NULL)
    df <- switch_results_filtered()
    sel <- input$switch_table_rows_selected
    if (is.null(df) ||
        is.null(sel) || length(sel) == 0)
      return(NULL)
    gene_name <- df$gene[sel]
    tx1 <- identity(df$dominant_tx1[sel])
    tx2 <- identity(df$dominant_tx2[sel])
    expressed_tx <- identity(active_avg()$transcript_id[active_avg()$gene == gene_name])
    g_exons <- tx_exons_df %>% filter(gene == gene_name, transcript_id %in% expressed_tx)
    if (nrow(g_exons) == 0)
      return(NULL)
    strand <- g_exons$strand[1]
    slots <- g_exons %>% distinct(start, end) %>% arrange(start)
    slots$slot <- 0L
    cs <- 0L
    ce <- -1L
    for (i in seq_len(nrow(slots))) {
      if (slots$start[i] > ce) {
        cs <- cs + 1L
        ce <- slots$end[i]
      }
      else
        ce <- max(ce, slots$end[i])
      slots$slot[i] <- cs
    }
    if (strand == "-") {
      max_s <- max(slots$slot)
      slots$slot <- max_s - slots$slot + 1L
    }
    g_exons <- g_exons %>% left_join(slots, by = c("start", "end"))
    pair <- g_exons %>% filter(transcript_id %in% c(tx1, tx2))
    if (nrow(pair) == 0)
      return(NULL)
    slots_iso1 <- unique(pair$slot[pair$transcript_id == tx1])
    slots_iso2 <- unique(pair$slot[pair$transcript_id == tx2])
    shared_slots <- intersect(slots_iso1, slots_iso2)
    diff_slots <- setdiff(union(slots_iso1, slots_iso2), shared_slots)
    len1 <- if (!is.null(tx_annot_df))
      tx_annot_df$tx_length[match(tx1, tx_annot_df$transcript_id)]
    else
      NA
    len2 <- if (!is.null(tx_annot_df))
      tx_annot_df$tx_length[match(tx2, tx_annot_df$transcript_id)]
    else
      NA
    lbl1 <- if (is.na(len1))
      paste0(tx1, " [", switch_label1(), "]")
    else
      paste0(tx1, " (", len1, " nt) [", switch_label1(), "]")
    lbl2 <- if (is.na(len2))
      paste0(tx2, " [", switch_label2(), "]")
    else
      paste0(tx2, " (", len2, " nt) [", switch_label2(), "]")
    pair <- pair %>%
      mutate(
        status = ifelse(slot %in% shared_slots, "Shared", "Differs"),
        iso_label = ifelse(transcript_id == tx1, lbl1, lbl2),
        exon_len = end - start + 1
      )
    pair$iso_label <- factor(pair$iso_label, levels = c(lbl2, lbl1))
    pair$status <- factor(pair$status, levels = c("Shared", "Differs"))
    list(
      exons = pair,
      n_slots = max(slots$slot),
      gene = gene_name,
      n_diff = length(diff_slots)
    )
  })
  
  output$switch_exon_plot <- renderPlotly({
    res <- switch_exon_data()
    validate(
      need(
        !is.null(res),
        "Click a row in the switching table above to compare exon structures."
      )
    )
    df <- res$exons
    title_txt <- paste0(
      "Exon Structure Comparison: ",
      res$gene,
      " (",
      res$n_diff,
      " differing exon slots of ",
      res$n_slots,
      " total)"
    )
    step <- max(1, floor(res$n_slots / 20))
    brks <- seq(1, res$n_slots, by = step)
    p <- ggplot(df, aes(
      x = slot,
      y = iso_label,
      fill = status,
      text = paste0(
        "Isoform: ",
        transcript_id,
        "<br>Exon ",
        exon_number,
        " of this isoform",
        "<br>Gene slot: ",
        slot,
        " of ",
        res$n_slots,
        "<br>",
        status,
        "<br>Length: ",
        exon_len,
        " bp",
        "<br>",
        chr,
        ":",
        start,
        "-",
        end
      )
    )) +
      geom_tile(
        width = 0.9,
        height = 0.7,
        color = "white",
        linewidth = 0.3
      ) +
      geom_text(
        aes(label = exon_number),
        color = "white",
        size = 3,
        fontface = "bold"
      ) +
      scale_fill_manual(
        values = c(Shared = "#225ea8", Differs = "#c51b8a"),
        drop = FALSE,
        name = NULL
      ) +
      scale_x_continuous(breaks = brks, expand = c(0.01, 0.01)) +
      labs(title = title_txt, x = "Exon slot (5' to 3')", y = NULL) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(size = 13, face = "bold"),
        legend.position = "bottom"
      )
    ggplotly(p, source = "switch_exons", tooltip = "text") %>% event_register("plotly_click") %>% layout(
      margin = list(b = 120),
      legend = list(
        orientation = "h",
        y = -0.45,
        x = 0.5,
        xanchor = "center"
      )
    )
  })
  
  output$dl_switch_table <- downloadHandler(
    filename = function()
      paste0(
        "isoform_switches_",
        switch_label1(),
        "_vs_",
        switch_label2(),
        "_minTPM",
        input$min_expr,
        "_score",
        input$min_switch_score,
        if (input$switch_mode == "tissue")
          paste0("_", input$switch_confidence)
        else
          "",
        if (HAS_TX_ANNOT &&
            !is.null(input$switch_biotype_filter) &&
            input$switch_biotype_filter != "All")
          paste0("_", input$switch_biotype_filter)
        else
          "",
        ".csv"
      ),
    content = function(file) {
      df <- switch_results()
      if (!is.null(df))
        write.csv(df, file, row.names = FALSE)
    }
  )
  
  output$dl_switch_plot_pdf <- downloadHandler(
    filename = function()
      paste0("switch_", switch_label1(), "_vs_", switch_label2(), ".pdf"),
    content = function(file) {
      df <- switch_results_filtered()
      sel <- input$switch_table_rows_selected
      if (is.null(df) || is.null(sel) || length(sel) == 0) {
        showNotification("Select a row in the table first.", type = "warning")
        return()
      }
      gene_name <- df$gene[sel]
      lbl1 <- switch_label1()
      lbl2 <- switch_label2()
      gdata <- fetch_quants("SELECT * FROM quants WHERE gene = ?", params = list(gene_name))
      if (nrow(gdata) == 0)
        return()
      gdata$display_label <- paste0(gdata$gene, " | ", gdata$transcript_id)
      cols1 <- intersect(
        active_meta() %>% filter(Tissue == input$switch_tissue1) %>% pull(sample_label),
        names(gdata)
      )
      cols2 <- intersect(
        active_meta() %>% filter(Tissue == input$switch_tissue2) %>% pull(sample_label),
        names(gdata)
      )
      if (length(cols1) == 0 || length(cols2) == 0)
        return()
      plot_df <- tibble(
        transcript = rep(gdata$display_label, 2),
        group = rep(c(lbl1, lbl2), each = nrow(gdata)),
        avg_tpm = c(
          rowMeans(gdata[, cols1, drop = FALSE], na.rm = TRUE),
          rowMeans(gdata[, cols2, drop = FALSE], na.rm = TRUE)
        )
      )
      top5 <- plot_df %>% group_by(transcript) %>% summarise(total = sum(avg_tpm), .groups = "drop") %>%
        arrange(desc(total)) %>% head(5) %>% pull(transcript)
      plot_df <- plot_df %>% mutate(isoform_group = ifelse(transcript %in% top5, transcript, "Other"))
      groups <- unique(plot_df$isoform_group)
      non_other <- groups[groups != "Other"]
      pal_base <- get_colors(length(non_other))
      pal <- if ("Other" %in% groups)
        setNames(c(pal_base, "#79aaff"), c(non_other, "Other"))
      else
        setNames(pal_base, non_other)
      p <- ggplot(plot_df, aes(x = group, y = avg_tpm, fill = isoform_group)) +
        geom_bar(stat = "identity", position = "fill") +
        labs(
          title = paste("Isoform Switch:", gene_name),
          x = NULL,
          y = "Proportion",
          fill = "Isoform"
        ) +
        theme_minimal() + scale_fill_manual(values = pal) +
        theme(
          plot.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 14),
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 15)
        )
      ggsave(file, p, width = 16, height = 7)
    }
  )
  
  # === QC ===
  
  qc_data <- eventReactive(input$run_qc, {
    ensure_db()
    df <- fetch_union("SELECT * FROM qc_stats") %>%
      left_join(active_meta() %>% select(sample_label, CellLine, Tissue),
                by = "sample_label")
    if (length(input$qc_tissue) > 0)
      df <- df %>% filter(Tissue %in% input$qc_tissue)
    df
  })
  
  output$qc_plot <- renderPlotly({
    df <- qc_data()
    if (is.null(df) || nrow(df) == 0)
      return(NULL)
    metric <- input$qc_metric
    metric_label <- switch(
      metric,
      n_isoforms = "Unique Isoforms",
      n_genes = "Unique Genes",
      iso_per_gene = "Avg Iso/Gene"
    )
    p <- ggplot(df, aes(
      x = reorder(Tissue, .data[[metric]], FUN = median),
      y = .data[[metric]],
      fill = Tissue,
      text = paste0(
        "Sample: ",
        CellLine,
        "\nTissue: ",
        Tissue,
        "\n",
        metric_label,
        ": ",
        round(.data[[metric]], 2)
      )
    )) +
      geom_boxplot(
        outlier.shape = NA,
        show.legend = FALSE,
        alpha = 0.75,
        linewidth = 0.4,
        color = "#333333"
      ) +
      geom_jitter(
        aes(color = Tissue),
        width = 0.25,
        size = 1.5,
        alpha = 0.6,
        show.legend = FALSE
      ) +
      scale_fill_manual(values = TISSUE_COLORS, labels = legend_label) +
      scale_color_manual(values = TISSUE_COLORS, labels = legend_label) + coord_flip() + scale_x_discrete(labels = legend_label) +
      labs(title = metric_label, x = NULL, y = metric_label) +
      theme_minimal(base_size = 13) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 14, face = "bold")
      )
    ggplotly(p, tooltip = "text") %>% layout(showlegend = FALSE)
  })
  
  output$qc_table <- renderDT({
    df <- qc_data()
    if (is.null(df))
      return(NULL)
    datatable(
      df %>% select(CellLine, Tissue, n_isoforms, n_genes, iso_per_gene),
      options = list(scrollX = TRUE, pageLength = 15),
      rownames = FALSE,
      colnames = c("Sample", "Tissue / Group", "Isoforms", "Genes", "Iso/Gene")
    ) %>% formatRound(columns = "iso_per_gene", digits = 2)
  })
  output$dl_qc_table <- downloadHandler(
    filename = function()
      paste0(
        "qc_",
        input$qc_metric,
        if (length(input$qc_tissue) > 0)
          paste0("_", paste(input$qc_tissue, collapse = "-"))
        else
          "_all",
        ".csv"
      ),
    content = function(file) {
      df <- qc_data()
      if (!is.null(df))
        write.csv(df, file, row.names = FALSE)
    }
  )
  
  qc_plot_obj <- reactive({
    df <- qc_data()
    if (is.null(df) || nrow(df) == 0)
      return(NULL)
    metric <- input$qc_metric
    metric_label <- switch(
      metric,
      n_isoforms = "Unique Isoforms",
      n_genes = "Unique Genes",
      iso_per_gene = "Avg Iso/Gene"
    )
    ggplot(df, aes(
      x = reorder(Tissue, .data[[metric]], FUN = median),
      y = .data[[metric]],
      fill = Tissue
    )) +
      geom_boxplot(outlier.shape = NA,
                   show.legend = FALSE,
                   alpha = 0.75) +
      geom_jitter(
        aes(color = Tissue),
        width = 0.25,
        size = 1.5,
        alpha = 0.6,
        show.legend = FALSE
      ) +
      scale_fill_manual(values = TISSUE_COLORS, labels = legend_label) + scale_color_manual(values = TISSUE_COLORS, labels = legend_label) + coord_flip() + scale_x_discrete(labels = legend_label) +
      labs(title = metric_label, x = NULL, y = metric_label) +
      theme_minimal(base_size = 13) + theme(
        panel.grid.major.y = element_blank(),
        axis.text.y = element_text(size = 14, face = "bold")
      )
  })
  output$dl_qc_plot <- downloadHandler(
    filename = function()
      paste0("qc_", input$qc_metric, ".pdf"),
    content = function(file) {
      p <- qc_plot_obj()
      if (!is.null(p))
        ggsave(file, p, width = 16, height = 7)
    }
  )
  
  # === Isoform Diversity ===
  
  diversity_gene_data <- eventReactive(input$run_diversity, {
    ensure_db()
    req(input$div_gene)
    if (input$div_mode == "tissue") {
      gdata <- active_avg() %>% filter(gene == input$div_gene)
      if (nrow(gdata) == 0)
        return(NULL)
      tissue_cols <- intersect(names(gdata), active_tissues())
      tibble(
        label = tissue_cols,
        Tissue = tissue_cols,
        Shannon = sapply(tissue_cols, function(t)
          shannon_entropy(gdata[[t]])),
        n_expressed = sapply(tissue_cols, function(t)
          sum(gdata[[t]] > 0)),
        total_tpm = sapply(tissue_cols, function(t)
          sum(gdata[[t]])),
        dominant_isoform = sapply(tissue_cols, function(t) {
          v <- gdata[[t]]
          if (all(v == 0))
            NA_character_
          else
            gdata$display_label[which.max(v)]
        }),
        dominant_prop = sapply(tissue_cols, function(t) {
          v <- gdata[[t]]
          tot <- sum(v, na.rm = TRUE)
          if (tot == 0)
            0
          else
            max(v) / tot
        })
      ) %>% filter(total_tpm >= input$div_min_tpm)
    } else {
      df <- fetch_quants(
                       "SELECT * FROM quants WHERE gene = ?",
                       params = list(input$div_gene))
      if (nrow(df) == 0)
        return(NULL)
      df$display_label <- paste0(df$gene, " | ", df$transcript_id)
      sample_info <- active_meta() %>% select(sample_label, CellLine, Tissue)
      cols <- intersect(sample_info$sample_label, names(df))
      if (length(cols) == 0)
        return(NULL)
      tpm_mat <- as.matrix(df[, cols, drop = FALSE])
      result <- tibble(
        label = sample_info$CellLine[match(cols, sample_info$sample_label)],
        Tissue = sample_info$Tissue[match(cols, sample_info$sample_label)],
        Shannon = sapply(seq_along(cols), function(i)
          shannon_entropy(tpm_mat[, i])),
        n_expressed = sapply(seq_along(cols), function(i)
          sum(tpm_mat[, i] > 0)),
        total_tpm = sapply(seq_along(cols), function(i)
          sum(tpm_mat[, i])),
        dominant_isoform = sapply(seq_along(cols), function(i) {
          v <- tpm_mat[, i]
          if (all(v == 0))
            NA_character_
          else
            df$display_label[which.max(v)]
        }),
        dominant_prop = sapply(seq_along(cols), function(i) {
          v <- tpm_mat[, i]
          tot <- sum(v, na.rm = TRUE)
          if (tot == 0)
            0
          else
            max(v) / tot
        })
      )
      if (length(input$div_tissue_filter) > 0)
        result <- result %>% filter(Tissue %in% input$div_tissue_filter)
      result %>% filter(total_tpm >= input$div_min_tpm)
    }
  })
  
  output$diversity_plot <- renderPlotly({
    df <- diversity_gene_data()
    if (is.null(df))
      return(NULL)
    is_cellline <- (input$div_mode == "cellline")
    if (!is.null(tx_annot_df)) {
      dom_tx_ids <- identity(sub(".*\\| ", "", df$dominant_isoform))
      df$dom_exons <- tx_annot_df$n_exons[match(dom_tx_ids, tx_annot_df$transcript_id)]
    } else {
      df$dom_exons <- NA_integer_
    }
    n_iso <- length(unique(active_avg()$transcript_id[active_avg()$gene == input$div_gene]))
    ttl <- paste0("Isoform Diversity: ",
                  input$div_gene,
                  " (",
                  n_iso,
                  " isoforms)")
    if (input$div_plot_type == "bubble") {
      p <- ggplot(
        df,
        aes(
          x = Shannon,
          y = reorder(label, Shannon),
          size = n_expressed,
          color = if (is_cellline)
            Tissue
          else
            Shannon,
          text = paste0(
            ifelse(is_cellline, "Sample: ", "Tissue: "),
            label,
            ifelse(is_cellline, paste0("\nTissue: ", Tissue), ""),
            "\nEntropy: ",
            round(Shannon, 3),
            "\nIsoforms: ",
            n_expressed,
            "\nDominant: ",
            dominant_isoform,
            " (",
            round(dominant_prop * 100, 1),
            "% | ",
            ifelse(is.na(dom_exons), "?", as.character(dom_exons)),
            " exons)",
            "\nTotal TPM: ",
            round(total_tpm, 2)
          )
        )
      ) +
        geom_point(alpha = 0.85) +
        scale_size_continuous(range = c(3, if (is_cellline)
          8
          else
            14), name = "Isoforms") +
        scale_x_continuous(limits = c(0, 1)) +
        labs(title = ttl, x = "Normalized Entropy", y = NULL) +
        theme_minimal(base_size = 13) + theme(panel.grid.minor = element_blank())
      if (is_cellline) {
        p <- p + scale_color_manual(values = TISSUE_COLORS, labels = legend_label, name = "Tissue")
        if (nrow(df) > 50)
          p <- p + theme(axis.text.y = element_blank())
      } else {
        p <- p + scale_color_gradientn(colours = TROPICAL,
                                       limits = c(0, 1),
                                       name = "Entropy")
      }
    } else {
      if (is_cellline) {
        p <- ggplot(df,
                    aes(
                      x = reorder(label, Shannon),
                      y = Shannon,
                      fill = Tissue,
                      text = paste0(
                        "Sample: ",
                        label,
                        "\nTissue: ",
                        Tissue,
                        "\nEntropy: ",
                        round(Shannon, 3)
                      )
                    )) +
          geom_col() + scale_fill_manual(values = TISSUE_COLORS, labels = legend_label) +
          coord_flip() + scale_x_discrete(labels = legend_label) + scale_y_continuous(limits = c(0, 1)) +
          labs(title = ttl,
               x = NULL,
               y = "Normalized Entropy") +
          theme_minimal(base_size = 13)
        if (nrow(df) > 50)
          p <- p + theme(axis.text.y = element_blank())
      } else {
        p <- ggplot(df,
                    aes(
                      x = reorder(label, Shannon),
                      y = Shannon,
                      fill = Shannon,
                      text = paste0("Tissue: ", label, "\nEntropy: ", round(Shannon, 3))
                    )) +
          geom_col() +
          scale_fill_gradientn(colours = TROPICAL, limits = c(0, 1)) +
          coord_flip() + scale_y_continuous(limits = c(0, 1)) +
          labs(title = ttl,
               x = NULL,
               y = "Normalized Entropy") +
          theme_minimal(base_size = 13)
      }
    }
    ggplotly(p, tooltip = "text")
  })
  
  output$diversity_gene_table <- renderDT({
    df <- diversity_gene_data()
    if (is.null(df))
      return(NULL)
    display <- df %>% arrange(desc(Shannon))
    if (!is.null(tx_annot_df)) {
      dom_tx_ids <- identity(sub(".*\\| ", "", display$dominant_isoform))
      display$dom_exons <- tx_annot_df$n_exons[match(dom_tx_ids, tx_annot_df$transcript_id)]
      gene_tx_ids <- identity(active_avg()$transcript_id[active_avg()$gene == input$div_gene])
      gene_ec <- tx_annot_df$n_exons[match(gene_tx_ids, tx_annot_df$transcript_id)]
      gene_ec <- gene_ec[!is.na(gene_ec)]
      display$gene_exon_range <- if (length(gene_ec) > 0) {
        if (min(gene_ec) == max(gene_ec))
          as.character(min(gene_ec))
        else
          paste0(min(gene_ec), "-", max(gene_ec))
      } else
        NA_character_
    } else {
      display$dom_exons <- NA_integer_
      display$gene_exon_range <- NA_character_
    }
    if (input$div_mode == "cellline") {
      datatable(
        display %>% select(
          label,
          Tissue,
          Shannon,
          n_expressed,
          dominant_isoform,
          dom_exons,
          gene_exon_range,
          dominant_prop,
          total_tpm
        ),
        selection = "single",
        options = list(pageLength = 20, scrollX = TRUE),
        rownames = FALSE,
        colnames = c(
          "Sample",
          "Tissue / Group",
          "Entropy",
          "Isoforms",
          "Dominant",
          "Dom. Exons",
          "Gene Exon Range",
          "Dom. Prop.",
          "Total TPM"
        )
      ) %>%
        formatRound(
          columns = c("Shannon", "dominant_prop", "total_tpm"),
          digits = 3
        )
    } else {
      datatable(
        display %>% select(
          label,
          Shannon,
          n_expressed,
          dominant_isoform,
          dom_exons,
          gene_exon_range,
          dominant_prop,
          total_tpm
        ),
        selection = "single",
        options = list(
          pageLength = 26,
          scrollX = TRUE,
          dom = 't'
        ),
        rownames = FALSE,
        colnames = c(
          "Tissue / Group",
          "Entropy",
          "Isoforms",
          "Dominant",
          "Dom. Exons",
          "Gene Exon Range",
          "Dom. Prop.",
          "Total TPM"
        )
      ) %>%
        formatRound(
          columns = c("Shannon", "dominant_prop", "total_tpm"),
          digits = 3
        )
    }
  })
  
  output$diversity_range_table <- renderDT({
    iso_counts <- active_avg() %>% group_by(gene) %>%
      summarise(n_isoforms = n(), .groups = "drop")
    display_df <- diversity_range_df %>%
      left_join(iso_counts, by = "gene") %>%
      filter(n_isoforms >= 3) %>%
      select(
        gene,
        n_isoforms,
        entropy_range,
        most_diverse_tissue,
        max_entropy,
        least_diverse_tissue,
        min_entropy
      ) %>%
      arrange(desc(entropy_range))
    datatable(
      display_df,
      options = list(pageLength = 25, scrollX = TRUE),
      rownames = FALSE,
      colnames = c(
        "Gene",
        "Isoforms",
        "Range",
        "Most Diverse",
        "Max",
        "Least Diverse",
        "Min"
      ),
      caption = "Genes with 3+ isoforms ranked by normalized entropy range across tissues."
    ) %>%
      formatRound(
        columns = c("entropy_range", "max_entropy", "min_entropy"),
        digits = 3
      )
  })
  
  output$dl_diversity <- downloadHandler(
    filename = function()
      paste0(
        "diversity_",
        input$div_gene,
        "_minTPM",
        input$div_min_tpm,
        ".csv"
      ),
    content = function(file) {
      df <- diversity_gene_data()
      if (!is.null(df))
        write.csv(df, file, row.names = FALSE)
    }
  )
  
  output$dl_diversity_plot <- downloadHandler(
    filename = function()
      paste0(
        "diversity_",
        input$div_gene,
        "_minTPM",
        input$div_min_tpm,
        ".pdf"
      ),
    content = function(file) {
      df <- diversity_gene_data()
      if (is.null(df))
        return()
      if (input$div_mode == "cellline" && nrow(df) > 30) {
        showNotification(
          paste0(
            "PDF not available for ",
            nrow(df),
            " samples (limit: 30). Filter by tissue or group to reduce, or use the interactive plot."
          ),
          type = "warning"
        )
        return()
      }
      n_iso <- length(unique(active_avg()$transcript_id[active_avg()$gene == input$div_gene]))
      ttl <- paste0("Isoform Diversity: ",
                    input$div_gene,
                    " (",
                    n_iso,
                    " isoforms)")
      is_cl_pdf <- (input$div_mode == "cellline")
      if (input$div_plot_type == "bubble") {
        p <- ggplot(df,
                    aes(
                      x = Shannon,
                      y = reorder(label, Shannon),
                      size = n_expressed,
                      color = if (is_cl_pdf)
                        Tissue
                      else
                        Shannon
                    )) +
          geom_point(alpha = 0.85) + scale_size_continuous(range = c(3, if (is_cl_pdf)
            8
            else
              14), name = "Isoforms") +
          scale_x_continuous(limits = c(0, 1))
        if (is_cl_pdf) {
          p <- p + scale_color_manual(values = TISSUE_COLORS, labels = legend_label, name = "Tissue")
        } else {
          p <- p + scale_color_gradientn(colours = TROPICAL,
                                         limits = c(0, 1),
                                         name = "Entropy")
        }
        p <- p + labs(title = ttl, x = "Normalized Entropy", y = NULL) + theme_minimal(base_size = 13)
      } else {
        p <- ggplot(df, aes(
          x = reorder(label, Shannon),
          y = Shannon,
          fill = Shannon
        )) +
          geom_col() + scale_fill_gradientn(colours = TROPICAL, limits = c(0, 1)) +
          coord_flip() + scale_y_continuous(limits = c(0, 1)) +
          labs(title = ttl, x = NULL, y = "Normalized Entropy") + theme_minimal(base_size = 13)
      }
      ggsave(file, p, width = 16, height = 7)
    }
  )
  
  output$dl_diversity_ranking <- downloadHandler(
    filename = function()
      "diversity_global_ranking.csv",
    content = function(file) {
      iso_counts <- active_avg() %>% group_by(gene) %>% summarise(n_isoforms = n(), .groups = "drop")
      df <- diversity_range_df %>% left_join(iso_counts, by = "gene") %>% filter(n_isoforms >= 3) %>% arrange(desc(entropy_range))
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  exon_structure_data <- reactive({
    req(input$div_gene)
    if (!HAS_TX_EXONS || is.null(tx_exons_df))
      return(NULL)
    expressed_tx <- identity(active_avg()$transcript_id[active_avg()$gene == input$div_gene])
    if (length(expressed_tx) == 0)
      return(NULL)
    g_exons <- tx_exons_df %>%
      filter(gene == input$div_gene, transcript_id %in% expressed_tx)
    if (nrow(g_exons) == 0)
      return(NULL)
    strand <- g_exons$strand[1]
    slots <- g_exons %>% distinct(start, end) %>% arrange(start)
    slots$slot <- 0L
    cur_slot <- 0L
    cur_end <- -1L
    for (i in seq_len(nrow(slots))) {
      if (slots$start[i] > cur_end) {
        cur_slot <- cur_slot + 1L
        cur_end <- slots$end[i]
      } else {
        cur_end <- max(cur_end, slots$end[i])
      }
      slots$slot[i] <- cur_slot
    }
    if (strand == "-") {
      max_s <- max(slots$slot)
      slots$slot <- max_s - slots$slot + 1L
    }
    g_exons <- g_exons %>% left_join(slots, by = c("start", "end"))
    iso_counts <- g_exons %>%
      group_by(transcript_id) %>%
      summarise(n_ex = n(), .groups = "drop") %>%
      arrange(desc(n_ex))
    if (HAS_REFSEQ && !is.null(refseq_map_df)) {
      iso_counts$nm_id <- refseq_map_df$refseq[match(iso_counts$transcript_id,
                                                     refseq_map_df$ensembl_transcript_id)]
    } else {
      iso_counts$nm_id <- NA_character_
    }
    iso_counts$id_str <- ifelse(
      is.na(iso_counts$nm_id) | iso_counts$nm_id == "",
      iso_counts$transcript_id,
      paste0(iso_counts$transcript_id, " / ", iso_counts$nm_id)
    )
    if (!is.null(tx_annot_df)) {
      iso_counts$tx_length <- tx_annot_df$tx_length[match(iso_counts$transcript_id, tx_annot_df$transcript_id)]
      iso_counts$iso_label <- ifelse(
        is.na(iso_counts$tx_length),
        paste0(iso_counts$id_str, " (", iso_counts$n_ex, " ex)"),
        paste0(
          iso_counts$id_str,
          " (",
          iso_counts$n_ex,
          " ex, ",
          iso_counts$tx_length,
          " nt)"
        )
      )
    } else {
      iso_counts$iso_label <- paste0(iso_counts$id_str, " (", iso_counts$n_ex, " ex)")
    }
    g_exons <- g_exons %>%
      left_join(iso_counts %>% select(transcript_id, iso_label), by = "transcript_id") %>%
      mutate(iso_label = factor(iso_label, levels = rev(iso_counts$iso_label)))
    list(
      exons = g_exons,
      n_slots = max(slots$slot),
      n_iso = nrow(iso_counts),
      strand = strand
    )
  })
  
  output$diversity_exon_plot <- renderPlotly({
    res <- exon_structure_data()
    validate(
      need(
        !is.null(res),
        "No exon structure data available for this gene. Check that transcript_exons.csv.gz is present in shiny_data/."
      )
    )
    highlight_tx <- ""
    sel <- input$diversity_gene_table_rows_selected
    if (!is.null(sel) && length(sel) > 0) {
      dd <- tryCatch(
        diversity_gene_data(),
        error = function(e)
          NULL
      )
      if (!is.null(dd) && nrow(dd) >= sel[1]) {
        display <- dd %>% arrange(desc(Shannon))
        dom_label <- as.character(display$dominant_isoform[sel[1]])
        if (!is.na(dom_label) && nchar(dom_label) > 0) {
          highlight_tx <- identity(sub(".*\\| ", "", dom_label))
        }
      }
    }
    df <- res$exons
    df$tx_clean <- identity(df$transcript_id)
    df$exon_len <- df$end - df$start + 1
    df$fill_group <- ifelse(df$tx_clean == highlight_tx &
                              nchar(highlight_tx) > 0,
                            "highlight",
                            "normal")
    df$fill_group <- factor(df$fill_group, levels = c("normal", "highlight"))
    title_txt <- paste0(
      "Exon Structure: ",
      input$div_gene,
      " (",
      res$n_iso,
      " isoforms, ",
      res$n_slots,
      " unique exons, strand ",
      res$strand,
      ")"
    )
    if (nchar(highlight_tx) > 0)
      title_txt <- paste0(title_txt, "  |  ", highlight_tx)
    step <- max(1, floor(res$n_slots / 20))
    brks <- seq(1, res$n_slots, by = step)
    p <- ggplot(df,
                aes(
                  x = slot,
                  y = iso_label,
                  fill = fill_group,
                  text = paste0(
                    "Isoform: ",
                    transcript_id,
                    "<br>Gene slot: ",
                    slot,
                    " of ",
                    res$n_slots,
                    "<br>Exon in this isoform: #",
                    exon_number,
                    "<br>Length: ",
                    exon_len,
                    " bp",
                    "<br>Genomic: ",
                    chr,
                    ":",
                    start,
                    "-",
                    end
                  )
                )) +
      geom_tile(
        width = 0.9,
        height = 0.8,
        color = "white",
        linewidth = 0.2
      ) +
      geom_text(
        aes(label = exon_number),
        color = "white",
        size = 2.6,
        fontface = "bold"
      ) +
      scale_fill_manual(
        values = c(normal = "#225ea8", highlight = "#c51b8a"),
        drop = FALSE,
        guide = "none"
      ) +
      scale_x_continuous(breaks = brks, expand = c(0.01, 0.01)) +
      labs(title = title_txt, x = "Exon slot (5' to 3')", y = NULL) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 13, face = "bold")
      )
    ggplotly(p, source = "diversity_exons", tooltip = "text") %>%
      event_register("plotly_click") %>%
      layout(showlegend = FALSE)
  })
  
  observe({
    res <- exon_structure_data()
    if (is.null(res))
      return()
    tx_choices <- sort(unique(as.character(res$exons$transcript_id)))
    updateSelectizeInput(
      session,
      "div_zoom_tx",
      choices = tx_choices,
      selected = if (length(tx_choices) > 0)
        tx_choices[1]
      else
        "",
      server = TRUE
    )
  })
  
  selected_zoom_tx <- reactive({
    click <- event_data("plotly_click", source = "diversity_exons")
    res <- exon_structure_data()
    if (!is.null(click) && !is.null(click$y) && !is.null(res)) {
      lvls <- levels(res$exons$iso_label)
      y_val <- click$y
      tx_label <- NULL
      if (is.numeric(y_val)) {
        idx <- as.integer(round(y_val))
        if (idx >= 1 && idx <= length(lvls))
          tx_label <- lvls[idx]
      } else {
        s_y <- as.character(y_val)
        if (s_y %in% lvls) {
          tx_label <- s_y
        } else if (suppressWarnings(!is.na(as.integer(s_y)))) {
          idx <- as.integer(s_y)
          if (idx >= 1 && idx <= length(lvls))
            tx_label <- lvls[idx]
        }
      }
      if (!is.null(tx_label)) {
        tx <- sub(" .*$", "", tx_label)
        return(trimws(tx))
      }
    }
    if (!is.null(input$div_zoom_tx) &&
        nchar(input$div_zoom_tx) > 0) {
      return(input$div_zoom_tx)
    }
    NULL
  })
  
  output$diversity_zoom_plot <- renderPlotly({
    tx <- selected_zoom_tx()
    validate(need(
      !is.null(tx),
      "Click an isoform in the grid above or pick one from the dropdown."
    ))
    if (!HAS_TX_EXONS || is.null(tx_exons_df)) {
      return(
        plotly_empty(type = "scatter", mode = "markers") %>%
          layout(title = "No exon data available")
      )
    }
    ex <- tx_exons_df %>%
      filter(transcript_id == tx) %>%
      arrange(start)
    validate(need(nrow(ex) > 0, paste0("No exon data for ", tx)))
    strand <- ex$strand[1]
    chr <- ex$chr[1]
    gene_min <- min(ex$start)
    gene_max <- max(ex$end)
    span_bp <- gene_max - gene_min
    ex <- ex %>%
      mutate(
        exon_len = end - start + 1,
        next_start = lead(start),
        intron_to_next = next_start - end - 1,
        label = as.character(exon_number),
        tip = paste0(
          "Exon ",
          exon_number,
          " of ",
          nrow(ex),
          "<br>Length: ",
          format(exon_len, big.mark = ","),
          " bp",
          "<br>",
          chr,
          ":",
          format(start, big.mark = ","),
          "-",
          format(end, big.mark = ","),
          ifelse(
            is.na(intron_to_next),
            "",
            paste0(
              "<br>Intron to next: ",
              format(intron_to_next, big.mark = ","),
              " bp"
            )
          )
        )
      )
    title_txt <- paste0(
      "Single Transcript: ",
      tx,
      "  (",
      nrow(ex),
      " exons, ",
      format(span_bp, big.mark = ","),
      " bp span on ",
      chr,
      ", strand ",
      strand,
      ")"
    )
    p <- ggplot(ex) +
      geom_segment(
        aes(
          x = gene_min,
          xend = gene_max,
          y = 0,
          yend = 0
        ),
        color = "#888888",
        linewidth = 0.6
      ) +
      geom_rect(
        aes(
          xmin = start,
          xmax = end,
          ymin = -0.45,
          ymax = 0.45,
          text = tip
        ),
        fill = "#225ea8",
        color = "white",
        linewidth = 0.3
      ) +
      geom_text(
        aes(
          x = (start + end) / 2,
          y = 0,
          label = label
        ),
        color = "white",
        size = 3,
        fontface = "bold"
      ) +
      scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
      scale_x_continuous(
        labels = function(x)
          format(x, big.mark = ",", scientific = FALSE),
        expand = expansion(mult = 0.01)
      ) +
      labs(
        title = title_txt,
        x = paste0(chr, " genomic position"),
        y = NULL
      ) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, face = "bold")
      )
    ggplotly(p, tooltip = "text") %>% layout(showlegend = FALSE)
  })
  
  
  output$dl_diversity_exon_grid <- downloadHandler(
    filename = function()
      paste0("exon_structure_", input$div_gene, ".pdf"),
    content = function(file) {
      res <- exon_structure_data()
      if (is.null(res))
        return()
      p <- ggplot(res$exons, aes(x = slot, y = iso_label)) +
        geom_tile(
          aes(fill = "present"),
          color = "white",
          linewidth = 0.3,
          width = 0.9,
          height = 0.85
        ) +
        geom_text(
          aes(label = exon_number),
          color = "white",
          size = 2.6,
          fontface = "bold"
        ) +
        scale_fill_manual(values = c("present" = "#225ea8")) +
        scale_x_continuous(
          breaks = seq(1, res$n_slots, by = max(1, floor(res$n_slots / 20))),
          limits = c(0.5, res$n_slots + 0.5),
          expand = c(0, 0)
        ) +
        labs(
          title = paste0(
            "Exon Structure: ",
            input$div_gene,
            " (",
            length(unique(res$exons$iso_label)),
            " isoforms, ",
            res$n_slots,
            " unique exons, strand ",
            res$strand,
            ")"
          ),
          x = "Exon slot (5' to 3')",
          y = NULL
        ) +
        theme_minimal(base_size = 11) +
        theme(
          legend.position = "none",
          panel.grid = element_blank(),
          plot.title = element_text(
            face = "bold",
            size = 13,
            hjust = 0.5
          )
        )
      ggsave(
        file,
        p,
        width = 14,
        height = max(5, 0.35 * length(unique(
          res$exons$iso_label
        )) + 2),
        limitsize = FALSE
      )
    }
  )
  
  output$dl_diversity_exon_zoom <- downloadHandler(
    filename = function() {
      tx <- selected_zoom_tx()
      paste0("exon_detail_", if (is.null(tx))
        "transcript"
        else
          tx, ".pdf")
    },
    content = function(file) {
      tx <- selected_zoom_tx()
      if (is.null(tx) ||
          !HAS_TX_EXONS || is.null(tx_exons_df))
        return()
      ex <- tx_exons_df %>% filter(transcript_id == tx) %>% arrange(start)
      if (nrow(ex) == 0)
        return()
      strand <- ex$strand[1]
      chr <- ex$chr[1]
      gene_min <- min(ex$start)
      gene_max <- max(ex$end)
      span_bp <- gene_max - gene_min
      p <- ggplot(ex) +
        geom_segment(
          aes(
            x = gene_min,
            xend = gene_max,
            y = 0,
            yend = 0
          ),
          color = "#888888",
          linewidth = 0.6
        ) +
        geom_rect(
          aes(
            xmin = start,
            xmax = end,
            ymin = -0.45,
            ymax = 0.45
          ),
          fill = "#225ea8",
          color = "white",
          linewidth = 0.3
        ) +
        geom_text(
          aes(
            x = (start + end) / 2,
            y = 0,
            label = exon_number
          ),
          color = "white",
          size = 3,
          fontface = "bold"
        ) +
        scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
        scale_x_continuous(
          labels = function(x)
            format(x, big.mark = ",", scientific = FALSE),
          expand = expansion(mult = 0.01)
        ) +
        labs(
          title = paste0(
            "Single Transcript: ",
            tx,
            "  (",
            nrow(ex),
            " exons, ",
            format(span_bp, big.mark = ","),
            " bp span on ",
            chr,
            ", strand ",
            strand,
            ")"
          ),
          x = paste0(chr, " genomic position"),
          y = NULL
        ) +
        theme_minimal(base_size = 11) +
        theme(
          axis.text.y = element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          plot.title = element_text(size = 12, face = "bold")
        )
      ggsave(file, p, width = 14, height = 3.5)
    }
  )
  
  # === Tissue-Specific ===
  
  specific_results <- eventReactive(input$run_specific, {
    ensure_db()
    dataTableProxy("specific_table") %>% selectRows(NULL)
    tissue_cols <- intersect(names(active_avg()), active_tissues())
    base_df <- active_avg()
    if (!is.null(input$spec_gene) && nchar(input$spec_gene) > 0) {
      base_df <- base_df %>% filter(gene == input$spec_gene)
      if (nrow(base_df) == 0) {
        showNotification(paste0("Gene '", input$spec_gene, "' not found."), type = "warning")
        return(NULL)
      }
    }
    base_df %>%
      mutate(
        target_tpm = .data[[input$spec_tissue]],
        total_tpm = rowSums(across(all_of(tissue_cols)), na.rm = TRUE),
        specificity = ifelse(total_tpm > 0, target_tpm / total_tpm, 0)
      ) %>%
      filter(
        total_tpm > 0,
        specificity >= input$spec_threshold,
        target_tpm >= input$spec_min_tpm
      ) %>%
      arrange(desc(specificity)) %>%
      select(gene,
             transcript_id,
             display_label,
             target_tpm,
             total_tpm,
             specificity)
  })
  
  output$specific_table <- renderDT({
    df <- specific_results()
    if (is.null(df) || nrow(df) == 0) {
      if (!is.null(input$spec_gene) && nchar(input$spec_gene) > 0) {
        tcols <- intersect(names(active_avg()), active_tissues())
        diag <- active_avg() %>%
          filter(gene == input$spec_gene) %>%
          mutate(
            target_tpm = .data[[input$spec_tissue]],
            total_tpm = rowSums(across(all_of(tcols)), na.rm = TRUE),
            frac_in_target = ifelse(total_tpm > 0, target_tpm / total_tpm, 0)
          )
        if (nrow(diag) == 0) {
          return(datatable(
            data.frame(
              Result = paste0("Gene '", input$spec_gene, "' not found in dataset.")
            ),
            options = list(dom = 't'),
            rownames = FALSE
          ))
        }
        top_tissue <- sapply(seq_len(nrow(diag)), function(i) {
          vals <- as.numeric(diag[i, tcols])
          if (all(vals == 0, na.rm = TRUE))
            NA_character_
          else
            tcols[which.max(vals)]
        })
        top_frac <- sapply(seq_len(nrow(diag)), function(i) {
          vals <- as.numeric(diag[i, tcols])
          tot <- sum(vals, na.rm = TRUE)
          if (tot == 0)
            0
          else
            max(vals, na.rm = TRUE) / tot
        })
        diag_out <- data.frame(
          Isoform = diag$transcript_id,
          `Avg TPM in target` = round(diag$target_tpm, 2),
          `Fraction in target` = round(diag$frac_in_target, 3),
          `Dominant tissue` = top_tissue,
          `Fraction in dominant` = round(top_frac, 3),
          check.names = FALSE,
          stringsAsFactors = FALSE
        )
        diag_out <- diag_out[order(-diag_out$`Avg TPM in target`), ]
        header <- paste0(
          input$spec_gene,
          " does not meet cutoffs in ",
          input$spec_tissue,
          " (fraction >= ",
          input$spec_threshold,
          ", avg TPM >= ",
          input$spec_min_tpm,
          "). Per-isoform diagnostic below."
        )
        return(datatable(
          diag_out,
          options = list(
            pageLength = 10,
            scrollX = TRUE,
            dom = 't'
          ),
          rownames = FALSE,
          caption = header
        ))
      } else {
        msg <- paste0(
          "No tissue-specific isoforms found in ",
          input$spec_tissue,
          " at fraction >= ",
          input$spec_threshold,
          ", avg TPM >= ",
          input$spec_min_tpm,
          ". Try lowering either slider."
        )
        return(datatable(
          data.frame(Result = msg),
          options = list(dom = 't'),
          rownames = FALSE
        ))
      }
    }
    display_df <- df %>%
      select(gene,
             display_label,
             transcript_id,
             target_tpm,
             specificity) %>%
      rename(
        Gene = gene,
        Isoform = display_label,
        `Avg TPM` = target_tpm,
        `Tissue Fraction` = specificity
      )
    display_df <- join_annotations(display_df)
    if (HAS_TX_ANNOT &&
        !is.null(input$spec_biotype_filter) &&
        input$spec_biotype_filter != "All")
      display_df <- display_df %>% filter(Biotype == input$spec_biotype_filter)
    display_df <- display_df %>% select(-transcript_id)
    skip_cols <- c("Gene", "Isoform", "Biotype", "UniProt", "Pfam", "InterPro")
    int_cols <- intersect(c("Exons", "Transcript Length (bp)", "CDS (bp)"),
                          names(display_df))
    round_cols <- setdiff(names(display_df), c(skip_cols, int_cols))
    datatable(
      display_df,
      selection = "single",
      options = list(pageLength = 20, scrollX = TRUE),
      rownames = FALSE
    ) %>%
      formatRound(columns = round_cols, digits = 3) %>%
      formatRound(columns = int_cols, digits = 0)
  })
  
  output$specific_plot <- renderPlotly({
    df <- specific_results()
    sel <- input$specific_table_rows_selected
    if (is.null(df) || is.null(sel) || length(sel) == 0)
      return(NULL)
    tx_id <- df$transcript_id[sel]
    tx_label <- df$display_label[sel]
    if (input$spec_detail_mode == "tissue") {
      tissue_cols <- intersect(names(active_avg()), active_tissues())
      row_data <- active_avg() %>% filter(transcript_id == tx_id)
      if (nrow(row_data) == 0)
        return(NULL)
      tpm_vals <- as.numeric(unlist(row_data[1, tissue_cols]))
      plot_df <- tibble(Tissue = tissue_cols, TPM = tpm_vals)
      total <- sum(plot_df$TPM, na.rm = TRUE)
      p <- ggplot(plot_df,
                  aes(
                    x = reorder(Tissue, TPM),
                    y = TPM,
                    fill = Tissue,
                    text = paste0(
                      "Tissue: ",
                      Tissue,
                      "\nAvg TPM: ",
                      round(TPM, 3),
                      "\nFraction: ",
                      round(if (total > 0)
                        TPM / total
                        else
                          0, 3)
                    )
                  )) +
        geom_col() + scale_fill_manual(values = TISSUE_COLORS, labels = legend_label) + coord_flip() + scale_x_discrete(labels = legend_label) +
        labs(
          title = paste("Expression:", tx_label),
          x = NULL,
          y = "Avg TPM"
        ) +
        theme_minimal(base_size = 13) + theme(legend.position = "none")
    } else {
      row_df <- fetch_quants("SELECT * FROM quants WHERE transcript_id = ?",
                           params = list(tx_id))
      if (nrow(row_df) == 0)
        return(NULL)
      sample_info <- active_meta() %>% select(sample_label, CellLine, Tissue)
      cols <- intersect(sample_info$sample_label, names(row_df))
      vals <- as.numeric(row_df[1, cols])
      plot_df <- tibble(
        CellLine = sample_info$CellLine[match(cols, sample_info$sample_label)],
        Tissue = sample_info$Tissue[match(cols, sample_info$sample_label)],
        TPM = vals
      )
      p <- ggplot(plot_df,
                  aes(
                    x = reorder(Tissue, TPM, FUN = median),
                    y = TPM,
                    fill = Tissue,
                    text = paste0(
                      "Tissue: ",
                      Tissue,
                      "\nSample: ",
                      CellLine,
                      "\nTPM: ",
                      round(TPM, 3)
                    )
                  )) +
        geom_boxplot(outlier.shape = NA, alpha = 0.7) +
        geom_jitter(width = 0.2,
                    size = 1.2,
                    alpha = 0.5) +
        scale_fill_manual(values = TISSUE_COLORS, labels = legend_label) + coord_flip() + scale_x_discrete(labels = legend_label) +
        labs(
          title = paste("Distribution:", tx_label),
          x = NULL,
          y = "TPM"
        ) +
        theme_minimal(base_size = 13) + theme(legend.position = "none")
    }
    ggplotly(p, tooltip = "text")
  })
  
  output$specific_detail_header <- renderText({
    if (input$spec_detail_mode == "tissue")
      "Tissue Average Expression"
    else
      "Per-Cell-Line Expression"
  })
  
  output$specific_detail_table <- renderDT({
    df <- specific_results()
    sel <- input$specific_table_rows_selected
    if (is.null(df) || is.null(sel) || length(sel) == 0)
      return(NULL)
    tx_id <- df$transcript_id[sel]
    if (input$spec_detail_mode == "tissue") {
      tissue_cols <- intersect(names(active_avg()), active_tissues())
      row_data <- active_avg() %>% filter(transcript_id == tx_id)
      if (nrow(row_data) == 0)
        return(NULL)
      tpm_vals <- as.numeric(unlist(row_data[1, tissue_cols]))
      total <- sum(tpm_vals, na.rm = TRUE)
      detail_df <- tibble(
        Tissue = tissue_cols,
        Avg_TPM = tpm_vals,
        Fraction = if (total > 0)
          tpm_vals / total
        else
          rep(0, length(tpm_vals))
      ) %>% arrange(desc(Avg_TPM))
      datatable(
        detail_df,
        options = list(
          pageLength = 26,
          scrollX = TRUE,
          dom = 't'
        ),
        rownames = FALSE,
        colnames = c("Tissue", "Avg TPM", "Fraction")
      ) %>% formatRound(columns = c("Avg_TPM", "Fraction"),
                        digits = 3)
    } else {
      row_df <- fetch_quants("SELECT * FROM quants WHERE transcript_id = ?",
                           params = list(tx_id))
      if (nrow(row_df) == 0)
        return(NULL)
      sample_info <- active_meta() %>% select(sample_label, CellLine, Tissue)
      cols <- intersect(sample_info$sample_label, names(row_df))
      vals <- as.numeric(row_df[1, cols])
      detail_df <- tibble(
        CellLine = sample_info$CellLine[match(cols, sample_info$sample_label)],
        Tissue = sample_info$Tissue[match(cols, sample_info$sample_label)],
        TPM = vals
      ) %>% arrange(desc(TPM))
      datatable(
        detail_df,
        options = list(pageLength = 20, scrollX = TRUE),
        rownames = FALSE,
        colnames = c("Sample", "Tissue / Group", "TPM")
      ) %>% formatRound(columns = "TPM", digits = 3)
    }
  })
  
  specific_exon_data <- reactive({
    if (!HAS_TX_EXONS || is.null(tx_exons_df))
      return(NULL)
    df <- specific_results()
    sel <- input$specific_table_rows_selected
    if (is.null(df) ||
        is.null(sel) || length(sel) == 0)
      return(NULL)
    gene_name <- df$gene[sel]
    highlight_tx <- identity(df$transcript_id[sel])
    expressed_tx <- identity(active_avg()$transcript_id[active_avg()$gene == gene_name])
    g_exons <- tx_exons_df %>% filter(gene == gene_name, transcript_id %in% expressed_tx)
    if (nrow(g_exons) == 0)
      return(NULL)
    strand <- g_exons$strand[1]
    slots <- g_exons %>% distinct(start, end) %>% arrange(start)
    slots$slot <- 0L
    cs <- 0L
    ce <- -1L
    for (i in seq_len(nrow(slots))) {
      if (slots$start[i] > ce) {
        cs <- cs + 1L
        ce <- slots$end[i]
      }
      else
        ce <- max(ce, slots$end[i])
      slots$slot[i] <- cs
    }
    if (strand == "-") {
      max_s <- max(slots$slot)
      slots$slot <- max_s - slots$slot + 1L
    }
    g_exons <- g_exons %>% left_join(slots, by = c("start", "end"))
    iso_counts <- g_exons %>%
      group_by(transcript_id) %>%
      summarise(n_ex = n(), .groups = "drop") %>%
      arrange(desc(n_ex))
    if (HAS_REFSEQ && !is.null(refseq_map_df)) {
      iso_counts$nm_id <- refseq_map_df$refseq[match(iso_counts$transcript_id,
                                                     refseq_map_df$ensembl_transcript_id)]
    } else {
      iso_counts$nm_id <- NA_character_
    }
    iso_counts$id_str <- ifelse(
      is.na(iso_counts$nm_id) | iso_counts$nm_id == "",
      iso_counts$transcript_id,
      paste0(iso_counts$transcript_id, " / ", iso_counts$nm_id)
    )
    if (!is.null(tx_annot_df)) {
      iso_counts$tx_length <- tx_annot_df$tx_length[match(iso_counts$transcript_id, tx_annot_df$transcript_id)]
      iso_counts$iso_label <- ifelse(
        is.na(iso_counts$tx_length),
        paste0(iso_counts$id_str, " (", iso_counts$n_ex, " ex)"),
        paste0(
          iso_counts$id_str,
          " (",
          iso_counts$n_ex,
          " ex, ",
          iso_counts$tx_length,
          " nt)"
        )
      )
    } else {
      iso_counts$iso_label <- paste0(iso_counts$id_str, " (", iso_counts$n_ex, " ex)")
    }
    g_exons <- g_exons %>%
      left_join(iso_counts %>% select(transcript_id, iso_label), by = "transcript_id") %>%
      mutate(iso_label = factor(iso_label, levels = rev(iso_counts$iso_label)))
    list(
      exons = g_exons,
      n_slots = max(slots$slot),
      n_iso = nrow(iso_counts),
      gene = gene_name,
      highlight_tx = highlight_tx
    )
  })
  
  output$specific_exon_plot <- renderPlotly({
    res <- specific_exon_data()
    validate(
      need(
        !is.null(res),
        "Click a row in the table above to see this isoform's exon structure in context."
      )
    )
    df <- res$exons
    df$tx_clean <- identity(df$transcript_id)
    df$exon_len <- df$end - df$start + 1
    df$fill_group <- ifelse(df$tx_clean == res$highlight_tx, "highlight", "normal")
    df$fill_group <- factor(df$fill_group, levels = c("normal", "highlight"))
    title_txt <- paste0(
      "Exon Structure: ",
      res$gene,
      " (",
      res$n_iso,
      " isoforms, ",
      res$n_slots,
      " unique exon slots)  |  highlighted: ",
      res$highlight_tx
    )
    step <- max(1, floor(res$n_slots / 20))
    brks <- seq(1, res$n_slots, by = step)
    p <- ggplot(df,
                aes(
                  x = slot,
                  y = iso_label,
                  fill = fill_group,
                  text = paste0(
                    "Isoform: ",
                    transcript_id,
                    "<br>Exon ",
                    exon_number,
                    " of this isoform",
                    "<br>Gene slot: ",
                    slot,
                    " of ",
                    res$n_slots,
                    "<br>Length: ",
                    exon_len,
                    " bp",
                    "<br>",
                    chr,
                    ":",
                    start,
                    "-",
                    end
                  )
                )) +
      geom_tile(
        width = 0.9,
        height = 0.8,
        color = "white",
        linewidth = 0.2
      ) +
      geom_text(
        aes(label = exon_number),
        color = "white",
        size = 2.6,
        fontface = "bold"
      ) +
      scale_fill_manual(
        values = c(normal = "#225ea8", highlight = "#c51b8a"),
        drop = FALSE,
        guide = "none"
      ) +
      scale_x_continuous(breaks = brks, expand = c(0.01, 0.01)) +
      labs(title = title_txt, x = "Exon slot (5' to 3')", y = NULL) +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 9),
        plot.title = element_text(size = 12, face = "bold")
      )
    ggplotly(p, source = "specific_exons", tooltip = "text") %>% event_register("plotly_click") %>% layout(showlegend = FALSE)
  })
  
  output$dl_specific <- downloadHandler(
    filename = function()
      paste0(
        "tissue_specific_",
        input$spec_tissue,
        "_thresh",
        input$spec_threshold,
        "_minTPM",
        input$spec_min_tpm,
        if (HAS_TX_ANNOT &&
            !is.null(input$spec_biotype_filter) &&
            input$spec_biotype_filter != "All")
          paste0("_", input$spec_biotype_filter)
        else
          "",
        ".csv"
      ),
    content = function(file) {
      df <- specific_results()
      if (is.null(df))
        return()
      df <- df %>%
        select(gene, transcript_id, target_tpm, total_tpm, specificity) %>%
        join_annotations()
      if (HAS_TX_ANNOT &&
          !is.null(input$spec_biotype_filter) &&
          input$spec_biotype_filter != "All")
        df <- df %>% filter(Biotype == input$spec_biotype_filter)
      write.csv(df, file, row.names = FALSE)
    }
  )
  
  output$dl_specific_plot <- downloadHandler(
    filename = function()
      paste0("tissue_specific_", input$spec_tissue, ".pdf"),
    content = function(file) {
      df <- specific_results()
      sel <- input$specific_table_rows_selected
      if (is.null(df) || is.null(sel) || length(sel) == 0) {
        showNotification("Select a row in the table first.", type = "warning")
        return()
      }
      tx_id <- df$transcript_id[sel]
      tx_label <- df$display_label[sel]
      tissue_cols <- intersect(names(active_avg()), active_tissues())
      row_data <- active_avg() %>% filter(transcript_id == tx_id)
      if (nrow(row_data) == 0)
        return()
      tpm_vals <- as.numeric(unlist(row_data[1, tissue_cols]))
      plot_df <- tibble(Tissue = tissue_cols, TPM = tpm_vals)
      p <- ggplot(plot_df, aes(
        x = reorder(Tissue, TPM),
        y = TPM,
        fill = Tissue
      )) +
        geom_col() + scale_fill_manual(values = TISSUE_COLORS, labels = legend_label) + coord_flip() + scale_x_discrete(labels = legend_label) +
        labs(
          title = paste("Expression:", tx_label),
          x = NULL,
          y = "Avg TPM"
        ) +
        theme_minimal(base_size = 13) + theme(legend.position = "none")
      ggsave(file, p, width = 16, height = 7)
    }
  )
  
  
  # -------- Tissue-Specific zoom --------
  observe({
    res <- specific_exon_data()
    if (is.null(res))
      return()
    tx_choices <- sort(unique(as.character(res$exons$transcript_id)))
    updateSelectizeInput(
      session,
      "spec_zoom_tx",
      choices = tx_choices,
      selected = if (length(tx_choices) > 0)
        tx_choices[1]
      else
        "",
      server = TRUE
    )
  })
  
  selected_spec_zoom_tx <- reactive({
    click <- event_data("plotly_click", source = "specific_exons")
    res <- specific_exon_data()
    if (!is.null(click) && !is.null(click$y) && !is.null(res)) {
      lvls <- levels(res$exons$iso_label)
      y_val <- click$y
      tx_label <- NULL
      if (is.numeric(y_val)) {
        idx <- as.integer(round(y_val))
        if (idx >= 1 && idx <= length(lvls))
          tx_label <- lvls[idx]
      } else {
        s_y <- as.character(y_val)
        if (s_y %in% lvls)
          tx_label <- s_y
      }
      if (!is.null(tx_label)) {
        tx <- sub(" .*$", "", tx_label)
        return(trimws(tx))
      }
    }
    if (!is.null(input$spec_zoom_tx) &&
        nchar(input$spec_zoom_tx) > 0) {
      return(input$spec_zoom_tx)
    }
    NULL
  })
  
  render_tx_detail <- function(tx) {
    if (is.null(tx) || !HAS_TX_EXONS || is.null(tx_exons_df)) {
      return(
        plotly_empty(type = "scatter", mode = "markers") %>%
          layout(title = "No exon data available")
      )
    }
    ex <- tx_exons_df %>% filter(transcript_id == tx) %>% arrange(start)
    if (nrow(ex) == 0) {
      return(plotly_empty(type = "scatter", mode = "markers") %>%
               layout(title = paste0("No exon data for ", tx)))
    }
    strand <- ex$strand[1]
    chr <- ex$chr[1]
    gene_min <- min(ex$start)
    gene_max <- max(ex$end)
    span_bp <- gene_max - gene_min
    ex <- ex %>% mutate(
      exon_len = end - start + 1,
      next_start = lead(start),
      intron_to_next = next_start - end - 1,
      tip = paste0(
        "Exon ",
        exon_number,
        " of ",
        nrow(ex),
        "<br>Length: ",
        format(exon_len, big.mark = ","),
        " bp",
        "<br>",
        chr,
        ":",
        format(start, big.mark = ","),
        "-",
        format(end, big.mark = ","),
        ifelse(
          is.na(intron_to_next),
          "",
          paste0(
            "<br>Intron to next: ",
            format(intron_to_next, big.mark = ","),
            " bp"
          )
        )
      )
    )
    title_txt <- paste0(
      "Single Transcript: ",
      tx,
      "  (",
      nrow(ex),
      " exons, ",
      format(span_bp, big.mark = ","),
      " bp span on ",
      chr,
      ", strand ",
      strand,
      ")"
    )
    p <- ggplot(ex) +
      geom_segment(
        aes(
          x = gene_min,
          xend = gene_max,
          y = 0,
          yend = 0
        ),
        color = "#888888",
        linewidth = 0.6
      ) +
      geom_rect(
        aes(
          xmin = start,
          xmax = end,
          ymin = -0.45,
          ymax = 0.45,
          text = tip
        ),
        fill = "#225ea8",
        color = "white",
        linewidth = 0.3
      ) +
      geom_text(
        aes(
          x = (start + end) / 2,
          y = 0,
          label = exon_number
        ),
        color = "white",
        size = 3,
        fontface = "bold"
      ) +
      scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
      scale_x_continuous(
        labels = function(x)
          format(x, big.mark = ",", scientific = FALSE),
        expand = expansion(mult = 0.01)
      ) +
      labs(
        title = title_txt,
        x = paste0(chr, " genomic position"),
        y = NULL
      ) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, face = "bold")
      )
    ggplotly(p, tooltip = "text") %>% layout(showlegend = FALSE)
  }
  
  render_tx_detail_ggplot <- function(tx) {
    if (is.null(tx) ||
        !HAS_TX_EXONS || is.null(tx_exons_df))
      return(NULL)
    ex <- tx_exons_df %>% filter(transcript_id == tx) %>% arrange(start)
    if (nrow(ex) == 0)
      return(NULL)
    strand <- ex$strand[1]
    chr <- ex$chr[1]
    gene_min <- min(ex$start)
    gene_max <- max(ex$end)
    span_bp <- gene_max - gene_min
    ggplot(ex) +
      geom_segment(
        aes(
          x = gene_min,
          xend = gene_max,
          y = 0,
          yend = 0
        ),
        color = "#888888",
        linewidth = 0.6
      ) +
      geom_rect(
        aes(
          xmin = start,
          xmax = end,
          ymin = -0.45,
          ymax = 0.45
        ),
        fill = "#225ea8",
        color = "white",
        linewidth = 0.3
      ) +
      geom_text(
        aes(
          x = (start + end) / 2,
          y = 0,
          label = exon_number
        ),
        color = "white",
        size = 3,
        fontface = "bold"
      ) +
      scale_y_continuous(limits = c(-1.1, 1.1), expand = c(0, 0)) +
      scale_x_continuous(
        labels = function(x)
          format(x, big.mark = ",", scientific = FALSE),
        expand = expansion(mult = 0.01)
      ) +
      labs(
        title = paste0(
          "Single Transcript: ",
          tx,
          "  (",
          nrow(ex),
          " exons, ",
          format(span_bp, big.mark = ","),
          " bp span on ",
          chr,
          ", strand ",
          strand,
          ")"
        ),
        x = paste0(chr, " genomic position"),
        y = NULL
      ) +
      theme_minimal(base_size = 11) +
      theme(
        axis.text.y = element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(size = 12, face = "bold")
      )
  }
  
  output$spec_zoom_plot <- renderPlotly({
    tx <- selected_spec_zoom_tx()
    validate(need(
      !is.null(tx),
      "Click an isoform in the grid above or pick one from the dropdown."
    ))
    render_tx_detail(tx)
  })
  
  output$dl_specific_exon_zoom <- downloadHandler(
    filename = function() {
      tx <- selected_spec_zoom_tx()
      paste0("exon_detail_", if (is.null(tx))
        "transcript"
        else
          tx, ".pdf")
    },
    content = function(file) {
      p <- render_tx_detail_ggplot(selected_spec_zoom_tx())
      if (!is.null(p))
        ggsave(file, p, width = 14, height = 3.5)
    }
  )
  
  # -------- Switching zoom --------
  observe({
    res <- switch_exon_data()
    if (is.null(res))
      return()
    tx_choices <- sort(unique(as.character(res$exons$transcript_id)))
    updateSelectizeInput(
      session,
      "switch_zoom_tx",
      choices = tx_choices,
      selected = if (length(tx_choices) > 0)
        tx_choices[1]
      else
        "",
      server = TRUE
    )
  })
  
  selected_switch_zoom_tx <- reactive({
    click <- event_data("plotly_click", source = "switch_exons")
    res <- switch_exon_data()
    if (!is.null(click) && !is.null(click$y) && !is.null(res)) {
      lvls <- levels(res$exons$iso_label)
      y_val <- click$y
      tx_label <- NULL
      if (is.numeric(y_val)) {
        idx <- as.integer(round(y_val))
        if (idx >= 1 && idx <= length(lvls))
          tx_label <- lvls[idx]
      } else {
        s_y <- as.character(y_val)
        if (s_y %in% lvls)
          tx_label <- s_y
      }
      if (!is.null(tx_label)) {
        tx <- sub(" .*$", "", tx_label)
        return(trimws(tx))
      }
    }
    if (!is.null(input$switch_zoom_tx) &&
        nchar(input$switch_zoom_tx) > 0) {
      return(input$switch_zoom_tx)
    }
    NULL
  })
  
  output$switch_zoom_plot <- renderPlotly({
    tx <- selected_switch_zoom_tx()
    validate(need(
      !is.null(tx),
      "Pick a transcript from the dropdown or click a row above."
    ))
    render_tx_detail(tx)
  })
  
  output$dl_switch_exon_zoom <- downloadHandler(
    filename = function() {
      tx <- selected_switch_zoom_tx()
      paste0("exon_detail_", if (is.null(tx))
        "transcript"
        else
          tx, ".pdf")
    },
    content = function(file) {
      p <- render_tx_detail_ggplot(selected_switch_zoom_tx())
      if (!is.null(p))
        ggsave(file, p, width = 14, height = 3.5)
    }
  )
  
  # === Clustering ===
  
  cluster_data <- eventReactive(input$run_cluster, {
    ensure_db()
    done <- log_timing("Clustering")
    on.exit(done())
    use_gene_level <- HAS_GENE_TPM &&
      !is.null(input$clust_level) && input$clust_level == "gene"
    use_3d <- input$clust_3d
    n_comps <- if (use_3d)
      3
    else
      2
    gbm_tmat <- NULL
    gbm_info <- NULL
    if (use_gene_level) {
      gmat <- gene_tpm_df
      sample_data_cols <- intersect(sample_cols, names(gmat))
      if (length(sample_data_cols) == 0) {
        showNotification("No matching sample columns.", type = "error")
        return(NULL)
      }
      num_mat <- as.matrix(gmat[, sample_data_cols, drop = FALSE])
      num_mat[is.na(num_mat)] <- 0
      gene_vars <- apply(num_mat, 1, var, na.rm = TRUE)
      top_idx <- order(gene_vars, decreasing = TRUE)[1:min(input$clust_top_var, nrow(num_mat))]
      num_mat <- num_mat[top_idx, , drop = FALSE]
      gene_ids <- if ("gene_id" %in% names(gmat))
        gmat$gene_id[top_idx]
      else if ("gene" %in% names(gmat))
        gmat$gene[top_idx]
      else
        gmat[[1]][top_idx]
      id_clean <- sub("\\..*", "", gene_ids)
      syms <- tryCatch(
        AnnotationDbi::mapIds(
          org.Hs.eg.db,
          keys = id_clean,
          column = "SYMBOL",
          keytype = "ENSEMBL",
          multiVals = "first"
        ),
        error = function(e)
          setNames(rep(NA, length(id_clean)), id_clean)
      )
      rownames(num_mat) <- ifelse(is.na(syms), gene_ids, syms)
      if (input$clust_tissue_avg) {
        tissue_mat <- sapply(tissue_list, function(t) {
          t_cols <- intersect(
            active_meta() %>% filter(Tissue == t) %>% pull(sample_label),
            colnames(num_mat)
          )
          if (length(t_cols) == 0)
            rep(NA, nrow(num_mat))
          else
            rowMeans(num_mat[, t_cols, drop = FALSE], na.rm = TRUE)
        })
        tmat <- t(tissue_mat)
        tmat[is.na(tmat)] <- 0
        tmat <- log2(tmat + 1)
        sample_info <- tibble(
          sample = rownames(tmat),
          CellLine = rownames(tmat),
          Tissue = rownames(tmat)
        )
      } else {
        active_slabels <- active_meta()$sample_label[active_meta()$Tissue %in% active_tissues()]
        keep_cols <- colnames(num_mat) %in% active_slabels
        tmat <- t(num_mat[, keep_cols, drop = FALSE])
        tmat <- log2(tmat + 1)
        sample_info <- active_meta() %>%
          filter(Tissue %in% active_tissues()) %>%
          select(sample_label, CellLine, Tissue) %>%
          rename(sample = sample_label)
      }
      level_label <- "Gene-Level"
      
      if (!is.null(input$dataset_mode) &&
          active_id() != REF_ID) {
        showNotification(
          "Gene-level clustering shows reference cell lines only. Use Isoform level to project query samples.",
          type = "warning",
          duration = 6
        )
      }
    } else {
      top_tx <- if (input$clust_gene == "__ALL__")
        dbGetQuery(
          db,
          "SELECT transcript_id FROM isoform_variance ORDER BY variance DESC LIMIT ?",
          params = list(input$clust_top_var)
        )$transcript_id
      else
        dbGetQuery(
          db,
          "SELECT transcript_id FROM isoform_variance WHERE gene = ? ORDER BY variance DESC LIMIT ?",
          params = list(input$clust_gene, input$clust_top_var)
        )$transcript_id
      if (length(top_tx) == 0)
        return(NULL)
      if (input$clust_tissue_avg) {
        sub <- active_avg() %>% filter(transcript_id %in% top_tx)
        tissue_cols <- intersect(names(sub), active_tissues())
        tmat <- t(as.matrix(sub[, tissue_cols, drop = FALSE]))
        colnames(tmat) <- paste0(sub$gene, " | ", sub$transcript_id)
        tmat[is.na(tmat)] <- 0
        tmat <- log2(tmat + 1)
        sample_info <- tibble(
          sample = rownames(tmat),
          CellLine = rownames(tmat),
          Tissue = rownames(tmat)
        )
      } else {
        placeholders <- paste(rep("?", length(top_tx)), collapse = ",")
        mat_df <- fetch_quants(paste0("SELECT * FROM quants WHERE transcript_id IN (",
            placeholders,
            ")"
          ),
          params = as.list(top_tx)
        )
        ccle_slabels <- active_meta()$sample_label[active_meta()$Tissue %in% CCLE_TISSUES]
        ccle_cols <- intersect(ccle_slabels, names(mat_df))
        tmat <- t(as.matrix(mat_df[, ccle_cols, drop = FALSE]))
        colnames(tmat) <- paste0(mat_df$gene, " | ", mat_df$transcript_id)
        tmat[is.na(tmat)] <- 0
        tmat <- log2(tmat + 1)
        sample_info <- active_meta() %>%
          filter(Tissue %in% CCLE_TISSUES) %>%
          select(sample_label, CellLine, Tissue) %>%
          rename(sample = sample_label)
        gbm_tmat <- NULL
        gbm_info <- NULL
        if (!is.null(input$dataset_mode) &&
            active_id() != REF_ID) {
          gbm_slabels <- active_meta()$sample_label[active_meta()$Tissue %in% query_tissues()]
          gbm_cols <- intersect(gbm_slabels, names(mat_df))
          if (length(gbm_cols) > 0) {
            gbm_tmat <- t(as.matrix(mat_df[, gbm_cols, drop = FALSE]))
            colnames(gbm_tmat) <- paste0(mat_df$gene, " | ", mat_df$transcript_id)
            gbm_tmat[is.na(gbm_tmat)] <- 0
            gbm_tmat <- log2(gbm_tmat + 1)
            gbm_info <- active_meta() %>%
              filter(Tissue %in% query_tissues()) %>%
              select(sample_label, CellLine, Tissue) %>%
              rename(sample = sample_label)
          }
        }
      }
      level_label <- "Isoform-Level"
    }
    if (length(input$clust_tissue_filter) > 0 &&
        !input$clust_tissue_avg) {
      keep_samples <- sample_info %>%
        filter(Tissue %in% input$clust_tissue_filter) %>% pull(sample)
      keep_idx <- which(rownames(tmat) %in% keep_samples)
      need_n <- if (input$clust_reduction == "pca") 3 else 5
      if (length(keep_idx) < need_n) {
        showNotification(paste0("Clustering needs at least ", need_n, " samples for ",
                                toupper(input$clust_reduction), "; the selected tissues have ",
                                length(keep_idx), ". Add more tissues or groups."),
                         type = "warning", duration = 8)
        return(NULL)
      }
      tmat <- tmat[keep_idx, , drop = FALSE]
      sample_info <- sample_info %>% filter(sample %in% rownames(tmat))
    }
    if (length(input$clust_tissue_filter) > 0 &&
        input$clust_tissue_avg) {
      keep_idx <- which(rownames(tmat) %in% input$clust_tissue_filter)
      if (length(keep_idx) < 3) {
        showNotification("Group-average clustering needs at least 3 tissues or groups selected.", type = "warning", duration = 8)
        return(NULL)
      }
      tmat <- tmat[keep_idx, , drop = FALSE]
      sample_info <- sample_info %>% filter(sample %in% rownames(tmat))
    }
    min_n <- if (input$clust_reduction == "pca") 3 else 5
    if (nrow(tmat) < min_n) {
      showNotification(paste0("Clustering needs at least ", min_n, " samples for ",
                              toupper(input$clust_reduction), "; the current selection has ",
                              nrow(tmat), ". Add more tissues or groups."),
                       type = "warning", duration = 8)
      return(NULL)
    }
    keep_feat <- apply(tmat, 2, function(v) isTRUE(var(v) > 0))
    if (sum(keep_feat) < 2) {
      showNotification("Too few variable features in this selection to cluster.", type = "warning")
      return(NULL)
    }
    tmat <- tmat[, keep_feat, drop = FALSE]
    if (!is.null(gbm_tmat) && ncol(gbm_tmat) == length(keep_feat)) gbm_tmat <- gbm_tmat[, keep_feat, drop = FALSE]
    if (input$clust_reduction == "pca") {
      pca <- prcomp(tmat,
                    scale. = (input$clust_method == "cor"),
                    center = TRUE)
      n_keep <- min(max(n_comps, 5), ncol(pca$x))
      embed_df <- as.data.frame(pca$x[, 1:n_keep])
      names(embed_df)[1:2] <- c("Dim1", "Dim2")
      if (use_3d && ncol(pca$x) >= 3)
        names(embed_df)[3] <- "Dim3"
      embed_df$sample <- rownames(embed_df)
      if (!is.null(gbm_tmat) && !is.null(gbm_info)) {
        gbm_scaled <- scale(gbm_tmat,
                            center = pca$center,
                            scale = if (input$clust_method == "cor")
                              pca$scale
                            else
                              FALSE)
        gbm_proj <- as.data.frame(gbm_scaled %*% pca$rotation[, 1:n_keep])
        names(gbm_proj)[1:2] <- c("Dim1", "Dim2")
        if (use_3d &&
            ncol(gbm_proj) >= 3)
          names(gbm_proj)[3] <- "Dim3"
        gbm_proj$sample <- rownames(gbm_tmat)
        embed_df <- rbind(embed_df, gbm_proj)
        sample_info <- rbind(sample_info, gbm_info)
      }
      ve <- summary(pca)$importance[2, 1:min(10, ncol(pca$x))] * 100
      embed_df <- embed_df %>% left_join(sample_info, by = "sample")
      embed_df$is_gbm <- embed_df$Tissue %in% query_tissues()
      n_pcs_for_mag <- if (use_3d)
        3
      else
        2
      n_pcs_avail <- min(n_pcs_for_mag, ncol(pca$rotation))
      loadings_mat <- pca$rotation[, 1:n_pcs_avail, drop = FALSE]
      feat_mag <- sqrt(rowSums(loadings_mat^2))
      top_n_feat <- min(100, length(feat_mag))
      top_idx_feat <- order(feat_mag, decreasing = TRUE)[1:top_n_feat]
      top_features <- data.frame(
        feature = rownames(loadings_mat)[top_idx_feat],
        PC1_loading = round(loadings_mat[top_idx_feat, 1], 4),
        PC2_loading = if (n_pcs_avail >= 2)
          round(loadings_mat[top_idx_feat, 2], 4)
        else
          NA,
        PC3_loading = if (n_pcs_avail >= 3)
          round(loadings_mat[top_idx_feat, 3], 4)
        else
          NA,
        magnitude = round(feat_mag[top_idx_feat], 4),
        stringsAsFactors = FALSE
      )
      if (!use_3d)
        top_features$PC3_loading <- NULL
      list(
        embed_df = embed_df,
        var_explained = ve,
        pca_obj = pca,
        top_features = top_features,
        ax1 = paste0("PC1 (", round(ve[1], 1), "%)"),
        ax2 = paste0("PC2 (", round(ve[2], 1), "%)"),
        ax3 = if (use_3d && length(ve) >= 3)
          paste0("PC3 (", round(ve[3], 1), "%)")
        else
          NULL,
        method = "pca",
        tissue_avg = input$clust_tissue_avg,
        is_3d = use_3d,
        level = level_label
      )
    } else {
      umap_model <- umap(
        tmat,
        n_neighbors = min(input$umap_neighbors, nrow(tmat) - 1),
        min_dist = input$umap_min_dist,
        n_components = n_comps,
        metric = "cosine",
        ret_model = TRUE
      )
      embed_df <- as.data.frame(umap_model$embedding)
      names(embed_df)[1:2] <- c("Dim1", "Dim2")
      if (use_3d && ncol(embed_df) >= 3)
        names(embed_df)[3] <- "Dim3"
      embed_df$sample <- rownames(tmat)
      if (!is.null(gbm_tmat) && !is.null(gbm_info)) {
        gbm_proj <- as.data.frame(umap_transform(gbm_tmat, umap_model))
        names(gbm_proj)[1:2] <- c("Dim1", "Dim2")
        if (use_3d &&
            ncol(gbm_proj) >= 3)
          names(gbm_proj)[3] <- "Dim3"
        gbm_proj$sample <- rownames(gbm_tmat)
        embed_df <- rbind(embed_df, gbm_proj)
        sample_info <- rbind(sample_info, gbm_info)
      }
      embed_df <- embed_df %>% left_join(sample_info, by = "sample")
      embed_df$is_gbm <- embed_df$Tissue %in% query_tissues()
      feat_var_u <- apply(tmat, 2, var, na.rm = TRUE)
      feat_mean_u <- colMeans(tmat, na.rm = TRUE)
      top_n_feat_u <- min(100, length(feat_var_u))
      top_idx_feat_u <- order(feat_var_u, decreasing = TRUE)[1:top_n_feat_u]
      top_features <- data.frame(
        feature = colnames(tmat)[top_idx_feat_u],
        variance = round(feat_var_u[top_idx_feat_u], 4),
        mean_log2tpm = round(feat_mean_u[top_idx_feat_u], 4),
        stringsAsFactors = FALSE
      )
      list(
        embed_df = embed_df,
        var_explained = NULL,
        top_features = top_features,
        ax1 = "UMAP 1",
        ax2 = "UMAP 2",
        ax3 = if (use_3d)
          "UMAP 3"
        else
          NULL,
        method = "umap",
        tissue_avg = input$clust_tissue_avg,
        is_3d = use_3d,
        level = level_label
      )
    }
  })
  
  output$cluster_plot <- renderPlotly({
    res <- cluster_data()
    if (is.null(res))
      return(NULL)
    df <- res$embed_df
    pt_size <- if (res$tissue_avg)
      input$clust_dot_size * 2
    else
      input$clust_dot_size
    title <- paste0(
      ifelse(res$method == "pca", "PCA", "UMAP"),
      ifelse(res$is_3d, " 3D", ""),
      " \u2014 ",
      ifelse(res$tissue_avg, "Tissue Averages", "Samples"),
      " [",
      res$level,
      "]"
    )
    if (res$is_3d && "Dim3" %in% names(df)) {
      plot_ly(
        df,
        x = ~ Dim1,
        y = ~ Dim2,
        z = ~ Dim3,
        color = ~ legend_label(Tissue),
        colors = setNames(TISSUE_COLORS[sort(unique(df$Tissue))], legend_label(sort(unique(df$Tissue)))),
        text = ~ paste0(
          ifelse(res$tissue_avg, "Tissue: ", "Sample: "),
          CellLine,
          "\nTissue: ",
          Tissue
        ),
        hoverinfo = "text",
        type = "scatter3d",
        mode = "markers",
        marker = list(
          size = if (res$tissue_avg)
            input$clust_dot_size * 2
          else
            input$clust_dot_size,
          opacity = 0.8
        )
      ) %>%
        layout(
          legend = list(font = list(size = 13), itemsizing = "constant"),
          title = title,
          scene = list(
            xaxis = list(
              title = res$ax1,
              zeroline = FALSE,
              tickfont = list(size = 12)
            ),
            yaxis = list(
              title = res$ax2,
              zeroline = FALSE,
              tickfont = list(size = 12)
            ),
            zaxis = list(
              title = res$ax3,
              zeroline = FALSE,
              tickfont = list(size = 12)
            )
          )
        )
    } else {
      ccle_df <- df[!df$is_gbm, ]
      gbm_df <- df[df$is_gbm, ]
      p <- ggplot(ccle_df, aes(
        x = Dim1,
        y = Dim2,
        color = Tissue,
        text = paste0(
          ifelse(res$tissue_avg, "Tissue: ", "Sample: "),
          CellLine,
          "\nTissue: ",
          Tissue,
          "\n",
          res$ax1,
          ": ",
          round(Dim1, 2),
          "\n",
          res$ax2,
          ": ",
          round(Dim2, 2)
        )
      )) +
        geom_point(size = pt_size, alpha = 0.8) +
        scale_color_manual(values = TISSUE_COLORS, labels = legend_label)
      if (nrow(gbm_df) > 0) {
        p <- p + geom_point(
          data = gbm_df,
          aes(
            x = Dim1,
            y = Dim2,
            text = paste0(
              "Patient: ",
              CellLine,
              "\nCondition: ",
              Tissue,
              "\n",
              res$ax1,
              ": ",
              round(Dim1, 2),
              "\n",
              res$ax2,
              ": ",
              round(Dim2, 2)
            )
          ),
          color = "#c51b8a",
          fill = "#c51b8a",
          shape = 18,
          size = pt_size * 1.3,
          alpha = 0.95,
          inherit.aes = FALSE
        )
      }
      p <- p +
        labs(
          title = title,
          x = res$ax1,
          y = res$ax2,
          color = "Tissue"
        ) +
        theme_minimal() +
        guides(color = guide_legend(override.aes = list(size = 5))) +
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 13))
      if (res$tissue_avg)
        p <- p + geom_text(
          aes(label = Tissue),
          size = 2.5,
          vjust = -1.2,
          show.legend = FALSE
        )
      ggplotly(p, tooltip = "text") %>% toWebGL() %>%
        layout(
          legend = list(font = list(size = 13), itemsizing = "constant"),
          xaxis = list(zeroline = FALSE),
          yaxis = list(zeroline = FALSE)
        )
    }
  })
  
  output$scree_plot <- renderPlotly({
    res <- cluster_data()
    if (is.null(res) || is.null(res$var_explained))
      return(NULL)
    ve <- res$var_explained
    df <- tibble(PC = factor(paste0("PC", seq_along(ve)), levels = paste0("PC", seq_along(ve))),
                 Variance = ve)
    p <- ggplot(df, aes(
      x = PC,
      y = Variance,
      text = paste0(PC, ": ", round(Variance, 1), "%")
    )) +
      geom_col(fill = "#225ea8") +
      labs(title = "Variance Explained", x = NULL, y = "% Variance") +
      theme_minimal()
    ggplotly(p, tooltip = "text")
  })
  
  # Task 37: Clustering PDF download
  output$cluster_top_features_table <- renderDT({
    res <- cluster_data()
    if (is.null(res) || is.null(res$top_features))
      return(NULL)
    datatable(
      res$top_features,
      options = list(
        pageLength = 15,
        scrollX = TRUE,
        order = list(list(ncol(
          res$top_features
        ) - 1, "desc"))
      ),
      rownames = FALSE,
      caption = paste0(
        "Top ",
        nrow(res$top_features),
        " features driving ",
        toupper(res$method),
        " (",
        res$level,
        ")"
      )
    )
  })
  
  output$dl_cluster_top_features <- downloadHandler(
    filename = function() {
      res <- cluster_data()
      method_lbl <- if (is.null(res))
        "clustering"
      else
        res$method
      paste0("top_features_", method_lbl, ".csv")
    },
    content = function(file) {
      res <- cluster_data()
      if (is.null(res) || is.null(res$top_features)) {
        write.csv(data.frame(message = "No clustering results available"),
                  file,
                  row.names = FALSE)
        return()
      }
      write.csv(res$top_features, file, row.names = FALSE)
    }
  )
  
  output$dl_cluster_plot <- downloadHandler(
    filename = function()
      paste0(
        "clustering_",
        input$clust_reduction,
        if (input$clust_gene != "__ALL__")
          paste0("_", input$clust_gene)
        else
          "_allgenes",
        ".pdf"
      ),
    content = function(file) {
      res <- cluster_data()
      if (is.null(res) || res$is_3d)
        return()
      df <- res$embed_df
      pt_size <- if (res$tissue_avg)
        input$clust_dot_size * 2
      else
        input$clust_dot_size
      p <- ggplot(df, aes(x = Dim1, y = Dim2, color = Tissue)) +
        geom_point(size = pt_size, alpha = 0.8) +
        scale_color_manual(values = TISSUE_COLORS, labels = legend_label) +
        labs(
          title = paste0(
            ifelse(res$method == "pca", "PCA", "UMAP"),
            " -- ",
            ifelse(res$tissue_avg, "Tissue Averages", "Samples"),
            " [",
            res$level,
            "]"
          ),
          x = res$ax1,
          y = res$ax2,
          color = "Tissue"
        ) +
        theme_minimal() +
        guides(color = guide_legend(override.aes = list(size = 5))) +
        theme(legend.text = element_text(size = 12),
              legend.title = element_text(size = 10))
      if (res$tissue_avg)
        p <- p + geom_text(
          aes(label = Tissue),
          size = 2.5,
          vjust = -1.2,
          show.legend = FALSE
        )
      ggsave(file, p, width = 12, height = 8)
    }
  )
  
  # Task 38: Top features driving clustering
  output$cluster_features_table <- renderDT({
    res <- cluster_data()
    if (is.null(res))
      return(NULL)
    if (res$method == "pca" && !is.null(res$pca_obj)) {
      pca_obj <- res$pca_obj
      loads <- pca_obj$rotation[, 1:min(2, ncol(pca_obj$rotation)), drop = FALSE]
      if (nrow(loads) == 0)
        return(NULL)
      feat_names <- rownames(loads)
      if (is.null(feat_names))
        feat_names <- paste0("Feature_", seq_len(nrow(loads)))
      feat_df <- data.frame(
        feature = feat_names,
        PC1_loading = loads[, 1],
        stringsAsFactors = FALSE
      )
      if (ncol(loads) >= 2) {
        feat_df$PC2_loading <- loads[, 2]
        feat_df$importance <- abs(loads[, 1]) + abs(loads[, 2])
      } else {
        feat_df$PC2_loading <- NA_real_
        feat_df$importance <- abs(loads[, 1])
      }
      feat_df <- feat_df[order(-feat_df$importance), ]
      feat_df <- head(feat_df, 50)
      datatable(
        feat_df,
        options = list(pageLength = 20, scrollX = TRUE),
        rownames = FALSE
      ) %>%
        formatRound(
          columns = c("PC1_loading", "PC2_loading", "importance"),
          digits = 4
        )
    } else {
      NULL
    }
  })
  
  # === Tissue Correlation ===
  
  corr_data <- eventReactive(input$run_corr, {
    ensure_db()
    top_tx <- if (input$corr_gene == "__ALL__")
      dbGetQuery(
        db,
        "SELECT transcript_id FROM isoform_variance ORDER BY variance DESC LIMIT ?",
        params = list(input$corr_top_var)
      )$transcript_id
    else
      dbGetQuery(
        db,
        "SELECT transcript_id FROM isoform_variance WHERE gene = ? ORDER BY variance DESC LIMIT ?",
        params = list(input$corr_gene, input$corr_top_var)
      )$transcript_id
    if (length(top_tx) == 0)
      return(NULL)
    if (HAS_TX_ANNOT &&
        !is.null(input$corr_biotype_filter) &&
        input$corr_biotype_filter != "All") {
      bt_tx <- tx_annot_df$transcript_id[tx_annot_df$biotype == input$corr_biotype_filter]
      top_tx_clean <- identity(top_tx)
      keep <- top_tx_clean %in% bt_tx
      top_tx <- top_tx[keep]
      if (length(top_tx) == 0)
        return(NULL)
    }
    if (input$corr_mode == "tissue") {
      sub <- active_avg() %>% filter(transcript_id %in% top_tx)
      tissue_cols <- intersect(names(sub), active_tissues())
      if (length(input$corr_tissue_filter) > 0)
        tissue_cols <- intersect(tissue_cols, input$corr_tissue_filter)
      if (length(tissue_cols) < 2)
        return(NULL)
      mat <- log2(as.matrix(sub[, tissue_cols, drop = FALSE]) + 1)
      list(
        cor_mat = cor(mat, method = input$corr_method),
        mode = "tissue",
        labels = tissue_cols
      )
    } else {
      placeholders <- paste(rep("?", length(top_tx)), collapse = ",")
      mat_df <- fetch_quants(paste0("SELECT * FROM quants WHERE transcript_id IN (",
          placeholders,
          ")"
        ),
        params = as.list(top_tx)
      )
      sample_info <- active_meta() %>% select(sample_label, CellLine, Tissue)
      if (length(input$corr_tissue_filter) > 0)
        sample_info <- sample_info %>% filter(Tissue %in% input$corr_tissue_filter)
      col_order <- intersect(sample_info$sample_label, names(mat_df))
      if (length(col_order) < 2)
        return(NULL)
      mat <- log2(as.matrix(mat_df[, col_order, drop = FALSE]) + 1)
      cor_mat <- cor(mat, method = input$corr_method)
      cl_names <- sample_info$CellLine[match(col_order, sample_info$sample_label)]
      cl_tissues <- sample_info$Tissue[match(col_order, sample_info$sample_label)]
      rownames(cor_mat) <- cl_names
      colnames(cor_mat) <- cl_names
      tissue_order <- order(cl_tissues, cl_names)
      cor_mat <- cor_mat[tissue_order, tissue_order]
      cl_tissues <- cl_tissues[tissue_order]
      list(
        cor_mat = cor_mat,
        mode = "cellline",
        labels = rownames(cor_mat),
        tissues = cl_tissues
      )
    }
  })
  
  output$corr_heatmap <- renderPlotly({
    res <- corr_data()
    if (is.null(res))
      return(NULL)
    cor_mat <- res$cor_mat
    if (res$mode == "tissue") {
      hc <- hclust(as.dist(1 - cor_mat))
      ord <- hc$labels[hc$order]
      cor_vals <- cor_mat[upper.tri(cor_mat)]
      cor_low  <- floor(min(cor_vals, na.rm = TRUE) * 10) / 10
      cor_mid  <- median(cor_vals, na.rm = TRUE)
      cor_df <- as.data.frame(as.table(cor_mat))
      names(cor_df) <- c("Tissue1", "Tissue2", "Correlation")
      cor_df$Tissue1 <- factor(cor_df$Tissue1, levels = ord)
      cor_df$Tissue2 <- factor(cor_df$Tissue2, levels = ord)
      p <- ggplot(cor_df,
                  aes(
                    x = Tissue1,
                    y = Tissue2,
                    fill = Correlation,
                    text = paste0(
                      Tissue1,
                      " vs ",
                      Tissue2,
                      "\n",
                      input$corr_method,
                      " r = ",
                      round(Correlation, 3)
                    )
                  )) +
        geom_tile(color = "white", linewidth = 0.3) +
        scale_fill_gradient2(
          low = "#225ea8",
          mid = "white",
          high = "#b2182b",
          midpoint = cor_mid,
          limits = c(cor_low, 1)
        ) +
        labs(
          title = paste0("Tissue Correlation (", input$corr_method, ")"),
          x = NULL,
          y = NULL
        ) +
        theme_minimal(base_size = 13) +
        theme(axis.text.x = element_text(
          angle = 45,
          hjust = 1,
          size = 13
        ),
        panel.grid = element_blank())
      ggplotly(p, tooltip = "text")
    } else {
      n <- nrow(cor_mat)
      cor_vals_cl <- cor_mat[upper.tri(cor_mat)]
      cor_low_cl  <- floor(min(cor_vals_cl, na.rm = TRUE) * 10) / 10
      plot_ly(
        x = rownames(cor_mat),
        y = rownames(cor_mat),
        z = cor_mat,
        type = "heatmap",
        colorscale = list(c(0, "#225ea8"), c(0.5, "#ffffff"), c(1, "#b2182b")),
        zmin = cor_low_cl,
        zmax = 1,
        hovertemplate = "%{x} vs %{y}<br>r = %{z:.3f}<extra></extra>"
      ) %>%
        layout(
          title = paste0("Sample Correlation (", input$corr_method, ", n=", n, ")"),
          xaxis = list(
            showticklabels = n <= 50,
            tickangle = 45,
            tickfont = list(size = 9)
          ),
          yaxis = list(
            showticklabels = n <= 50,
            tickfont = list(size = 9)
          )
        )
    }
  })
  
  output$corr_table <- renderDT({
    res <- corr_data()
    if (is.null(res))
      return(NULL)
    cor_df <- as.data.frame(as.table(res$cor_mat))
    names(cor_df) <- c("Sample1", "Sample2", "Correlation")
    cor_df <- cor_df %>%
      filter(as.character(Sample1) < as.character(Sample2)) %>%
      arrange(Correlation)
    if (res$mode == "cellline" && !is.null(res$tissues)) {
      tissue_lookup <- setNames(res$tissues, res$labels)
      cor_df <- cor_df %>%
        mutate(Tissue1 = tissue_lookup[as.character(Sample1)], Tissue2 = tissue_lookup[as.character(Sample2)]) %>%
        select(Sample1, Tissue1, Sample2, Tissue2, Correlation)
    }
    datatable(cor_df,
              options = list(pageLength = 20, scrollX = TRUE),
              rownames = FALSE) %>%
      formatRound(columns = "Correlation", digits = 3)
  })
  
  output$dl_corr_plot <- downloadHandler(
    filename = function()
      paste0(
        "correlation_",
        input$corr_mode,
        "_",
        input$corr_method,
        ".pdf"
      ),
    content = function(file) {
      res <- corr_data()
      if (is.null(res))
        return()
      if (res$mode == "tissue") {
        cor_mat <- res$cor_mat
        hc <- hclust(as.dist(1 - cor_mat))
        ord <- hc$labels[hc$order]
        cor_vals <- cor_mat[upper.tri(cor_mat)]
        cor_low <- floor(min(cor_vals, na.rm = TRUE) * 10) / 10
        cor_mid <- median(cor_vals, na.rm = TRUE)
        cor_df <- as.data.frame(as.table(cor_mat))
        names(cor_df) <- c("Tissue1", "Tissue2", "Correlation")
        cor_df$Tissue1 <- factor(cor_df$Tissue1, levels = ord)
        cor_df$Tissue2 <- factor(cor_df$Tissue2, levels = ord)
        p <- ggplot(cor_df, aes(x = Tissue1, y = Tissue2, fill = Correlation)) +
          geom_tile(color = "white", linewidth = 0.3) +
          scale_fill_gradient2(
            low = "#225ea8",
            mid = "white",
            high = "#b2182b",
            midpoint = cor_mid,
            limits = c(cor_low, 1)
          ) +
          labs(
            title = paste0("Tissue Correlation (", input$corr_method, ")"),
            x = NULL,
            y = NULL
          ) +
          theme_minimal(base_size = 13) +
          theme(axis.text.x = element_text(
            angle = 45,
            hjust = 1,
            size = 10
          ),
          panel.grid = element_blank())
        ggsave(file, p, width = 10, height = 9)
      }
    }
  )
  
  output$dl_corr_table <- downloadHandler(
    filename = function()
      paste0(
        "correlation_",
        input$corr_mode,
        "_",
        input$corr_method,
        "_top",
        input$corr_top_var,
        ".csv"
      ),
    content = function(file) {
      res <- corr_data()
      if (is.null(res))
        return()
      cor_df <- as.data.frame(as.table(res$cor_mat))
      names(cor_df) <- c("Sample1", "Sample2", "Correlation")
      cor_df <- cor_df %>%
        filter(as.character(Sample1) < as.character(Sample2)) %>%
        arrange(Correlation)
      if (res$mode == "cellline" && !is.null(res$tissues)) {
        tissue_lookup <- setNames(res$tissues, res$labels)
        cor_df <- cor_df %>%
          mutate(Tissue1 = tissue_lookup[as.character(Sample1)], Tissue2 = tissue_lookup[as.character(Sample2)]) %>%
          select(Sample1, Tissue1, Sample2, Tissue2, Correlation)
      }
      write.csv(cor_df, file, row.names = FALSE)
    }
  )
  
  # === Differential Expression ===
  
  de_results <- reactiveVal(NULL)
  de_active_method <- reactiveVal("wilcoxon")
  
  output$de_method_info <- renderUI({
    method <- de_active_method()
    if (method == "wilcoxon")
      return(NULL)
    label <- if (method == "deseq2_gene")
      "DESeq2 Gene-Level (precomputed from tximport gene counts)"
    else
      "DESeq2 Transcript-Level, isoforms (precomputed from Salmon counts)"
    tags$div(style = "background: #eef4ff; padding: 8px 12px; border-radius: 6px; margin-bottom: 10px; font-size: 13px;", tags$strong("Method: "), label)
  })
  
  observeEvent(input$run_de, {
    ensure_db()
    done <- log_timing(paste("Differential expression", input$de_method))
    on.exit(done())
    if (input$de_mode == "tissue" && identical(input$de_tissue1, input$de_tissue2)) {
      showNotification("Pick two different groups.", type = "warning")
      return()
    }
    method <- input$de_method
    if (method %in% c("deseq2_gene", "deseq2_tx")) {
      if (input$de_mode == "cellline") {
        showNotification("DESeq2 requires replicates. Use Wilcoxon.", type = "warning")
        return()
      }
      withProgress(message = "Querying DESeq2 results...", {
        in_q <- c(input$de_tissue1, input$de_tissue2) %in% query_tissues()
        if (any(in_q) && !all(in_q)) {
          showNotification("DESeq2 results exist only within one dataset. Use Wilcoxon to compare across datasets.",
                           type = "warning", duration = 8)
          de_results(NULL)
          return()
        }
        con <- if (all(in_q)) deseq_con(active_id()) else deseq_db
        if (is.null(con)) {
          showNotification("No DESeq2 database for this dataset.", type = "error")
          de_results(NULL)
          return()
        }
        t1_db <- if (all(in_q)) raw_group(input$de_tissue1) else unrename_tissue(input$de_tissue1)
        t2_db <- if (all(in_q)) raw_group(input$de_tissue2) else unrename_tissue(input$de_tissue2)
        table_name <- if (method == "deseq2_gene")
          "gene_results"
        else
          "transcript_results"
        setProgress(value = 0.3,
                    detail = paste0("Fetching from ", table_name, "..."))
        res <- dbGetQuery(
          con,
          paste0(
            "SELECT * FROM ",
            table_name,
            " WHERE group1 = ? AND group2 = ?"
          ),
          params = list(t1_db, t2_db)
        )
        flipped <- FALSE
        if (nrow(res) == 0) {
          res <- dbGetQuery(
            con,
            paste0(
              "SELECT * FROM ",
              table_name,
              " WHERE group1 = ? AND group2 = ?"
            ),
            params = list(t2_db, t1_db)
          )
          if (nrow(res) > 0)
            flipped <- TRUE
        }
        if (nrow(res) == 0) {
          showNotification("No DESeq2 results for this pair.",
                           type = "warning",
                           duration = 8)
          de_results(NULL)
          return()
        }
        setProgress(value = 0.7,
                    detail = paste0(nrow(res), " results loaded"))
        num_cols <- c("baseMean",
                      "log2FoldChange",
                      "lfcSE",
                      "stat",
                      "pvalue",
                      "padj")
        for (col in intersect(num_cols, names(res)))
          res[[col]] <- as.numeric(res[[col]])
        res <- res %>% filter(!is.na(padj))
        if (flipped && "log2FoldChange" %in% names(res))
          res$log2FoldChange <- -res$log2FoldChange
        if ("log2FoldChange" %in% names(res))
          res <- res %>% dplyr::rename(log2FC = log2FoldChange)
        if ("feature_id" %in% names(res) && method == "deseq2_gene")
          res <- res %>% dplyr::rename(gene_id = feature_id, gene = gene_symbol)
        else if ("feature_id" %in% names(res) &&
                 method == "deseq2_tx")
          res <- res %>% dplyr::rename(transcript_id = feature_id, gene = gene_symbol)
        res <- res %>%
          select(-any_of(c(
            "tissue1", "tissue2", "group1", "group2", "level", "gene_id_clean"
          ))) %>%
          arrange(padj, desc(abs(log2FC)))
        de_results(res)
        de_active_method(method)
        showNotification(paste0("DESeq2: ", nrow(res), " results loaded."), type = "message")
      })
      return()
    }
    de_active_method("wilcoxon")
    withProgress(message = "Running DE analysis...", {
      if (input$de_mode == "tissue") {
        t1_samples <- active_meta() %>%
          filter(Tissue == input$de_tissue1) %>% pull(sample_label)
        t2_samples <- active_meta() %>%
          filter(Tissue == input$de_tissue2) %>% pull(sample_label)
        cols1 <- intersect(t1_samples, quant_fields())
        cols2 <- intersect(t2_samples, quant_fields())
        if (length(cols1) < 2 || length(cols2) < 2) {
          de_results(NULL)
          return()
        }
        setProgress(value = 0.1, detail = "Fetching expression data...")
        full_df <- fetch_quant_cols(c(cols1, cols2))
        mat1 <- as.matrix(full_df[, cols1, drop = FALSE])
        mat2 <- as.matrix(full_df[, cols2, drop = FALSE])
        mean1 <- rowMeans(mat1, na.rm = TRUE)
        mean2 <- rowMeans(mat2, na.rm = TRUE)
        keep <- expr_ok(mean1, mean2, input$de_min_tpm, input$de_expr_rule)
        full_df <- full_df[keep, ]
        mat1 <- mat1[keep, , drop = FALSE]
        mat2 <- mat2[keep, , drop = FALSE]
        mean1 <- mean1[keep]
        mean2 <- mean2[keep]
        setProgress(value = 0.3,
                    detail = paste0("Testing ", nrow(full_df), " transcripts..."))
        log2fc <- log2((mean1 + 1) / (mean2 + 1))
        pvals <- sapply(seq_len(nrow(mat1)), function(i)
          tryCatch(
            wilcox.test(mat1[i, ], mat2[i, ], exact = FALSE)$p.value,
            error = function(e)
              NA_real_
          ))
        setProgress(value = 0.8, detail = "Adjusting p-values...")
        res <- data.frame(
          gene = full_df$gene,
          transcript_id = full_df$transcript_id,
          mean_1 = mean1,
          mean_2 = mean2,
          log2FC = log2fc,
          pvalue = pvals,
          stringsAsFactors = FALSE
        )
        res$padj <- p.adjust(res$pvalue, method = "BH")
        res <- res %>% arrange(padj, desc(abs(log2FC)))
        names(res)[names(res) == "mean_1"] <- paste0("mean_", input$de_tissue1)
        names(res)[names(res) == "mean_2"] <- paste0("mean_", input$de_tissue2)
        de_results(res)
      } else {
        cl1_label <- active_meta() %>%
          filter(CellLine == input$de_cl1) %>% pull(sample_label) %>% head(1)
        cl2_label <- active_meta() %>%
          filter(CellLine == input$de_cl2) %>% pull(sample_label) %>% head(1)
        if (length(cl1_label) == 0 || length(cl2_label) == 0) {
          de_results(NULL)
          return()
        }
        all_df <- fetch_quant_cols(c(cl1_label, cl2_label))
        names(all_df)[3:4] <- c("tpm_1", "tpm_2")
        all_df <- all_df %>%
          filter(expr_ok(tpm_1, tpm_2, input$de_min_tpm, input$de_expr_rule)) %>%
          mutate(log2FC = log2((tpm_1 + 1) / (tpm_2 + 1))) %>%
          arrange(desc(abs(log2FC)))
        names(all_df)[names(all_df) == "tpm_1"] <- paste0("TPM_", input$de_cl1)
        names(all_df)[names(all_df) == "tpm_2"] <- paste0("TPM_", input$de_cl2)
        de_results(all_df)
      }
    })
  })
  
  output$de_plot <- renderPlotly({
    res <- de_results()
    if (is.null(res) || nrow(res) == 0)
      return(NULL)
    method <- de_active_method()
    if ("padj" %in% names(res) &&
        (method != "wilcoxon" || input$de_mode == "tissue")) {
      fc_cut <- input$de_volcano_fc
      pval_cut <- input$de_volcano_pval
      res <- res %>% mutate(
        sig_category = case_when(
          padj < pval_cut & abs(log2FC) >= fc_cut ~ "Significant",
          padj < pval_cut ~ "Sig. p-value only",
          TRUE ~ "Not significant"
        )
      )
      color_map <- c(
        "Significant" = BERRY_SIG,
        "Sig. p-value only" = BERRY_LIGHT,
        "Not significant" = GREY_NONSIG
      )
      id_col <- if ("transcript_id" %in% names(res))
        "transcript_id"
      else if ("gene_id" %in% names(res))
        "gene_id"
      else
        NULL
      gene_col <- if ("gene" %in% names(res))
        "gene"
      else if ("gene_name" %in% names(res))
        "gene_name"
      else
        NULL
      res$hover_text <- paste0(
        if (!is.null(gene_col))
          paste0("Gene: ", res[[gene_col]], "\n")
        else
          "",
        if (!is.null(id_col))
          paste0("ID: ", res[[id_col]], "\n")
        else
          "",
        "log2FC: ",
        round(res$log2FC, 2),
        "\nAdj. p: ",
        signif(res$padj, 3)
      )
      method_label <- switch(method,
                             deseq2_gene = "DESeq2 Gene",
                             deseq2_tx = "DESeq2 Transcript",
                             "Wilcoxon Transcript")
      p <- ggplot(res,
                  aes(
                    x = log2FC,
                    y = -log10(padj),
                    color = sig_category,
                    text = hover_text
                  )) +
        geom_point(size = input$de_legend_size * 0.4, alpha = 0.7) +
        guides(color = guide_legend(override.aes = list(size = input$de_legend_size))) +
        scale_color_manual(values = color_map) +
        geom_hline(
          yintercept = -log10(pval_cut),
          linetype = "dashed",
          color = "grey50",
          linewidth = 0.3
        ) +
        geom_vline(
          xintercept = c(-fc_cut, fc_cut),
          linetype = "dashed",
          color = "grey50",
          linewidth = 0.3
        ) +
        labs(
          title = paste0(
            method_label,
            " DE: ",
            input$de_tissue1,
            " vs ",
            input$de_tissue2
          ),
          x = "log2 Fold Change",
          y = "-log10(Adj. p-value)",
          color = NULL
        ) +
        theme_minimal(base_size = 13) + theme(legend.position = "bottom")
      ggplotly(p, tooltip = "text") %>% toWebGL() %>%
        layout(
          legend = list(
            font = list(size = 14),
            orientation = "h",
            y = -0.25,
            itemsizing = "constant"
          ),
          xaxis = list(zeroline = FALSE),
          yaxis = list(zeroline = FALSE),
          margin = list(b = 200)
        )
    } else {
      mean_col1 <- names(res)[3]
      mean_col2 <- names(res)[4]
      res <- res %>% mutate(
        avg_expr = (!!sym(mean_col1) +!!sym(mean_col2)) / 2,
        highlight = ifelse(abs(log2FC) >= 1, "FC >= 2", "FC < 2")
      )
      p <- ggplot(res,
                  aes(
                    x = log2(avg_expr + 1),
                    y = log2FC,
                    color = highlight,
                    text = paste0(
                      "Gene: ",
                      gene,
                      "\nTranscript: ",
                      transcript_id,
                      "\nlog2FC: ",
                      round(log2FC, 2)
                    )
                  )) +
        geom_point(size = 1.5, alpha = 0.7) +
        scale_color_manual(values = c("FC >= 2" = BERRY_SIG, "FC < 2" = GREY_NONSIG)) +
        geom_hline(
          yintercept = 0,
          color = "grey30",
          linewidth = 0.3
        ) +
        labs(
          title = paste("FC:", input$de_cl1, "vs", input$de_cl2),
          x = "log2(Mean TPM + 1)",
          y = "log2 Fold Change",
          color = NULL
        ) +
        theme_minimal(base_size = 13) + theme(legend.position = "bottom")
      ggplotly(p, tooltip = "text") %>%
        layout(
          legend = list(font = list(size = 13)),
          xaxis = list(zeroline = FALSE),
          yaxis = list(zeroline = FALSE)
        )
    }
  })
  
  output$de_table <- renderDT({
    res <- de_results()
    if (is.null(res) || nrow(res) == 0)
      return(NULL)
    display_res <- res %>% select(-any_of(c("sig_category", "hover_text")))
    display_res <- join_annotations(display_res)
    # Force transcript_id (Isoform ID) and gene (Gene name) to columns 1 and 2
    lead_cols <- intersect(c("transcript_id", "gene"), names(display_res))
    if (length(lead_cols) > 0) {
      other_cols <- setdiff(names(display_res), lead_cols)
      display_res <- display_res[, c(lead_cols, other_cols), drop = FALSE]
    }
    if (HAS_TX_ANNOT &&
        !is.null(input$de_biotype_filter) &&
        input$de_biotype_filter != "All")
      display_res <- display_res %>% filter(Biotype == input$de_biotype_filter)
    skip_cols <- c(
      "gene",
      "gene_name",
      "gene_id",
      "transcript_id",
      "RefSeq",
      "Biotype",
      "UniProt",
      "Pfam",
      "InterPro"
    )
    int_cols <- intersect(c("Exons", "Transcript Length (bp)", "CDS (bp)"),
                          names(display_res))
    round_cols <- setdiff(names(display_res), c(skip_cols, int_cols))
    datatable(
      display_res,
      selection = "single",
      options = list(pageLength = 20, scrollX = TRUE),
      rownames = FALSE
    ) %>%
      formatRound(columns = round_cols, digits = 4) %>%
      formatRound(columns = int_cols, digits = 0)
  })
  
  output$dl_de_table <- downloadHandler(
    filename = function() {
      prefix <- switch(de_active_method(),
                       deseq2_gene = "DESeq2_gene",
                       deseq2_tx = "DESeq2_tx",
                       "DE")
      if (input$de_mode == "tissue")
        paste0(prefix,
               "_",
               input$de_tissue1,
               "_vs_",
               input$de_tissue2,
               ".csv")
      else
        paste0("FC_", input$de_cl1, "_vs_", input$de_cl2, ".csv")
    },
    content = function(file) {
      res <- de_results()
      if (!is.null(res))
        write.csv(res %>% select(-any_of(c(
          "sig_category", "hover_text"
        ))), file, row.names = FALSE)
    }
  )
  
  output$dl_de_sig <- downloadHandler(
    filename = function() {
      prefix <- switch(de_active_method(),
                       deseq2_gene = "DESeq2_gene_sig",
                       deseq2_tx = "DESeq2_tx_sig",
                       "DE_sig")
      if (input$de_mode == "tissue")
        paste0(
          prefix,
          "_",
          input$de_tissue1,
          "_vs_",
          input$de_tissue2,
          "_fc",
          input$de_volcano_fc,
          "_p",
          input$de_volcano_pval,
          ".csv"
        )
      else
        paste0("FC_sig_", input$de_cl1, "_vs_", input$de_cl2, ".csv")
    },
    content = function(file) {
      res <- de_results()
      if (is.null(res))
        return()
      if ("padj" %in% names(res)) {
        sig <- res %>%
          filter(padj <= input$de_volcano_pval,
                 abs(log2FC) >= input$de_volcano_fc) %>%
          select(-any_of(c("sig_category", "hover_text")))
      } else {
        sig <- res %>%
          filter(abs(log2FC) >= input$de_volcano_fc) %>%
          select(-any_of(c("sig_category", "hover_text")))
      }
      write.csv(sig, file, row.names = FALSE)
    }
  )
  
  output$dl_de_plot <- downloadHandler(
    filename = function()
      paste0(
        "DE_volcano_",
        input$de_tissue1,
        "_vs_",
        input$de_tissue2,
        ".pdf"
      ),
    content = function(file) {
      res <- de_results()
      if (is.null(res) ||
          nrow(res) == 0 || !"padj" %in% names(res))
        return()
      fc_cut <- input$de_volcano_fc
      pval_cut <- input$de_volcano_pval
      res <- res %>% mutate(
        sig_category = case_when(
          padj < pval_cut & abs(log2FC) >= fc_cut ~ "Significant",
          padj < pval_cut ~ "Sig. p-value only",
          TRUE ~ "Not significant"
        )
      )
      color_map <- c(
        "Significant" = BERRY_SIG,
        "Sig. p-value only" = BERRY_LIGHT,
        "Not significant" = GREY_NONSIG
      )
      p <- ggplot(res, aes(
        x = log2FC,
        y = -log10(padj),
        color = sig_category
      )) +
        geom_point(size = 1.5, alpha = 0.7) + scale_color_manual(values = color_map) +
        geom_hline(
          yintercept = -log10(input$de_volcano_pval),
          linetype = "dashed",
          color = "grey50"
        ) +
        geom_vline(
          xintercept = c(-input$de_volcano_fc, input$de_volcano_fc),
          linetype = "dashed",
          color = "grey50"
        ) +
        labs(
          title = paste0("DE: ", input$de_tissue1, " vs ", input$de_tissue2),
          x = "log2 Fold Change",
          y = "-log10(Adj. p-value)",
          color = NULL
        ) +
        theme_minimal(base_size = 15) + theme(
          plot.title = element_text(size = 18, face = "bold"),
          axis.text = element_text(size = 13),
          axis.title = element_text(size = 14),
          legend.position = "bottom",
          legend.text = element_text(size = 14),
          legend.title = element_text(size = 15)
        )
      ggsave(file, p, width = 16, height = 7)
    }
  )
  
  # === KEGG ===
  
  kegg_results_up <- reactiveVal(NULL)
  kegg_results_down <- reactiveVal(NULL)
  kegg_gene_list <- reactiveVal(NULL)
  kegg_source_label <- reactiveVal("")
  kegg_tissue1_label <- reactiveVal("")
  kegg_tissue2_label <- reactiveVal("")
  
  output$kegg_source_info <- renderUI({
    lbl <- kegg_source_label()
    if (nchar(lbl) == 0)
      return(helpText("Run DE or Switching first."))
    tags$div(style = "background: #f0f4ff; padding: 10px; border-radius: 6px; margin-bottom: 10px;", tags$strong("Source: "), tags$span(lbl))
  })
  output$kegg_up_label <- renderText({
    lbl <- kegg_tissue1_label()
    if (nchar(lbl) > 0)
      paste0("Enriched in ", lbl)
    else
      ""
  })
  output$kegg_down_label <- renderText({
    lbl <- kegg_tissue2_label()
    if (nchar(lbl) > 0)
      paste0("Enriched in ", lbl)
    else
      ""
  })
  
  observeEvent(input$run_kegg, {
    withProgress(message = "Running KEGG enrichment...", {
      if (input$kegg_source == "de") {
        res <- de_results()
        if (is.null(res)) {
          showNotification("Run DE first.", type = "warning")
          return()
        }
        if (!"padj" %in% names(res)) {
          showNotification("KEGG from DE requires Tissue vs Tissue.", type = "warning")
          return()
        }
        sig <- res %>% filter(padj < input$kegg_pval_cutoff,
                              abs(log2FC) >= input$kegg_fc_cutoff)
        gene_col <- if ("gene" %in% names(sig))
          "gene"
        else if ("gene_name" %in% names(sig))
          "gene_name"
        else
          NULL
        if (is.null(gene_col)) {
          showNotification("No gene column.", type = "error")
          return()
        }
        genes_up <- sig %>% filter(log2FC > 0) %>%
          pull(!!sym(gene_col)) %>% unique()
        genes_down <- sig %>% filter(log2FC < 0) %>%
          pull(!!sym(gene_col)) %>% unique()
        kegg_tissue1_label(input$de_tissue1)
        kegg_tissue2_label(input$de_tissue2)
        method_label <- switch(
          de_active_method(),
          deseq2_gene = "DESeq2 gene",
          deseq2_tx = "DESeq2 transcript",
          "Wilcoxon"
        )
        kegg_source_label(
          paste0(
            method_label,
            " DE: ",
            length(genes_up),
            " up in ",
            input$de_tissue1,
            ", ",
            length(genes_down),
            " up in ",
            input$de_tissue2
          )
        )
        universe_gene_col <- if ("gene" %in% names(res))
          "gene"
        else
          "gene_name"
        universe_symbols <- unique(res[[universe_gene_col]])
        universe_symbols <- universe_symbols[!grepl("^ENSG", universe_symbols)]
      } else {
        sw <- tryCatch(
          switch_results(),
          error = function(e)
            NULL
        )
        if (is.null(sw) || nrow(sw) == 0) {
          showNotification("Run switching first.", type = "warning")
          return()
        }
        genes_up <- unique(sw$gene)
        genes_down <- character(0)
        kegg_tissue1_label(paste0(
          input$switch_tissue1,
          " vs ",
          input$switch_tissue2,
          " (switching)"
        ))
        kegg_tissue2_label("")
        kegg_source_label(paste0("Switching: ", length(genes_up), " genes"))
        universe_symbols <- unique(active_avg()$gene)
        universe_symbols <- universe_symbols[!grepl("^ENSG", universe_symbols)]
      }
      genes_up <- genes_up[!grepl("^ENSG", genes_up)]
      genes_down <- genes_down[!grepl("^ENSG", genes_down)]
      setProgress(value = 0.3, detail = "Mapping to Entrez IDs...")
      entrez_map <- symbols_to_entrez(unique(c(genes_up, genes_down)))
      universe_entrez <- symbols_to_entrez(universe_symbols)
      kegg_gene_list(setNames(names(entrez_map), entrez_map))
      entrez_up <- entrez_map[names(entrez_map) %in% genes_up]
      entrez_down <- entrez_map[names(entrez_map) %in% genes_down]
      setProgress(value = 0.5, detail = "KEGG enrichment (up)...")
      kegg_results_up(if (length(entrez_up) > 0)
        tryCatch(
          enrichKEGG(
            gene = entrez_up,
            organism = "hsa",
            universe = as.character(universe_entrez),
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.2
          ),
          error = function(e)
            NULL
        )
        else
          NULL)
      setProgress(value = 0.7, detail = "KEGG enrichment (down)...")
      kegg_results_down(if (length(entrez_down) > 0)
        tryCatch(
          enrichKEGG(
            gene = entrez_down,
            organism = "hsa",
            universe = as.character(universe_entrez),
            pvalueCutoff = 0.05,
            qvalueCutoff = 0.2
          ),
          error = function(e)
            NULL
        )
        else
          NULL)
      if (is.null(kegg_results_up()) &&
          is.null(kegg_results_down()))
        showNotification("No enriched KEGG pathways.", type = "warning")
    })
  })
  
  make_kegg_dot <- function(kegg, direction_label, top_n = 20) {
    if (is.null(kegg) || nrow(as.data.frame(kegg)) == 0)
      return(NULL)
    kdf <- as.data.frame(kegg) %>% arrange(p.adjust) %>% head(top_n) %>%
      mutate(
        GeneRatio_num = sapply(GeneRatio, function(x) {
          parts <- as.numeric(strsplit(x, "/")[[1]])
          parts[1] / parts[2]
        }),
        Description = factor(Description, levels = rev(Description))
      )
    p <- ggplot(
      kdf,
      aes(
        x = GeneRatio_num,
        y = Description,
        size = Count,
        color = p.adjust,
        text = paste0(
          "Pathway: ",
          Description,
          "\nRatio: ",
          GeneRatio,
          "\nCount: ",
          Count,
          "\nAdj. p: ",
          signif(p.adjust, 3)
        )
      )
    ) +
      geom_point() +
      scale_color_gradientn(colours = rev(TROPICAL)) +
      scale_size_continuous(range = c(3, 10)) +
      labs(
        title = direction_label,
        x = "Gene Ratio",
        y = NULL,
        color = "Adj. p",
        size = "Count"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 12),
        axis.text.x = element_text(size = 11),
        axis.title = element_text(size = 13)
      )
    ggplotly(p, tooltip = "text")
  }
  
  output$kegg_dot_up <- renderPlotly({
    make_kegg_dot(
      kegg_results_up(),
      paste0("Enriched in ", kegg_tissue1_label()),
      input$kegg_top_n
    )
  })
  output$kegg_dot_down <- renderPlotly({
    make_kegg_dot(
      kegg_results_down(),
      paste0("Enriched in ", kegg_tissue2_label()),
      input$kegg_top_n
    )
  })
  
  output$kegg_table <- renderDT({
    up_df <- if (!is.null(kegg_results_up()))
      as.data.frame(kegg_results_up()) %>% mutate(direction = kegg_tissue1_label())
    else
      NULL
    down_df <- if (!is.null(kegg_results_down()))
      as.data.frame(kegg_results_down()) %>% mutate(direction = kegg_tissue2_label())
    else
      NULL
    kdf <- bind_rows(up_df, down_df)
    if (is.null(kdf) || nrow(kdf) == 0)
      return(NULL)
    kdf <- kdf %>%
      select(direction,
             ID,
             Description,
             GeneRatio,
             BgRatio,
             pvalue,
             p.adjust,
             qvalue,
             Count) %>%
      arrange(p.adjust)
    datatable(
      kdf,
      selection = "single",
      options = list(pageLength = 20, scrollX = TRUE),
      rownames = FALSE,
      colnames = c(
        "Enriched In",
        "ID",
        "Pathway",
        "Gene Ratio",
        "BgRatio",
        "p-value",
        "Adj. p",
        "q-value",
        "Count"
      )
    ) %>%
      formatSignif(columns = c("pvalue", "p.adjust", "qvalue"),
                   digits = 3)
  })
  
  output$kegg_gene_table <- renderDT({
    up_df <- if (!is.null(kegg_results_up()))
      as.data.frame(kegg_results_up()) %>% mutate(direction = kegg_tissue1_label())
    else
      NULL
    down_df <- if (!is.null(kegg_results_down()))
      as.data.frame(kegg_results_down()) %>% mutate(direction = kegg_tissue2_label())
    else
      NULL
    kdf <- bind_rows(up_df, down_df)
    sel <- input$kegg_table_rows_selected
    if (is.null(kdf) || is.null(sel) || length(sel) == 0)
      return(NULL)
    kdf <- kdf %>% arrange(p.adjust)
    entrez_ids <- strsplit(kdf$geneID[sel], "/")[[1]]
    entrez_to_sym <- kegg_gene_list()
    symbols <- ifelse(entrez_ids %in% names(entrez_to_sym),
                      entrez_to_sym[entrez_ids],
                      entrez_ids)
    datatable(
      tibble(
        Entrez_ID = entrez_ids,
        Gene_Symbol = symbols,
        Enriched_In = kdf$direction[sel]
      ),
      options = list(pageLength = 50, scrollX = TRUE),
      rownames = FALSE
    )
  })
  
  output$dl_kegg_table <- downloadHandler(
    filename = "kegg_enrichment.csv",
    content = function(file) {
      up_df <- if (!is.null(kegg_results_up()))
        as.data.frame(kegg_results_up()) %>% mutate(direction = kegg_tissue1_label())
      else
        NULL
      down_df <- if (!is.null(kegg_results_down()))
        as.data.frame(kegg_results_down()) %>% mutate(direction = kegg_tissue2_label())
      else
        NULL
      kdf <- bind_rows(up_df, down_df)
      if (!is.null(kdf) && nrow(kdf) > 0)
        write.csv(kdf, file, row.names = FALSE)
    }
  )
  
  output$dl_kegg_plot <- downloadHandler(
    filename = "kegg_enrichment.pdf",
    content = function(file) {
      make_kegg_gg <- function(kegg, label, top_n = 20) {
        if (is.null(kegg) || nrow(as.data.frame(kegg)) == 0)
          return(NULL)
        kdf <- as.data.frame(kegg) %>% arrange(p.adjust) %>% head(top_n) %>%
          mutate(
            GeneRatio_num = sapply(GeneRatio, function(x) {
              parts <- as.numeric(strsplit(x, "/")[[1]])
              parts[1] / parts[2]
            }),
            Description = factor(Description, levels = rev(Description))
          )
        ggplot(kdf,
               aes(
                 x = GeneRatio_num,
                 y = Description,
                 size = Count,
                 color = p.adjust
               )) +
          geom_point() + scale_color_gradientn(colours = rev(TROPICAL)) +
          scale_size_continuous(range = c(3, 10)) +
          labs(
            title = label,
            x = "Gene Ratio",
            y = NULL,
            color = "Adj. p",
            size = "Count"
          ) +
          theme_minimal(base_size = 14) +
          theme(
            plot.title = element_text(size = 16, face = "bold"),
            axis.text = element_text(size = 13),
            axis.title = element_text(size = 14),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 13)
          )
      }
      p1 <- make_kegg_gg(
        kegg_results_up(),
        paste0("Enriched in ", kegg_tissue1_label()),
        input$kegg_top_n
      )
      p2 <- make_kegg_gg(
        kegg_results_down(),
        paste0("Enriched in ", kegg_tissue2_label()),
        input$kegg_top_n
      )
      plots <- Filter(Negate(is.null), list(p1, p2))
      if (length(plots) == 0)
        return()
      pdf(file, width = 16, height = 7 * length(plots))
      for (p in plots)
        print(p)
      dev.off()
    }
  )
  
  # === Drug Response ===
  
  drug_data <- eventReactive(input$run_drug, {
    ensure_db()
    done <- log_timing("Drug response")
    on.exit(done())
    req(input$drug_select)
    withProgress(message = "Computing drug response embedding...", {
      use_gene <- HAS_GENE_TPM &&
        !is.null(input$drug_level) && input$drug_level == "gene"
      use_3d <- !is.null(input$drug_3d) && input$drug_3d
      n_comps <- if (use_3d)
        3
      else
        2
      if (use_gene) {
        gmat <- gene_tpm_df
        sample_data_cols <- intersect(sample_cols, names(gmat))
        num_mat <- as.matrix(gmat[, sample_data_cols, drop = FALSE])
        num_mat[is.na(num_mat)] <- 0
        gene_vars <- apply(num_mat, 1, var, na.rm = TRUE)
        top_idx <- order(gene_vars, decreasing = TRUE)[1:min(input$drug_top_var, nrow(num_mat))]
        num_mat <- num_mat[top_idx, , drop = FALSE]
        col_order <- sample_data_cols
        tmat <- t(num_mat)
        level_label <- "Gene-Level"
      } else {
        top_tx <- dbGetQuery(
          db,
          "SELECT transcript_id FROM isoform_variance ORDER BY variance DESC LIMIT ?",
          params = list(input$drug_top_var)
        )$transcript_id
        if (length(top_tx) == 0)
          return(NULL)
        placeholders <- paste(rep("?", length(top_tx)), collapse = ",")
        mat_df <- fetch_quants(paste0("SELECT * FROM quants WHERE transcript_id IN (",
            placeholders,
            ")"
          ),
          params = as.list(top_tx)
        )
        col_order <- intersect(active_meta()$sample_label, names(mat_df))
        tmat <- t(as.matrix(mat_df[, col_order, drop = FALSE]))
        level_label <- "Isoform-Level"
      }
      tmat[is.na(tmat)] <- 0
      tmat <- log2(tmat + 1)
      sample_info <- active_meta() %>%
        select(sample_label, CellLine, Tissue) %>%
        rename(sample = sample_label)
      ccle_slabels_drug <- active_meta()$sample_label[active_meta()$Tissue %in% CCLE_TISSUES]
      col_order <- intersect(intersect(col_order, ccle_slabels_drug),
                             rownames(tmat))
      tmat_ccle <- tmat[col_order, , drop = FALSE]
      sample_info_ccle <- sample_info %>% filter(sample %in% col_order)
      gbm_tmat_drug <- NULL
      gbm_info_drug <- NULL
      if (!is.null(input$dataset_mode) &&
          active_id() != REF_ID) {
        gbm_slabels_drug <- active_meta()$sample_label[active_meta()$Tissue %in% query_tissues()]
        gbm_cols_drug <- intersect(gbm_slabels_drug, rownames(tmat))
        if (length(gbm_cols_drug) > 0) {
          gbm_tmat_drug <- tmat[gbm_cols_drug, , drop = FALSE]
          gbm_info_drug <- sample_info %>% filter(sample %in% gbm_cols_drug)
        }
      }
      tmat <- tmat_ccle
      sample_info <- sample_info_ccle
      if (length(input$drug_tissue_filter) > 0) {
        sample_info <- sample_info %>% filter(Tissue %in% input$drug_tissue_filter)
        col_order <- intersect(sample_info$sample, rownames(tmat))
        tmat <- tmat[col_order, , drop = FALSE]
      }
      need_n <- if (input$drug_reduction == "umap") 5 else 3
      if (nrow(tmat) < need_n) {
        showNotification(paste0("The embedding needs at least ", need_n, " cell lines for ",
                                toupper(input$drug_reduction), "; the selected tissues have ",
                                nrow(tmat), ". Add more tissues."),
                         type = "warning", duration = 8)
        return(NULL)
      }
      keep_feat <- apply(tmat, 2, function(v) isTRUE(var(v) > 0))
      tmat <- tmat[, keep_feat, drop = FALSE]
      if (!is.null(gbm_tmat_drug) && ncol(gbm_tmat_drug) == length(keep_feat)) gbm_tmat_drug <- gbm_tmat_drug[, keep_feat, drop = FALSE]
      setProgress(value = 0.3, detail = "Running embedding...")
      if (input$drug_reduction == "umap") {
        umap_drug_model <- umap(
          tmat,
          n_neighbors = min(15, nrow(tmat) - 1),
          min_dist = 0.1,
          n_components = n_comps,
          metric = "cosine",
          ret_model = TRUE
        )
        embed_df <- as.data.frame(umap_drug_model$embedding)
        names(embed_df)[1:2] <- c("Dim1", "Dim2")
        if (use_3d &&
            ncol(embed_df) >= 3)
          names(embed_df)[3] <- "Dim3"
        embed_df$sample <- rownames(tmat)
        ax1 <- "UMAP 1"
        ax2 <- "UMAP 2"
        ax3 <- if (use_3d)
          "UMAP 3"
        else
          NULL
      } else {
        pca <- prcomp(tmat, scale. = TRUE, center = TRUE)
        n_keep <- min(max(n_comps, 2), ncol(pca$x))
        embed_df <- as.data.frame(pca$x[, 1:n_keep])
        names(embed_df)[1:2] <- c("Dim1", "Dim2")
        if (use_3d &&
            ncol(embed_df) >= 3)
          names(embed_df)[3] <- "Dim3"
        embed_df$sample <- rownames(embed_df)
        ve <- summary(pca)$importance[2, 1:min(3, ncol(pca$x))] * 100
        ax1 <- paste0("PC1 (", round(ve[1], 1), "%)")
        ax2 <- paste0("PC2 (", round(ve[2], 1), "%)")
        ax3 <- if (use_3d &&
                   length(ve) >= 3)
          paste0("PC3 (", round(ve[3], 1), "%)")
        else
          NULL
      }
      if (!is.null(gbm_tmat_drug) && !is.null(gbm_info_drug)) {
        if (input$drug_reduction == "umap") {
          gbm_drug_proj <- as.data.frame(umap_transform(gbm_tmat_drug, umap_drug_model))
        } else {
          gbm_drug_scaled <- scale(gbm_tmat_drug,
                                   center = pca$center,
                                   scale = pca$scale)
          gbm_drug_proj <- as.data.frame(gbm_drug_scaled %*% pca$rotation[, 1:n_comps])
        }
        names(gbm_drug_proj)[1:2] <- c("Dim1", "Dim2")
        if (use_3d &&
            ncol(gbm_drug_proj) >= 3)
          names(gbm_drug_proj)[3] <- "Dim3"
        gbm_drug_proj$sample <- rownames(gbm_tmat_drug)
        embed_df <- rbind(embed_df, gbm_drug_proj)
        sample_info <- rbind(sample_info, gbm_info_drug)
      }
      setProgress(value = 0.6, detail = "Joining drug response...")
      active_gdsc <- if (input$drug_dataset == "gdsc1" &&
                         !is.null(gdsc1_response_df))
        gdsc1_response_df
      else
        gdsc_response_df
      drug_sub <- active_gdsc %>%
        filter(DRUG_NAME == input$drug_select) %>%
        select(any_of(c("CellLine", "LN_IC50", "AUC", "MAX_CONC", "PUTATIVE_TARGET", "PATHWAY_NAME"))) %>%
        distinct(CellLine, .keep_all = TRUE)
      df <- embed_df %>%
        left_join(sample_info, by = "sample") %>%
        left_join(drug_sub, by = "CellLine")
      med_auc <- median(df$AUC, na.rm = TRUE)
      tissue_medians <- df %>% filter(!is.na(AUC)) %>%
        group_by(Tissue) %>%
        summarise(tissue_med_auc = median(AUC, na.rm = TRUE),
                  .groups = "drop")
      df <- df %>% left_join(tissue_medians, by = "Tissue") %>%
        mutate(
          sensitivity = case_when(
            is.na(AUC) ~ "No data",
            AUC <= tissue_med_auc ~ "Sensitive",
            TRUE ~ "Resistant"
          ),
          sensitivity_global = case_when(
            is.na(AUC) ~ "No data",
            AUC <= med_auc ~ "Sensitive",
            TRUE ~ "Resistant"
          )
        )
      target <- if (nrow(drug_sub) > 0)
        drug_sub$PUTATIVE_TARGET[1]
      else
        "Unknown"
      pathway <- if (nrow(drug_sub) > 0)
        drug_sub$PATHWAY_NAME[1]
      else
        "Unknown"
      df$is_gbm <- df$Tissue %in% query_tissues()
      list(
        df = df,
        ax1 = ax1,
        ax2 = ax2,
        ax3 = if (exists("ax3"))
          ax3
        else
          NULL,
        is_3d = use_3d,
        level = level_label,
        drug = input$drug_select,
        med_auc = med_auc,
        target = target,
        pathway = pathway,
        method = input$drug_reduction,
        n_with_data = sum(!is.na(df$AUC)),
        dataset = toupper(input$drug_dataset)
      )
    })
  })
  
  output$drug_plot <- renderPlotly({
    res <- drug_data()
    if (is.null(res))
      return(NULL)
    df <- res$df
    subtitle <- paste0(
      res$drug,
      " (",
      res$target,
      ") \u2014 ",
      res$n_with_data,
      "/",
      nrow(df),
      " with data | ",
      res$dataset,
      " [",
      res$level,
      "]"
    )
    if (res$is_3d && "Dim3" %in% names(df)) {
      if (input$drug_metric == "auc") {
        df_has <- df %>% filter(!is.na(AUC))
        df_no <- df %>% filter(is.na(AUC))
        plot_ly() %>%
          add_trace(
            data = df_no,
            x = ~ Dim1,
            y = ~ Dim2,
            z = ~ Dim3,
            type = "scatter3d",
            mode = "markers",
            marker = list(
              color = "#bdbdbd",
              size = 3,
              opacity = 0.4
            ),
            name = "No data",
            text = ~ paste0("Sample: ", CellLine, "\nTissue: ", Tissue),
            hoverinfo = "text"
          ) %>%
          add_trace(
            data = df_has,
            x = ~ Dim1,
            y = ~ Dim2,
            z = ~ Dim3,
            type = "scatter3d",
            mode = "markers",
            marker = list(
              color = ~ AUC,
              colorscale = list(c(0, "#ffffd9"), c(0.5, "#41b6c4"), c(1, "#081d58")),
              size = 4,
              opacity = 0.8,
              showscale = TRUE,
              colorbar = list(title = "AUC", len = 0.6)
            ),
            showlegend = FALSE,
            text = ~ paste0(
              "Sample: ",
              CellLine,
              "\nTissue: ",
              Tissue,
              "\nAUC: ",
              round(AUC, 3)
            ),
            hoverinfo = "text"
          ) %>%
          layout(
            title = list(text = subtitle, font = list(size = 12)),
            scene = list(
              xaxis = list(title = res$ax1),
              yaxis = list(title = res$ax2),
              zaxis = list(title = res$ax3)
            )
          )
      } else {
        sens_col <- if (input$drug_metric == "binary_global")
          "sensitivity_global"
        else
          "sensitivity"
        df$sens_display <- df[[sens_col]]
        plot_ly(
          df,
          x = ~ Dim1,
          y = ~ Dim2,
          z = ~ Dim3,
          color = ~ sens_display,
          colors = c(
            "Sensitive" = "#225ea8",
            "Resistant" = "#c51b8a",
            "No data" = "#e8e8e8"
          ),
          type = "scatter3d",
          mode = "markers",
          marker = list(size = 4, opacity = 0.8),
          text = ~ paste0(
            "Sample: ",
            CellLine,
            "\nTissue: ",
            Tissue,
            ifelse(is.na(AUC), "\nNo data", paste0("\nAUC: ", round(AUC, 3)))
          ),
          hoverinfo = "text"
        ) %>%
          layout(
            title = list(text = subtitle, font = list(size = 12)),
            scene = list(
              xaxis = list(title = res$ax1),
              yaxis = list(title = res$ax2),
              zaxis = list(title = res$ax3)
            )
          )
      }
    } else if (input$drug_metric == "auc") {
      df_has <- df %>% filter(!is.na(AUC))
      df_no <- df %>% filter(is.na(AUC))
      p <- plot_ly() %>%
        add_trace(
          data = df_no,
          x = ~ Dim1,
          y = ~ Dim2,
          type = "scatter",
          mode = "markers",
          marker = list(
            color = "#bdbdbd",
            size = input$drug_legend_size,
            opacity = 0.4
          ),
          name = "No data",
          text = ~ paste0("Sample: ", CellLine, "\nTissue: ", Tissue),
          hoverinfo = "text"
        ) %>%
        add_trace(
          data = df_has,
          x = ~ Dim1,
          y = ~ Dim2,
          type = "scatter",
          mode = "markers",
          marker = list(
            color = ~ AUC,
            colorscale = list(c(0, "#ffffd9"), c(0.5, "#41b6c4"), c(1, "#081d58")),
            size = input$drug_legend_size + 1,
            opacity = 0.8,
            showscale = TRUE,
            colorbar = list(
              title = "AUC",
              len = 0.6,
              y = 0.5
            )
          ),
          showlegend = FALSE,
          text = ~ paste0(
            "Sample: ",
            CellLine,
            "\nTissue: ",
            Tissue,
            "\nAUC: ",
            round(AUC, 3)
          ),
          hoverinfo = "text"
        )
      if ("is_gbm" %in% names(df) && any(df$is_gbm)) {
        gbm_auc_df <- df[df$is_gbm, ]
        p <- p %>% add_trace(
          data = gbm_auc_df,
          x = ~ Dim1,
          y = ~ Dim2,
          type = "scatter",
          mode = "markers",
          marker = list(
            color = "#c51b8a",
            symbol = "diamond",
            size = input$drug_legend_size + 2,
            opacity = 0.95
          ),
          name = paste0(DATASETS$display_name[DATASETS$dataset_id == active_id()], " sample"),
          text = ~ paste0(
            "Patient: ",
            CellLine,
            "\nCondition: ",
            Tissue,
            "\nNo drug data"
          ),
          hoverinfo = "text"
        )
      }
      p <- p %>%
        layout(
          title = list(text = subtitle, font = list(size = 12)),
          legend = list(font = list(size = 14), itemsizing = "constant"),
          xaxis = list(
            title = res$ax1,
            zeroline = FALSE,
            tickfont = list(size = 14),
            itemsizing = "constant"
          ),
          yaxis = list(
            title = res$ax2,
            zeroline = FALSE,
            tickfont = list(size = 14),
            itemsizing = "constant"
          )
        )
    } else {
      sens_col <- if (input$drug_metric == "binary_global")
        "sensitivity_global"
      else
        "sensitivity"
      df$sens_display <- df[[sens_col]]
      p <- ggplot(df,
                  aes(
                    x = Dim1,
                    y = Dim2,
                    color = sens_display,
                    shape = sens_display,
                    text = paste0(
                      "Sample: ",
                      CellLine,
                      "\nTissue: ",
                      Tissue,
                      ifelse(is.na(AUC), "\nNo data", paste0("\nAUC: ", round(AUC, 3)))
                    )
                  )) +
        geom_point(size = input$drug_legend_size * 0.4,
                   alpha = 0.8) +
        scale_color_manual(values = c(
          "Sensitive" = "#225ea8",
          "Resistant" = "#c51b8a",
          "No data" = "#e8e8e8"
        )) +
        scale_shape_manual(
          values = c(
            "Sensitive" = 16,
            "Resistant" = 17,
            "No data" = 4
          ),
          name = "Response"
        ) +
        labs(
          title = subtitle,
          x = res$ax1,
          y = res$ax2,
          color = "Response",
          shape = "Response"
        ) +
        theme_minimal() +
        guides(
          color = guide_legend(override.aes = list(size = input$drug_legend_size)),
          shape = guide_legend(override.aes = list(size = input$drug_legend_size))
        )
      if (!is.null(input$drug_split) &&
          input$drug_split && !res$is_3d) {
        p <- p %+% (df %>% dplyr::filter(sens_display != "No data")) + facet_wrap(~
                                                                                    sens_display, ncol = 2)
      }
      if ("is_gbm" %in% names(df) && any(df$is_gbm)) {
        gbm_drug_df <- df[df$is_gbm, ]
        p <- p + geom_point(
          data = gbm_drug_df,
          aes(
            x = Dim1,
            y = Dim2,
            text = paste0(
              "Patient: ",
              CellLine,
              "\nCondition: ",
              Tissue,
              ifelse(is.na(AUC), "\nNo drug data", paste0("\nAUC: ", round(AUC, 3)))
            )
          ),
          color = "#c51b8a",
          shape = 18,
          size = 3.5,
          alpha = 0.95,
          inherit.aes = FALSE
        )
      }
      ggplotly(p, tooltip = "text") %>% toWebGL() %>%
        layout(
          legend = list(font = list(size = 14), itemsizing = "constant"),
          xaxis = list(zeroline = FALSE),
          yaxis = list(zeroline = FALSE)
        )
    }
  })
  
  output$drug_table <- renderDT({
    res <- drug_data()
    if (is.null(res))
      return(NULL)
    sens_col <- if (input$drug_metric == "binary_global")
      "sensitivity_global"
    else
      "sensitivity"
    df <- res$df %>% filter(!is.na(AUC)) %>%
      mutate(Response = dplyr::recode(.data[[sens_col]],
                                      Sensitive = "Below median AUC",
                                      Resistant = "Above median AUC"))
    has_max <- "MAX_CONC" %in% names(df)
    if (has_max)
      df <- df %>% mutate(IC50_beyond = ifelse(LN_IC50 > log(MAX_CONC), "Yes", "No"))
    df <- df %>%
      select(any_of(c("CellLine", "Tissue", "AUC", "LN_IC50", "MAX_CONC", "IC50_beyond", "tissue_med_auc", "Response"))) %>%
      arrange(AUC)
    nice <- c(CellLine = "Cell Line", Tissue = "Tissue", AUC = "AUC", LN_IC50 = "LN(IC50)",
              MAX_CONC = "Max tested conc (uM)", IC50_beyond = "IC50 above max conc",
              tissue_med_auc = "Tissue Median AUC", Response = "Relative response")
    datatable(
      df,
      options = list(pageLength = 20, scrollX = TRUE),
      rownames = FALSE,
      colnames = unname(nice[names(df)])
    ) %>%
      formatRound(columns = c("AUC", "LN_IC50", "tissue_med_auc"),
                  digits = 3) %>%
      formatStyle("Response", backgroundColor = styleEqual(c("Below median AUC", "Above median AUC"), c("#dce8f5", "#f0e0e0")))
  })
  
  output$dl_drug_plot <- downloadHandler(
    filename = function()
      paste0("drug_", input$drug_dataset, "_", input$drug_select, ".pdf"),
    content = function(file) {
      res <- drug_data()
      if (is.null(res))
        return()
      df <- res$df
      if (input$drug_metric == "auc") {
        p <- ggplot(df %>% filter(!is.na(AUC)), aes(x = Dim1, y = Dim2, color = AUC)) +
          geom_point(size = 2.5, alpha = 0.8) +
          scale_color_gradientn(colours = TROPICAL) +
          labs(
            title = paste0(res$drug, " (", res$target, ") -- ", res$dataset),
            x = res$ax1,
            y = res$ax2,
            color = "AUC"
          ) +
          theme_minimal()
      } else {
        sens_col <- if (input$drug_metric == "binary_global")
          "sensitivity_global"
        else
          "sensitivity"
        df$sens_display <- df[[sens_col]]
        p <- ggplot(df, aes(x = Dim1, y = Dim2, color = sens_display)) +
          geom_point(size = 2.5, alpha = 0.8) +
          scale_color_manual(values = c(
            "Sensitive" = "#225ea8",
            "Resistant" = "#c51b8a",
            "No data" = "#e8e8e8"
          )) +
          labs(
            title = paste0(res$drug, " -- ", res$dataset),
            x = res$ax1,
            y = res$ax2,
            color = "Response"
          ) +
          theme_minimal()
      }
      ggsave(file, p, width = 10, height = 8)
    }
  )
  
  output$dl_drug_table <- downloadHandler(
    filename = function()
      paste0(
        "drug_response_",
        input$drug_dataset,
        "_",
        input$drug_select,
        ".csv"
      ),
    content = function(file) {
      res <- drug_data()
      if (is.null(res))
        return()
      df <- res$df %>% filter(!is.na(AUC)) %>%
        mutate(
          response_tissue = dplyr::recode(sensitivity, Sensitive = "Below median AUC", Resistant = "Above median AUC"),
          response_global = dplyr::recode(sensitivity_global, Sensitive = "Below median AUC", Resistant = "Above median AUC")
        )
      if ("MAX_CONC" %in% names(df))
        df <- df %>% mutate(IC50_above_max_conc = ifelse(LN_IC50 > log(MAX_CONC), "Yes", "No"))
      write.csv(
        df %>%
          select(any_of(c("CellLine", "Tissue", "AUC", "LN_IC50", "MAX_CONC", "IC50_above_max_conc",
                          "tissue_med_auc", "response_tissue", "response_global"))) %>%
          arrange(AUC),
        file,
        row.names = FALSE
      )
    }
  )
  
  # === Circos ===
  
  circos_plot_fn <- function() {
    tissue_col <- input$circos_tissue
    window_bp <- input$circos_window_kb * 1000
    chr_names <- c(paste0("chr", 1:22), "chrX", "chrY")
    chr_lengths <- c(
      248956422,
      242193529,
      198295559,
      190214555,
      181538259,
      170805979,
      159345973,
      145138636,
      138394717,
      133797422,
      135086622,
      133275309,
      114364328,
      107043718,
      101991189,
      90338345,
      83257441,
      80373285,
      58617616,
      64444167,
      46709983,
      50818468,
      156040895,
      57227415
    )
    genome_df <- data.frame(
      chr = chr_names,
      start = 0L,
      end = chr_lengths,
      stringsAsFactors = FALSE
    )
    our_txids <- unique(active_avg()$transcript_id)
    coords <- transcript_coords_df %>%
      filter(transcript_id %in% our_txids, chr %in% chr_names) %>%
      left_join(active_avg() %>% select(transcript_id, gene, all_of(tissue_col)), by = "transcript_id")
    coords$tpm <- coords[[tissue_col]]
    if (!is.null(input$circos_tissue_specific) &&
        input$circos_tissue_specific) {
      tissue_cols_all <- intersect(names(active_avg()), active_tissues())
      spec_df <- active_avg() %>%
        mutate(
          target = .data[[tissue_col]],
          total = rowSums(across(all_of(
            tissue_cols_all
          )), na.rm = TRUE),
          frac = ifelse(total > 0, target / total, 0)
        ) %>%
        filter(frac >= 0.5, target >= 1)
      coords <- coords %>% filter(transcript_id %in% spec_df$transcript_id)
    }
    coords$mid <- (coords$start + coords$end) / 2
    coords$bin_start <- floor(coords$mid / window_bp) * window_bp
    bin_expr <- coords %>%
      group_by(chr, bin_start) %>%
      summarise(
        avg_tpm = mean(tpm, na.rm = TRUE),
        n_tx = n(),
        .groups = "drop"
      ) %>%
      mutate(bin_end = bin_start + window_bp,
             avg_tpm_log = log2(avg_tpm + 1))
    bin_expr <- as.data.frame(bin_expr)
    max_val <- quantile(bin_expr$avg_tpm_log[bin_expr$avg_tpm_log > 0], 0.95, na.rm = TRUE)
    if (is.na(max_val) || max_val == 0)
      max_val <- 1
    bin_expr$avg_tpm_log <- pmin(bin_expr$avg_tpm_log, max_val)
    switch_data <- NULL
    if (input$circos_show_switches || input$circos_show_arcs) {
      sw <- tryCatch(
        switch_results(),
        error = function(e)
          NULL
      )
      if (is.null(sw) || nrow(sw) == 0)
        showNotification(
          "No switching data. Run Detect Switches first.",
          type = "warning",
          duration = 8
        )
      else {
        sw_top <- sw %>% arrange(desc(switch_score)) %>% head(input$circos_top_switch)
        if (!"dominant_tx1" %in% names(sw_top)) {
          sw_top$dominant_tx1 <- ""
          sw_top$dominant_tx2 <- ""
        }
        switch_data <- coords %>%
          filter(gene %in% sw_top$gene) %>%
          group_by(gene) %>%
          summarise(chr = chr[1],
                    mid = mean(mid),
                    .groups = "drop") %>%
          inner_join(sw_top %>% select(gene, switch_score, dominant_tx1, dominant_tx2),
                     by = "gene") %>%
          filter(chr %in% chr_names)
        switch_data <- as.data.frame(switch_data)
      }
    }
    par(mar = c(1, 1, 1, 1))
    circos.clear()
    circos.par(
      start.degree = 90,
      gap.after = c(rep(1, 23), 8),
      cell.padding = c(0, 0, 0, 0),
      track.margin = c(0.01, 0.01)
    )
    circos.genomicInitialize(genome_df,
                             plotType = c("axis", "labels"),
                             labels.cex = 0.55)
    expr_bed <- bin_expr[, c("chr", "bin_start", "bin_end", "avg_tpm_log")]
    circos.genomicTrack(
      expr_bed,
      ylim = c(0, max_val * 1.05),
      track.height = 0.32,
      bg.col = "#f8f9fa",
      bg.border = "#dee2e6",
      panel.fun = function(region, value, ...) {
        v <- value[[1]]
        col <- {
          pal <- switch(
            input$circos_palette,
            viridis = c("#440154", "#21918c", "#fde725"),
            heat = c("#bd0026", "#fd8d3c", "#ffffb2"),
            cool = c("#08306b", "#4292c6", "#f7fbff"),
            c("#225ea8", "#41b6c4", "#c51b8a")
          )
          hi_thresh <- input$circos_high_thresh / 100
          lo_thresh <- input$circos_low_thresh / 100
          ifelse(v > max_val * hi_thresh,
                 pal[3],
                 ifelse(v > max_val * lo_thresh, pal[2], pal[1]))
        }
        circos.genomicRect(
          region,
          value,
          ytop.column = 1,
          ybottom = 0,
          col = col,
          border = NA
        )
      }
    )
    if (!is.null(switch_data) &&
        nrow(switch_data) > 0 && input$circos_show_switches) {
      sw_bed <- data.frame(
        chr = switch_data$chr,
        start = as.integer(switch_data$mid),
        end = as.integer(switch_data$mid + 1),
        score = switch_data$switch_score,
        stringsAsFactors = FALSE
      )
      max_score <- max(sw_bed$score, na.rm = TRUE)
      sw_ramp <- colorRampPalette(c("#225ea8", "#FF4500", "#c51b8a"))
      sw_t1 <- tryCatch(
        input$switch_tissue1,
        error = function(e)
          "Tissue1"
      )
      sw_t2 <- tryCatch(
        input$switch_tissue2,
        error = function(e)
          "Tissue2"
      )
      sw_col1 <- TISSUE_COLORS[sw_t1]
      sw_col2 <- TISSUE_COLORS[sw_t2]
      if (is.na(sw_col1))
        sw_col1 <- "#225ea8"
      if (is.na(sw_col2))
        sw_col2 <- "#c51b8a"
      circos.genomicTrack(
        sw_bed[, c("chr", "start", "end", "score")],
        ylim = c(0, max_score * 1.1),
        track.height = 0.15,
        bg.col = "#fafafa",
        bg.border = "#dee2e6",
        panel.fun = function(region, value, ...) {
          s <- value[[1]]
          sn <- s / max_score
          pt_col <- ifelse(sn > 0.5, sw_col1, sw_col2)
          circos.genomicPoints(
            region,
            value,
            pch = 16,
            cex = 0.6 + sn * 1.2,
            col = pt_col
          )
        }
      )
    }
    if (!is.null(switch_data) &&
        nrow(switch_data) >= 2 && input$circos_show_arcs) {
      pairs <- switch_data[order(-switch_data$switch_score), ]
      np <- min(nrow(pairs), 25)
      if (np >= 2) {
        arc_ramp <- colorRampPalette(c("#c7e9b4", "#225ea8", "#c51b8a"))
        arc_cols <- arc_ramp(np - 1)
        for (i in seq_len(np - 1)) {
          score_avg <- (pairs$switch_score[i] + pairs$switch_score[i + 1]) / 2
          lwd_val <- 0.5 + (score_avg / max(pairs$switch_score)) * 3
          circos.link(
            pairs$chr[i],
            as.integer(pairs$mid[i]),
            pairs$chr[i + 1],
            as.integer(pairs$mid[i + 1]),
            col = adjustcolor(arc_cols[i], alpha.f = 0.5),
            lwd = lwd_val,
            rou = 0.3
          )
        }
      }
    }
    text(
      0,
      0,
      paste0("Isoform\nExpression\n", tissue_col),
      cex = 0.9,
      col = "#1a4d8f",
      font = 2
    )
    
    pal_legend <- switch(
      input$circos_palette,
      viridis = c("#fde725", "#21918c", "#440154"),
      heat = c("#ffffb2", "#fd8d3c", "#bd0026"),
      cool = c("#f7fbff", "#4292c6", "#08306b"),
      c("#c51b8a", "#41b6c4", "#225ea8")
    )
    legend(
      "bottomleft",
      legend = c("High expression", "Moderate", "Low expression"),
      fill = pal_legend,
      border = NA,
      cex = 0.9,
      title = "Expression Level",
      title.cex = 1.0,
      bg = "white",
      box.col = "grey80"
    )
    if (!is.null(switch_data) &&
        nrow(switch_data) > 0 && input$circos_show_switches) {
      sw_t1_lbl <- tryCatch(
        input$switch_tissue1,
        error = function(e)
          "Tissue 1"
      )
      sw_t2_lbl <- tryCatch(
        input$switch_tissue2,
        error = function(e)
          "Tissue 2"
      )
      sw_col1_leg <- TISSUE_COLORS[sw_t1_lbl]
      sw_col2_leg <- TISSUE_COLORS[sw_t2_lbl]
      if (is.na(sw_col1_leg))
        sw_col1_leg <- "#225ea8"
      if (is.na(sw_col2_leg))
        sw_col2_leg <- "#c51b8a"
      legend(
        "bottomright",
        legend = c(
          paste0("Favors ", sw_t1_lbl),
          paste0("Favors ", sw_t2_lbl)
        ),
        pch = 16,
        col = c(sw_col1_leg, sw_col2_leg),
        pt.cex = c(1.5, 1.5),
        cex = 0.9,
        title = paste0("Switches: ", sw_t1_lbl, " vs ", sw_t2_lbl),
        title.cex = 1.0,
        bg = "white",
        box.col = "grey80"
      )
    }
    circos.clear()
  }
  
  output$circos_plot <- renderPlot({
    input$run_circos
    isolate(circos_plot_fn())
  }, bg = "white")
  output$dl_circos <- downloadHandler(
    filename = function()
      paste0("circos_", input$circos_tissue, ".pdf"),
    content = function(file) {
      pdf(file, width = 10, height = 10)
      circos_plot_fn()
      dev.off()
    }
  )
  
  output$circos_chr_plot <- renderPlotly({
    validate(need(input$circos_chr_detail != "None",
                  "Pick a chromosome under 'Chromosome Detail' in the sidebar to see this view."))
    tissue_col <- input$circos_tissue
    sel_chr <- input$circos_chr_detail
    window_bp <- input$circos_window_kb * 1000
    our_txids <- unique(active_avg()$transcript_id)
    norm_chr <- function(x) ifelse(grepl("^chr", x), x, paste0("chr", x))
    coords <- transcript_coords_df %>%
      filter(transcript_id %in% our_txids, norm_chr(chr) == sel_chr) %>%
      left_join(active_avg() %>% select(transcript_id, gene, all_of(tissue_col)), by = "transcript_id")
    coords$tpm <- coords[[tissue_col]]
    coords$mid <- (coords$start + coords$end) / 2
    coords$bin_start <- floor(coords$mid / window_bp) * window_bp
    bin_data <- coords %>%
      group_by(bin_start) %>%
      summarise(
        avg_tpm = mean(tpm, na.rm = TRUE),
        n_tx = n(),
        top_genes = paste0(head(unique(gene), 5), collapse = ", "),
        top_isoforms = paste0(head(unique(transcript_id), 5), collapse = ", "),
        .groups = "drop"
      ) %>%
      mutate(
        bin_end = bin_start + window_bp,
        bin_mid = bin_start + window_bp / 2,
        pos_mb = round(bin_mid / 1e6, 2)
      )
    bin_data <- as.data.frame(bin_data)
    validate(need(nrow(bin_data) > 0,
                  paste0("No expressed transcripts found on ", sel_chr, " for this tissue or group.")))
    bin_data$avg_tpm[!is.finite(bin_data$avg_tpm)] <- 0
    bin_data$hover <- paste0(
      sel_chr, ": ", format(bin_data$bin_start, big.mark = ",", scientific = FALSE, trim = TRUE),
      "-", format(bin_data$bin_end, big.mark = ",", scientific = FALSE, trim = TRUE),
      "<br>Avg TPM: ", round(bin_data$avg_tpm, 2),
      "<br>Transcripts: ", bin_data$n_tx,
      "<br>Genes: ", bin_data$top_genes,
      "<br>Isoforms: ", bin_data$top_isoforms
    )
    p <- plot_ly(
      bin_data, x = ~pos_mb, y = ~avg_tpm, type = "bar",
      width = window_bp / 1e6 * 0.9,
      marker = list(color = "#225ea8"),
      text = ~hover, hoverinfo = "text", textposition = "none",
      name = "Avg TPM"
    )
    sw <- tryCatch(switch_results(), error = function(e) NULL)
    if (!is.null(sw) && nrow(sw) > 0 && isTRUE(input$circos_show_switches)) {
      sw_coords <- coords %>%
        filter(gene %in% sw$gene) %>%
        group_by(gene) %>%
        summarise(mid = mean(mid), .groups = "drop") %>%
        inner_join(sw %>% select(gene, switch_score), by = "gene")
      if (nrow(sw_coords) > 0) {
        sw_coords <- as.data.frame(sw_coords)
        sw_coords$pos_mb <- sw_coords$mid / 1e6
        sw_coords$y <- max(bin_data$avg_tpm) * 1.05
        p <- p %>% add_trace(
          data = sw_coords, x = ~pos_mb, y = ~y, type = "scatter", mode = "markers",
          marker = list(color = "#c51b8a", opacity = 0.7,
                        size = 6 + 10 * sw_coords$switch_score / max(sw_coords$switch_score)),
          text = ~paste0("Switch: ", gene, "<br>Score: ", round(switch_score, 3)),
          hoverinfo = "text", name = "Switching genes", inherit = FALSE
        )
      }
    }
    p %>% layout(
      title = list(text = paste0(sel_chr, " (", legend_label(tissue_col), ")")),
      xaxis = list(title = "Position (Mb)"),
      yaxis = list(title = "Avg TPM"),
      bargap = 0.05,
      showlegend = FALSE
    )
  })
  
  # === Circos supplementary tables ===
  
  circos_tab10a_data <- reactive({
    sw <- tryCatch(
      switch_results(),
      error = function(e)
        NULL
    )
    if (is.null(sw) || nrow(sw) == 0)
      return(NULL)
    coords_sub <- transcript_coords_df %>%
      group_by(gene_name) %>%
      summarise(
        chr = chr[1],
        position_mb = round(mean((start + end) / 2) / 1e6, 2),
        .groups = "drop"
      )
    sw %>%
      arrange(desc(switch_score)) %>%
      head(input$circos_top_switch) %>%
      left_join(coords_sub, by = c("gene" = "gene_name")) %>%
      select(any_of(
        c(
          "gene",
          "chr",
          "position_mb",
          "switch_score",
          "confidence",
          "dominant_tx1",
          "prop1",
          "dominant_tx2",
          "prop2"
        )
      ))
  })
  
  output$circos_tab10a <- renderDT({
    df <- circos_tab10a_data()
    if (is.null(df))
      return(datatable(
        data.frame(Note = "Run 'Detect Switches' in the Isoform Switching tab first."),
        options = list(dom = 't'),
        rownames = FALSE
      ))
    datatable(df,
              options = list(pageLength = 15, scrollX = TRUE),
              rownames = FALSE) %>%
      formatRound(columns = intersect(
        c("switch_score", "prop1", "prop2", "position_mb"),
        names(df)
      ), digits = 3)
  })
  
  output$dl_circos_tab10a <- downloadHandler(
    filename = function() {
      t1 <- tryCatch(
        input$switch_tissue1,
        error = function(e)
          "tissue1"
      )
      t2 <- tryCatch(
        input$switch_tissue2,
        error = function(e)
          "tissue2"
      )
      paste0("top_switching_genes_", t1, "_vs_", t2, ".csv")
    },
    content = function(file) {
      df <- circos_tab10a_data()
      if (!is.null(df))
        write.csv(df, file, row.names = FALSE)
    }
  )
  
  circos_tab10b_data <- reactive({
    tissue_col <- input$circos_tissue
    chr_order <- c(paste0("chr", 1:22), "chrX", "chrY")
    expr <- transcript_coords_df %>%
      left_join(active_avg() %>% select(transcript_id, all_of(tissue_col)), by = "transcript_id") %>%
      filter(chr %in% chr_order)
    expr$tpm <- expr[[tissue_col]]
    base <- expr %>%
      group_by(chr) %>%
      summarise(
        avg_tpm = round(mean(tpm, na.rm = TRUE), 2),
        transcript_count = n(),
        .groups = "drop"
      )
    sw <- tryCatch(
      switch_results(),
      error = function(e)
        NULL
    )
    sw_counts <- if (!is.null(sw) && nrow(sw) > 0) {
      sw %>% left_join(transcript_coords_df %>% distinct(gene_name, chr),
                       by = c("gene" = "gene_name")) %>%
        count(chr, name = "switching_events")
    } else {
      data.frame(chr = character(0), switching_events = integer(0))
    }
    base %>%
      left_join(sw_counts, by = "chr") %>%
      mutate(switching_events = ifelse(is.na(switching_events), 0, switching_events),
             tissue = tissue_col) %>%
      select(tissue, chr, avg_tpm, transcript_count, switching_events) %>%
      arrange(match(chr, chr_order))
  })
  
  output$circos_tab10b <- renderDT({
    df <- circos_tab10b_data()
    datatable(
      df,
      options = list(
        pageLength = 24,
        scrollX = TRUE,
        dom = 't'
      ),
      rownames = FALSE
    )
  })
  
  output$dl_circos_tab10b <- downloadHandler(
    filename = function()
      paste0("per_chromosome_summary_", input$circos_tissue, ".csv"),
    content = function(file) {
      df <- circos_tab10b_data()
      if (!is.null(df))
        write.csv(df, file, row.names = FALSE)
    }
  )
}

enableBookmarking(store = "url")

shinyApp(ui = ui, server = server)
