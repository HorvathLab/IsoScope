.libPaths(c("/home/azka.sajjad/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
library(data.table)
library(RSQLite)
library(DBI)

N_BOOT <- 1000
set.seed(42)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) stop("usage: Rscript 07_build_database.R <dataset_id>")
id <- args[1]

cfg <- fread("data/datasets.csv")
cfg <- cfg[dataset_id == id]
if (nrow(cfg) != 1) stop("unknown dataset: ", id)
outdir <- file.path("out", id)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
cat("dataset:", id, "|", cfg$display_name, "\n")

meta <- fread(cfg$samples_file)
if (!"sample_label" %in% names(meta)) stop("samples file needs sample_label")
if (!"group_label" %in% names(meta)) {
  if ("Tissue" %in% names(meta)) meta[, group_label := Tissue] else stop("no group_label")
}
if (!"subject_id" %in% names(meta)) {
  if ("CellLine" %in% names(meta)) meta[, subject_id := CellLine] else stop("no subject_id")
}
if (!"sample_id" %in% names(meta)) {
  if ("SRR_ID" %in% names(meta)) meta[, sample_id := SRR_ID] else stop("no sample_id")
}
if (id == "ccle") meta <- meta[!grepl("G207", group_label)]
cat("samples:", nrow(meta), "| subjects:", uniqueN(meta$subject_id),
    "| groups:", uniqueN(meta$group_label), "\n")

tpm <- fread(cfg$tpm_file)
id_col <- names(tpm)[1]
setnames(tpm, id_col, "transcript_id")
keep <- intersect(names(tpm), meta$sample_id)
if (length(keep) != nrow(meta)) stop("sample mismatch: matrix has ", length(keep), " of ", nrow(meta))
tpm <- tpm[, c("transcript_id", keep), with = FALSE]
setnames(tpm, keep, meta$sample_label[match(keep, meta$sample_id)])
cat("transcripts in matrix:", nrow(tpm), "\n")

tx2gene <- fread("data/tx2gene.csv")
setnames(tx2gene, c("transcript_id", "gene"))
tpm <- merge(tx2gene, tpm, by = "transcript_id")
setcolorder(tpm, c("gene", "transcript_id"))

labels <- meta$sample_label
mat <- as.matrix(tpm[, ..labels])
nz <- rowSums(mat) > 0
tpm <- tpm[nz]; mat <- mat[nz, , drop = FALSE]
genes <- tpm$gene; txids <- tpm$transcript_id
cat("transcripts retained:", nrow(tpm), "| genes:", uniqueN(genes), "\n")

db_path <- file.path(outdir, "quants.db")
if (file.exists(db_path)) file.remove(db_path)
db <- dbConnect(RSQLite::SQLite(), db_path)

dbWriteTable(db, "quants", tpm, overwrite = TRUE)
dbExecute(db, "CREATE INDEX idx_gene ON quants (gene)")
dbExecute(db, "CREATE INDEX idx_transcript ON quants (transcript_id)")

dbWriteTable(db, "samples", meta, overwrite = TRUE)

dbWriteTable(db, "dataset_meta", data.frame(
  dataset_id = id, display_name = cfg$display_name,
  is_reference = cfg$is_reference,
  n_samples = nrow(meta), n_subjects = uniqueN(meta$subject_id),
  n_groups = uniqueN(meta$group_label),
  quantifier = cfg$quantifier, annotation = cfg$annotation,
  annotation_version = cfg$annotation_version,
  design_formula = cfg$design_formula,
  source_accession = cfg$source_accession,
  remapped_from = cfg$remapped_from,
  mapping_rate = cfg$mapping_rate,
  n_transcripts = nrow(tpm), n_genes = uniqueN(genes)
), overwrite = TRUE)

grp <- split(seq_along(labels), meta$group_label)
cat("computing group averages\n")
avg <- data.table(gene = genes, transcript_id = txids)
for (g in names(grp)) avg[[g]] <- rowMeans(mat[, grp[[g]], drop = FALSE])
fwrite(avg, file.path(outdir, "group_avg.csv"))

cat("computing qc stats\n")
qc <- data.table(
  sample_label = labels,
  total_tpm = colSums(mat),
  n_isoforms = colSums(mat > 0),
  n_genes = sapply(seq_along(labels), function(i) uniqueN(genes[mat[, i] > 0]))
)
qc[, iso_per_gene := n_isoforms / pmax(n_genes, 1)]
dbWriteTable(db, "qc_stats", qc, overwrite = TRUE)

cat("computing isoform variance\n")
dbWriteTable(db, "isoform_variance", data.table(
  gene = genes, transcript_id = txids,
  variance = matrixStats::rowVars(mat)
), overwrite = TRUE)
dbExecute(db, "CREATE INDEX idx_var ON isoform_variance (variance DESC)")
dbExecute(db, "CREATE INDEX idx_gene_var ON isoform_variance (gene)")

cat("computing dominant isoforms, bootstrap n =", N_BOOT, "\n")
avgm <- as.matrix(avg[, names(grp), with = FALSE])
gidx <- split(seq_len(nrow(avg)), genes)
dom <- vector("list", length(gidx) * length(grp)); k <- 0L
for (gi in seq_along(gidx)) {
  rows <- gidx[[gi]]
  gname <- names(gidx)[gi]
  for (g in names(grp)) {
    cols <- grp[[g]]
    v <- avgm[rows, g]
    tot <- sum(v)
    if (tot == 0) next
    props <- v / tot
    top <- which.max(props)
    ci_lo <- NA_real_; ci_hi <- NA_real_
    if (length(cols) >= 3) {
      sub <- mat[rows, cols, drop = FALSE]
      ridx <- matrix(sample(seq_along(cols), N_BOOT * length(cols), replace = TRUE), nrow = N_BOOT)
      bp <- apply(ridx, 1, function(r) {
        a <- rowMeans(sub[, r, drop = FALSE]); s <- sum(a)
        if (s == 0) 0 else a[top] / s
      })
      ci_lo <- quantile(bp, 0.025, names = FALSE)
      ci_hi <- quantile(bp, 0.975, names = FALSE)
    }
    k <- k + 1L
    dom[[k]] <- list(gene = gname, group_label = g, dominant_tx = txids[rows][top],
                     dominant_prop = props[top], ci_lower = ci_lo, ci_upper = ci_hi,
                     n_samples = length(cols), total_tpm = tot)
  }
  if (gi %% 5000 == 0) cat("  ", gi, "/", length(gidx), "genes\n")
}
dom <- rbindlist(dom[seq_len(k)])
fwrite(dom, file.path(outdir, "dominant_isoforms.csv"))

cat("computing diversity\n")
div <- vector("list", length(gidx))
for (gi in seq_along(gidx)) {
  rows <- gidx[[gi]]
  e <- sapply(names(grp), function(g) {
    v <- avgm[rows, g]; v <- v[v > 0]
    if (length(v) <= 1) return(0)
    p <- v / sum(v)
    -sum(p * log2(p)) / log2(length(v))
  })
  div[[gi]] <- list(gene = names(gidx)[gi], entropy_range = max(e) - min(e),
                    most_diverse_group = names(e)[which.max(e)], max_entropy = max(e),
                    least_diverse_group = names(e)[which.min(e)], min_entropy = min(e))
}
div <- rbindlist(div)
setorder(div, -entropy_range)
fwrite(div, file.path(outdir, "diversity_range.csv"))

dbDisconnect(db)
cat("done:", db_path, "\n")
