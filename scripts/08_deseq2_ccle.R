.libPaths(c("/home/azka.sajjad/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
library(DESeq2)
library(parallel)
library(RSQLite)
library(DBI)

metadata_path    <- "data/sample_metadata.csv"
tx2gene_path     <- "data/tx2gene.csv"
gene_counts_path <- "out/salmon_gene_counts.csv"
tx_counts_path   <- "/data/lab/AZ/cell_line_project/salmon_counts_transcript.csv"
db_out_path      <- "out/deseq2_ccle.db"
n_cores          <- max(1, detectCores() - 2)

meta <- read.csv(metadata_path)
srr_col    <- grep("SRR", names(meta), value = TRUE)[1]
tissue_col <- grep("Tissue", names(meta), ignore.case = TRUE, value = TRUE)[1]
cat("using columns:", srr_col, "/", tissue_col, "\n")
meta$SRR_ID <- meta[[srr_col]]
meta$Tissue <- meta[[tissue_col]]

meta <- meta[!grepl("G207", meta$Tissue), ]
cat("CCLE samples:", nrow(meta), "| tissues:", length(unique(meta$Tissue)), "\n")

tx2gene <- read.csv(tx2gene_path)
tx2gene$tx_clean <- sub("\\..*", "", tx2gene$transcript_id)
tx_symbol_map <- setNames(tx2gene$gene_name, tx2gene$tx_clean)

tissues <- sort(unique(meta$Tissue))
tissue_pairs <- combn(tissues, 2, simplify = FALSE)
cat("pairs:", length(tissue_pairs), "| cores:", n_cores, "\n")

run_pair <- function(pair, mat, matched) {
  a <- pair[1]; b <- pair[2]
  s <- intersect(meta$SRR_ID[meta$Tissue %in% c(a, b)], matched)
  if (length(s) < 4) return(NULL)
  cnt <- mat[, s]
  cnt <- cnt[rowSums(cnt >= 10) >= 3, ]
  if (nrow(cnt) < 10) return(NULL)
  tf <- relevel(factor(meta$Tissue[match(s, meta$SRR_ID)]), ref = b)
  tryCatch({
    dds <- DESeq(DESeqDataSetFromMatrix(cnt, data.frame(Tissue = tf, row.names = s), ~Tissue), quiet = TRUE)
    res <- as.data.frame(results(dds, contrast = c("Tissue", a, b)))
    res$feature_id <- rownames(res)
    res$tissue1 <- a; res$tissue2 <- b
    res[!is.na(res$padj), ]
  }, error = function(e) NULL)
}

load_matrix <- function(path, id_col) {
  x <- read.csv(path, check.names = FALSE)
  ids <- x[[id_col]]
  x[[id_col]] <- NULL
  matched <- intersect(names(x), meta$SRR_ID)
  m <- round(as.matrix(x[, matched]))
  rownames(m) <- ids
  mode(m) <- "integer"
  rm(x); gc()
  m <- m[rowSums(m >= 10) >= 3, ]
  list(mat = m, matched = matched)
}

cat("\n--- gene level ---\n")
g <- load_matrix(gene_counts_path, "gene")
cat("features:", nrow(g$mat), "| samples:", length(g$matched), "\n")
gene_results <- mclapply(tissue_pairs, run_pair, mat = g$mat, matched = g$matched, mc.cores = n_cores)
gene_results <- do.call(rbind, gene_results[!sapply(gene_results, is.null)])
gene_results$gene_symbol <- gene_results$feature_id
gene_results$level <- "gene"
cat("gene results:", nrow(gene_results), "\n")
rm(g); gc()

cat("\n--- transcript level ---\n")
t <- load_matrix(tx_counts_path, "transcript_id")
cat("features:", nrow(t$mat), "| samples:", length(t$matched), "\n")
tx_results <- mclapply(tissue_pairs, run_pair, mat = t$mat, matched = t$matched, mc.cores = n_cores)
tx_results <- do.call(rbind, tx_results[!sapply(tx_results, is.null)])
tx_results$gene_symbol <- tx_symbol_map[sub("\\..*", "", tx_results$feature_id)]
tx_results$gene_symbol[is.na(tx_results$gene_symbol)] <- tx_results$feature_id[is.na(tx_results$gene_symbol)]
tx_results$level <- "transcript"
cat("transcript results:", nrow(tx_results), "\n")
rm(t); gc()

if (file.exists(db_out_path)) file.remove(db_out_path)
db <- dbConnect(RSQLite::SQLite(), db_out_path)
dbWriteTable(db, "gene_results", gene_results, overwrite = TRUE)
dbExecute(db, "CREATE INDEX idx_gene_pair ON gene_results (tissue1, tissue2)")
dbWriteTable(db, "transcript_results", tx_results, overwrite = TRUE)
dbExecute(db, "CREATE INDEX idx_tx_pair ON transcript_results (tissue1, tissue2)")
dbDisconnect(db)
cat("\ndone:", db_out_path, "\n")
