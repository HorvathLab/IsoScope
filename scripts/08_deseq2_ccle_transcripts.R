.libPaths(c("/home/azka.sajjad/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
library(DESeq2)
library(parallel)
library(RSQLite)
library(DBI)

metadata_path  <- "data/sample_metadata.csv"
tx2gene_path   <- "data/tx2gene.csv"
tx_counts_path <- "out/salmon_counts_transcript_versioned.csv"
db_out_path    <- "out/deseq2_ccle.db"
n_cores        <- max(1, detectCores() - 2)

meta <- read.csv(metadata_path)
meta$SRR_ID <- meta[[grep("SRR", names(meta), value = TRUE)[1]]]
meta$Tissue <- meta[[grep("Tissue", names(meta), ignore.case = TRUE, value = TRUE)[1]]]
meta <- meta[!grepl("G207", meta$Tissue), ]
cat("CCLE samples:", nrow(meta), "| tissues:", length(unique(meta$Tissue)), "\n")

tx2gene <- read.csv(tx2gene_path)
tx_symbol_map <- setNames(tx2gene$gene_name, tx2gene$transcript_id)

tissues <- sort(unique(meta$Tissue))
tissue_pairs <- combn(tissues, 2, simplify = FALSE)
cat("pairs:", length(tissue_pairs), "| cores:", n_cores, "\n")

x <- read.csv(tx_counts_path, check.names = FALSE)
ids <- x[["transcript_id"]]
x[["transcript_id"]] <- NULL
matched <- intersect(names(x), meta$SRR_ID)
mat <- round(as.matrix(x[, matched]))
rownames(mat) <- ids
mode(mat) <- "integer"
rm(x); gc()
mat <- mat[rowSums(mat >= 10) >= 3, ]
cat("features:", nrow(mat), "| samples:", length(matched), "\n")
cat("id format:", rownames(mat)[1], "\n")

run_pair <- function(pair) {
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
    res$group1 <- a; res$group2 <- b
    res[!is.na(res$padj), ]
  }, error = function(e) NULL)
}

tx_results <- mclapply(tissue_pairs, run_pair, mc.cores = n_cores)
tx_results <- do.call(rbind, tx_results[!sapply(tx_results, is.null)])
tx_results$gene_symbol <- tx_symbol_map[tx_results$feature_id]
tx_results$gene_symbol[is.na(tx_results$gene_symbol)] <- tx_results$feature_id[is.na(tx_results$gene_symbol)]
tx_results$level <- "transcript"
cat("transcript results:", nrow(tx_results), "\n")
cat("unmapped symbols:", sum(tx_results$gene_symbol == tx_results$feature_id), "\n")
rm(mat); gc()

db <- dbConnect(RSQLite::SQLite(), db_out_path)
dbExecute(db, "DROP TABLE IF EXISTS transcript_results")
dbWriteTable(db, "transcript_results", tx_results, overwrite = TRUE)
dbExecute(db, "CREATE INDEX idx_tx_pair ON transcript_results (group1, group2)")
dbExecute(db, "CREATE INDEX idx_tx_feature ON transcript_results (feature_id)")

dbExecute(db, "ALTER TABLE gene_results RENAME COLUMN tissue1 TO group1")
dbExecute(db, "ALTER TABLE gene_results RENAME COLUMN tissue2 TO group2")
dbExecute(db, "CREATE INDEX IF NOT EXISTS idx_gene_feature ON gene_results (feature_id)")
dbDisconnect(db)
cat("done:", db_out_path, "\n")
