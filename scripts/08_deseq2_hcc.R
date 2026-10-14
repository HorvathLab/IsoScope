.libPaths(c("/home/azka.sajjad/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
library(DESeq2)
library(RSQLite)
library(DBI)

EXCLUDE <- "SRR37349497"

meta <- read.csv("data/samples_hcc.csv")
meta <- meta[meta$sample_id != EXCLUDE, ]

tab <- table(meta$subject_id, meta$group_label)
paired <- rownames(tab)[tab[, "Tumor"] > 0 & tab[, "Adjacent"] > 0]
meta <- meta[meta$subject_id %in% paired, ]
cat("paired subjects:", length(paired), "| samples:", nrow(meta), "\n")
print(table(meta$subject_id, meta$group_label))

meta$subject_id  <- factor(meta$subject_id)
meta$group_label <- factor(meta$group_label, levels = c("Adjacent", "Tumor"))

run_level <- function(path, id_col, label) {
  x <- read.csv(path, check.names = FALSE)
  ids <- x[[id_col]]
  x[[id_col]] <- NULL
  m <- round(as.matrix(x[, meta$sample_id]))
  rownames(m) <- ids
  mode(m) <- "integer"
  m <- m[rowSums(m >= 10) >= 3, ]
  cat(label, "features:", nrow(m), "\n")
  dds <- DESeq(DESeqDataSetFromMatrix(m, meta, ~ subject_id + group_label))
  res <- as.data.frame(results(dds, contrast = c("group_label", "Tumor", "Adjacent")))
  res$feature_id <- rownames(res)
  res$group1 <- "Tumor"; res$group2 <- "Adjacent"; res$level <- label
  res <- res[!is.na(res$padj), ]
  cat(label, "results:", nrow(res), "| padj<0.05:", sum(res$padj < 0.05), "\n")
  res
}

gene_res <- run_level("out/hcc/salmon_gene_counts.csv", "gene", "gene")
gene_res$gene_symbol <- gene_res$feature_id

tx_res <- run_level("out/hcc/salmon_counts_transcript.csv", "transcript_id", "transcript")
tx2gene <- read.csv("data/tx2gene.csv")
map <- setNames(tx2gene$gene_name, tx2gene$transcript_id)
tx_res$gene_symbol <- map[tx_res$feature_id]
tx_res$gene_symbol[is.na(tx_res$gene_symbol)] <- tx_res$feature_id[is.na(tx_res$gene_symbol)]

db_path <- "out/hcc/deseq2_hcc.db"
if (file.exists(db_path)) file.remove(db_path)
db <- dbConnect(RSQLite::SQLite(), db_path)
dbWriteTable(db, "gene_results", gene_res, overwrite = TRUE)
dbWriteTable(db, "transcript_results", tx_res, overwrite = TRUE)
dbDisconnect(db)
cat("done:", db_path, "\n")
