.libPaths(c("/home/azka.sajjad/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
library(DESeq2)
library(RSQLite)
library(DBI)

meta <- read.csv("data/samples_pca.csv")
meta$group_label <- factor(meta$group_label)
cat("samples:", nrow(meta), "| donors:", length(unique(meta$subject_id)), "\n")
print(table(meta$group_label))

groups <- sort(levels(meta$group_label))
pairs <- combn(groups, 2, simplify = FALSE)
cat("contrasts:", length(pairs), "\n")

run_level <- function(path, id_col, label) {
  x <- read.csv(path, check.names = FALSE)
  ids <- x[[id_col]]
  x[[id_col]] <- NULL
  m <- round(as.matrix(x[, meta$sample_id]))
  rownames(m) <- ids
  mode(m) <- "integer"
  m <- m[rowSums(m >= 10) >= 3, ]
  cat(label, "features:", nrow(m), "| id:", rownames(m)[1], "\n")

  out <- list()
  for (p in pairs) {
    a <- p[1]; b <- p[2]
    s <- meta$sample_id[meta$group_label %in% c(a, b)]
    if (length(s) < 4) next
    cnt <- m[, s]
    cnt <- cnt[rowSums(cnt >= 10) >= 3, ]
    gf <- relevel(factor(meta$group_label[match(s, meta$sample_id)]), ref = b)
    r <- tryCatch({
      dds <- DESeq(DESeqDataSetFromMatrix(cnt, data.frame(group_label = gf, row.names = s), ~group_label), quiet = TRUE)
      res <- as.data.frame(results(dds, contrast = c("group_label", a, b)))
      res$feature_id <- rownames(res)
      res$group1 <- a; res$group2 <- b
      res[!is.na(res$padj), ]
    }, error = function(e) NULL)
    if (!is.null(r)) {
      cat("  ", a, "vs", b, ":", nrow(r), "results,", sum(r$padj < 0.05), "sig\n")
      out[[length(out) + 1]] <- r
    }
  }
  res <- do.call(rbind, out)
  res$level <- label
  res
}

gene_res <- run_level("out/pca/salmon_gene_counts.csv", "gene", "gene")
gene_res$gene_symbol <- gene_res$feature_id

tx_res <- run_level("out/pca/salmon_counts_transcript.csv", "transcript_id", "transcript")
tx2gene <- read.csv("data/tx2gene.csv")
map <- setNames(tx2gene$gene_name, tx2gene$transcript_id)
tx_res$gene_symbol <- map[tx_res$feature_id]
tx_res$gene_symbol[is.na(tx_res$gene_symbol)] <- tx_res$feature_id[is.na(tx_res$gene_symbol)]
cat("unmapped symbols:", sum(tx_res$gene_symbol == tx_res$feature_id), "\n")

db_path <- "out/pca/deseq2_pca.db"
if (file.exists(db_path)) file.remove(db_path)
db <- dbConnect(RSQLite::SQLite(), db_path)
dbWriteTable(db, "gene_results", gene_res, overwrite = TRUE)
dbWriteTable(db, "transcript_results", tx_res, overwrite = TRUE)
dbExecute(db, "CREATE INDEX idx_gene_pair ON gene_results (group1, group2)")
dbExecute(db, "CREATE INDEX idx_tx_pair ON transcript_results (group1, group2)")
dbDisconnect(db)
cat("done:", db_path, "\n")
