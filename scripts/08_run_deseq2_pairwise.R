library(DESeq2)
library(parallel)
library(RSQLite)
library(DBI)

metadata_path       <- "data/sample_metadata.csv"
tx2gene_path        <- "shiny_data/tx2gene.csv"
ensg_map_path       <- "data/ensg_to_symbol.csv"
gene_counts_path    <- "data/star_counts_gene.csv"
tx_counts_path      <- "data/salmon_counts_transcript.csv"
db_out_path         <- "shiny_data/deseq2.db"
n_cores             <- max(1, detectCores() - 1)

sample_metadata <- read.csv(metadata_path)
tx2gene <- read.csv(tx2gene_path)
tx2gene$transcript_id_clean <- sub("\\..*", "", tx2gene$transcript_id)
tx_symbol_map <- setNames(tx2gene$gene_name, tx2gene$transcript_id_clean)
ensg_map <- read.csv(ensg_map_path)
gene_symbol_map <- setNames(ensg_map$gene_name, ensg_map$gene_id)

tissues <- sort(unique(sample_metadata$Tissue))
tissue_pairs <- combn(tissues, 2, simplify = FALSE)

# --- Gene-level ---

gene_counts <- read.csv(gene_counts_path, check.names = FALSE)
gene_ids <- gene_counts$gene_id
gene_counts$gene_id <- NULL
matched_gene <- intersect(names(gene_counts), sample_metadata$SRR_ID)
gene_mat <- round(as.matrix(gene_counts[, matched_gene]))
rownames(gene_mat) <- gene_ids
mode(gene_mat) <- "integer"
rm(gene_counts); gc()
gene_mat <- gene_mat[rowSums(gene_mat >= 10) >= 3, ]

run_gene_pair <- function(pair) {
  tissue_a <- pair[1]; tissue_b <- pair[2]
  pair_samples <- intersect(sample_metadata$SRR_ID[sample_metadata$Tissue %in% c(tissue_a, tissue_b)], matched_gene)
  if (length(pair_samples) < 4) return(NULL)
  pair_counts <- gene_mat[, pair_samples]
  pair_counts <- pair_counts[rowSums(pair_counts >= 10) >= 3, ]
  tissue_factor <- relevel(factor(sample_metadata$Tissue[match(pair_samples, sample_metadata$SRR_ID)]), ref = tissue_b)
  tryCatch({
    dds <- DESeq(DESeqDataSetFromMatrix(pair_counts, data.frame(Tissue = tissue_factor, row.names = pair_samples), ~Tissue), quiet = TRUE)
    res <- as.data.frame(results(dds, contrast = c("Tissue", tissue_a, tissue_b)))
    res$feature_id <- rownames(res)
    res$tissue1 <- tissue_a; res$tissue2 <- tissue_b
    res[!is.na(res$padj), ]
  }, error = function(e) NULL)
}

cat("Gene-level:", nrow(gene_mat), "features |", length(tissue_pairs), "pairs |", n_cores, "cores\n")
gene_results <- mclapply(tissue_pairs, run_gene_pair, mc.cores = n_cores)
gene_results <- do.call(rbind, gene_results[!sapply(gene_results, is.null)])
gene_results$gene_id_clean <- sub("\\..*", "", gene_results$feature_id)
gene_results$gene_symbol <- gene_symbol_map[gene_results$gene_id_clean]
gene_results$gene_symbol[is.na(gene_results$gene_symbol)] <- gene_results$gene_id_clean[is.na(gene_results$gene_symbol)]
gene_results$level <- "gene"
cat("Gene-level:", nrow(gene_results), "results\n")
rm(gene_mat); gc()

# --- Transcript-level ---

tx_counts <- read.csv(tx_counts_path, check.names = FALSE)
tx_ids <- tx_counts$transcript_id
tx_counts$transcript_id <- NULL
matched_tx <- intersect(names(tx_counts), sample_metadata$SRR_ID)
tx_mat <- round(as.matrix(tx_counts[, matched_tx]))
rownames(tx_mat) <- tx_ids
mode(tx_mat) <- "integer"
rm(tx_counts); gc()
tx_mat <- tx_mat[rowSums(tx_mat >= 10) >= 3, ]

run_tx_pair <- function(pair) {
  tissue_a <- pair[1]; tissue_b <- pair[2]
  pair_samples <- intersect(sample_metadata$SRR_ID[sample_metadata$Tissue %in% c(tissue_a, tissue_b)], matched_tx)
  if (length(pair_samples) < 4) return(NULL)
  pair_counts <- tx_mat[, pair_samples]
  pair_counts <- pair_counts[rowSums(pair_counts >= 10) >= 3, ]
  tissue_factor <- relevel(factor(sample_metadata$Tissue[match(pair_samples, sample_metadata$SRR_ID)]), ref = tissue_b)
  tryCatch({
    dds <- DESeq(DESeqDataSetFromMatrix(pair_counts, data.frame(Tissue = tissue_factor, row.names = pair_samples), ~Tissue), quiet = TRUE)
    res <- as.data.frame(results(dds, contrast = c("Tissue", tissue_a, tissue_b)))
    res$feature_id <- rownames(res)
    res$tissue1 <- tissue_a; res$tissue2 <- tissue_b
    res[!is.na(res$padj), ]
  }, error = function(e) NULL)
}

cat("Transcript-level:", nrow(tx_mat), "features |", length(tissue_pairs), "pairs |", n_cores, "cores\n")
tx_results <- mclapply(tissue_pairs, run_tx_pair, mc.cores = n_cores)
tx_results <- do.call(rbind, tx_results[!sapply(tx_results, is.null)])
tx_results$gene_symbol <- tx_symbol_map[tx_results$feature_id]
tx_results$gene_symbol[is.na(tx_results$gene_symbol)] <- tx_results$feature_id[is.na(tx_results$gene_symbol)]
tx_results$level <- "transcript"
cat("Transcript-level:", nrow(tx_results), "results\n")

# --- Build SQLite database ---

if (file.exists(db_out_path)) file.remove(db_out_path)
db <- dbConnect(RSQLite::SQLite(), db_out_path)

dbWriteTable(db, "gene_results", gene_results, overwrite = TRUE)
dbExecute(db, "CREATE INDEX idx_gene_pair ON gene_results (tissue1, tissue2)")

dbWriteTable(db, "transcript_results", tx_results, overwrite = TRUE)
dbExecute(db, "CREATE INDEX idx_tx_pair ON transcript_results (tissue1, tissue2)")

dbDisconnect(db)
cat("Database:", db_out_path, "\n")
