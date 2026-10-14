.libPaths(c("/home/azka.sajjad/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
library(tximport)
library(data.table)

salmon_root <- "/data/lab/smart/software/varscope_reviewers/varscope_prjna682450"
outdir <- "out/gbm"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

meta <- fread("data/samples_gbm.csv")
files <- file.path(salmon_root, paste0(meta$sample_id, "_salmon"), "quant.sf", "quant.sf")
names(files) <- meta$sample_id
stopifnot(all(file.exists(files)))
message(sprintf("quant.sf files: %d", length(files)))

tx2gene <- fread("data/tx2gene.csv")
setnames(tx2gene, c("TXNAME", "GENEID"))
message(sprintf("tx2gene: %d transcripts, %d genes",
                nrow(tx2gene), uniqueN(tx2gene$GENEID)))

txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene,
                countsFromAbundance = "lengthScaledTPM")

message(sprintf("gene matrix: %d genes x %d samples",
                nrow(txi$counts), ncol(txi$counts)))

counts_dt <- cbind(data.table(gene = rownames(txi$counts)), as.data.table(txi$counts))
fwrite(counts_dt, file.path(outdir, "salmon_gene_counts.csv"))

tpm_dt <- cbind(data.table(gene = rownames(txi$abundance)), as.data.table(txi$abundance))
fwrite(tpm_dt, file.path(outdir, "salmon_gene_tpm.csv.gz"))

colsums <- colSums(txi$abundance)
message(sprintf("TPM column sums: min %.0f, max %.0f (expect 1e6)",
                min(colsums), max(colsums)))
message("done")
