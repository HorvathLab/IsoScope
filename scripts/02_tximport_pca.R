.libPaths(c("/home/azka.sajjad/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
library(tximport)
library(data.table)

salmon_root <- "/data/lab/smart/software/varscope_reviewers/varscope_prjna1497887"
tx2gene_path <- "data/tx2gene.csv"
outdir <- "out/pca"

dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

dirs <- list.files(salmon_root, pattern = "^SRR.*_salmon$", full.names = TRUE)
files <- file.path(dirs, "quant.sf", "quant.sf")
names(files) <- sub("_salmon$", "", basename(dirs))

exists_mask <- file.exists(files)
if (any(!exists_mask)) {
  message("missing quant.sf for: ", paste(names(files)[!exists_mask], collapse = ", "))
  files <- files[exists_mask]
}
message(sprintf("quant.sf files: %d", length(files)))

tx2gene <- fread(tx2gene_path)
setnames(tx2gene, c("TXNAME", "GENEID"))
message(sprintf("tx2gene: %d transcripts, %d genes",
                nrow(tx2gene), uniqueN(tx2gene$GENEID)))

first <- fread(files[1], nrows = 5)
message("quant.sf ids: ", paste(head(first$Name, 2), collapse = ", "))
message("tx2gene ids:  ", paste(head(tx2gene$TXNAME, 2), collapse = ", "))

txi <- tximport(files,
                type = "salmon",
                tx2gene = tx2gene,
                countsFromAbundance = "lengthScaledTPM")

message(sprintf("gene matrix: %d genes x %d samples",
                nrow(txi$counts), ncol(txi$counts)))

counts_dt <- data.table(gene = rownames(txi$counts))
counts_dt <- cbind(counts_dt, as.data.table(txi$counts))
fwrite(counts_dt, file.path(outdir, "salmon_gene_counts.csv"))

tpm_dt <- data.table(gene = rownames(txi$abundance))
tpm_dt <- cbind(tpm_dt, as.data.table(txi$abundance))
fwrite(tpm_dt, file.path(outdir, "salmon_gene_tpm.csv.gz"))

colsums <- colSums(txi$abundance)
message(sprintf("TPM column sums: min %.0f, max %.0f (expect 1e6)",
                min(colsums), max(colsums)))
message("done")
