library(data.table)

star_dir <- "data/star_output"
gtf_path <- "data/gencode.v49.basic.annotation.gtf"
outpath  <- "shiny_data/tpm_matrix.txt.gz"

files <- list.files(star_dir, pattern = "ReadsPerGene.out.tab$", full.names = TRUE)
samples <- gsub("_ReadsPerGene.out.tab", "", basename(files))

gtf <- read.table(gtf_path, sep = "\t", comment.char = "#")
genes <- gtf[gtf[, 3] == "gene", ]
gene_ids <- sub('gene_id ', '', genes[, 9])
gene_ids <- unlist(lapply(gene_ids, function(x) strsplit(x, ';')[[1]][1]))
gene_lengths <- setNames(genes[, 5] - genes[, 4], gene_ids)

first <- read.table(files[1], skip = 4)
counts_matrix <- matrix(nrow = nrow(first), ncol = length(files))
rownames(counts_matrix) <- first[, 1]
colnames(counts_matrix) <- samples

for (i in seq_along(files)) {
  counts_matrix[, i] <- read.table(files[i], skip = 4)[, 2]
}

common <- intersect(rownames(counts_matrix), names(gene_lengths))
counts_matrix <- counts_matrix[common, ]
lengths <- gene_lengths[common]

tpm <- ((counts_matrix / lengths) * 1e6) / sum(counts_matrix / lengths)

fwrite(data.table(gene_id = rownames(tpm), tpm), outpath, sep = "\t")
message(sprintf("Done: %d genes x %d samples", nrow(tpm), ncol(tpm)))
