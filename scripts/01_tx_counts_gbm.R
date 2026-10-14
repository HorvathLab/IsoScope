.libPaths(c("/home/azka.sajjad/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
library(data.table)

salmon_root <- "/data/lab/smart/software/varscope_reviewers/varscope_prjna682450"
outdir <- "out/gbm"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

meta <- fread("data/samples_gbm.csv")
files <- file.path(salmon_root, paste0(meta$sample_id, "_salmon"), "quant.sf", "quant.sf")
names(files) <- meta$sample_id
stopifnot(all(file.exists(files)))
message("quant.sf files: ", length(files))

counts <- NULL
tpms <- NULL
for (s in names(files)) {
  q <- fread(files[s], select = c("Name", "TPM", "NumReads"))
  if (is.null(counts)) {
    counts <- data.table(transcript_id = q$Name)
    tpms <- data.table(transcript_id = q$Name)
  }
  stopifnot(identical(counts$transcript_id, q$Name))
  counts[[s]] <- q$NumReads
  tpms[[s]] <- q$TPM
}

fwrite(counts, file.path(outdir, "salmon_counts_transcript.csv"))
fwrite(tpms, file.path(outdir, "salmon_tpm_transcript.csv.gz"))
message(sprintf("transcripts: %d x %d samples", nrow(counts), ncol(counts) - 1))
