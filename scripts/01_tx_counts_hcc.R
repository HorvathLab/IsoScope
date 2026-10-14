.libPaths(c("/home/azka.sajjad/R/x86_64-pc-linux-gnu-library/4.6", .libPaths()))
library(data.table)

salmon_root <- "/data/lab/smart/software/varscope_reviewers/varscope_prjna1427671"
outdir <- "out/hcc"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

dirs <- list.files(salmon_root, pattern = "^SRR.*_salmon$", full.names = TRUE)
files <- file.path(dirs, "quant.sf", "quant.sf")
names(files) <- sub("_salmon$", "", basename(dirs))
files <- files[file.exists(files)]
message("quant.sf files: ", length(files))

counts <- NULL
tpms <- NULL
for (s in names(files)) {
  q <- fread(files[s], select = c("Name", "TPM", "NumReads"))
  if (is.null(counts)) {
    counts <- data.table(transcript_id = q$Name)
    tpms <- data.table(transcript_id = q$Name)
  }
  counts[[s]] <- q$NumReads
  tpms[[s]] <- q$TPM
}

fwrite(counts, file.path(outdir, "salmon_counts_transcript.csv"))
fwrite(tpms, file.path(outdir, "salmon_tpm_transcript.csv.gz"))
message(sprintf("transcripts: %d x %d samples", nrow(counts), ncol(counts) - 1))
