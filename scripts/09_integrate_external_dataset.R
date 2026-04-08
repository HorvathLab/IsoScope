library(biomaRt)
library(data.table)

tx2gene_path <- "shiny_data/tx2gene.csv"
input_dir    <- "data/external_quant_sf"
output_dir   <- "data/remapped_quant_sf"
mapping_out  <- "shiny_data/enst_to_refseq.csv"

# --- Part 1: Build ENST-to-RefSeq mapping from BioMart ---

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

raw <- getBM(
  attributes = c("ensembl_transcript_id", "refseq_mrna", "transcript_mane_select"),
  mart = ensembl
)

mane <- raw[raw$transcript_mane_select != "", ]
mane$refseq <- gsub("\\..*", "", mane$transcript_mane_select)
mane <- mane[!duplicated(mane$ensembl_transcript_id), ]

non_mane <- raw[raw$transcript_mane_select == "" & raw$refseq_mrna != "", ]
non_mane$refseq <- non_mane$refseq_mrna
non_mane <- non_mane[!duplicated(non_mane$ensembl_transcript_id), ]
non_mane <- non_mane[!non_mane$ensembl_transcript_id %in% mane$ensembl_transcript_id, ]

mapping <- rbind(
  mane[, c("ensembl_transcript_id", "refseq")],
  non_mane[, c("ensembl_transcript_id", "refseq")]
)
mapping <- mapping[mapping$refseq != "", ]

write.csv(mapping, mapping_out, row.names = FALSE)
message(sprintf("Mapping: %d total (%d MANE Select, %d fallback)",
                nrow(mapping), nrow(mane), nrow(mapping) - nrow(mane)))

# --- Part 2: Remap external RefSeq quant.sf files to ENST ---

if (!dir.exists(input_dir)) {
  message("No external quant.sf directory found. Skipping remap step.")
  quit(save = "no")
}

if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

nm_to_enst <- setNames(mapping$ensembl_transcript_id, mapping$refseq)

tx2gene <- fread(tx2gene_path)
enst_versioned <- setNames(tx2gene$transcript_id, sub("\\..*", "", tx2gene$transcript_id))

quant_files <- list.files(input_dir, pattern = "quant\\.sf$", recursive = TRUE, full.names = TRUE)

for (qf in quant_files) {
  sample_name <- basename(dirname(qf))
  df <- fread(qf)

  df$nm_clean <- gsub("\\..*", "", df$Name)
  df$enst_base <- nm_to_enst[df$nm_clean]
  df <- df[!is.na(df$enst_base), ]

  df$enst_versioned <- enst_versioned[df$enst_base]
  df <- df[!is.na(df$enst_versioned), ]

  tpm_total <- sum(df$TPM)
  df$TPM <- df$TPM / tpm_total * 1e6

  out <- data.frame(
    Name = df$enst_versioned,
    Length = df$Length,
    EffectiveLength = df$EffectiveLength,
    TPM = df$TPM,
    NumReads = df$NumReads
  )

  outpath <- file.path(output_dir, sample_name, "quant.sf")
  dir.create(dirname(outpath), showWarnings = FALSE, recursive = TRUE)
  fwrite(out, outpath, sep = "\t")
  message(sprintf("  %s: %d mapped, TPM renormalized", sample_name, nrow(out)))
}

message(sprintf("Done: %d samples remapped", length(quant_files)))
