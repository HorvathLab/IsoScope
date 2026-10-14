library(biomaRt)
library(data.table)
library(dplyr)

outdir <- "shiny_data"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)
httr::set_config(httr::timeout(600))

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl", version = 115)

annot <- read.csv(file.path(outdir, "transcript_annotations.csv"), stringsAsFactors = FALSE)
our_ids <- unique(annot$transcript_id)

query_chunks <- function(attrs, ids, chunk_size = 3000) {
  chunks <- split(ids, ceiling(seq_along(ids) / chunk_size))
  results <- list()
  for (i in seq_along(chunks)) {
    results[[i]] <- tryCatch(
      getBM(attributes = attrs, filters = "ensembl_transcript_id",
            values = chunks[[i]], mart = ensembl),
      error = function(e) {
        sub_chunks <- split(chunks[[i]], ceiling(seq_along(chunks[[i]]) / (chunk_size / 2)))
        do.call(rbind, lapply(sub_chunks, function(sc) {
          Sys.sleep(2)
          tryCatch(getBM(attributes = attrs, filters = "ensembl_transcript_id",
                         values = sc, mart = ensembl), error = function(e2) data.frame())
        }))
      }
    )
    Sys.sleep(1)
  }
  do.call(rbind, results)
}

tx_info <- query_chunks(
  attrs = c("ensembl_transcript_id", "transcript_length", "cds_length", "uniprotswissprot"),
  ids = our_ids
)

pc_ids <- unique(tx_info$ensembl_transcript_id[!is.na(tx_info$cds_length) & tx_info$cds_length > 0])

pfam <- query_chunks(attrs = c("ensembl_transcript_id", "pfam"), ids = pc_ids)
pfam <- pfam[!is.na(pfam$pfam) & pfam$pfam != "", ]

interpro <- query_chunks(
  attrs = c("ensembl_transcript_id", "interpro", "interpro_description"),
  ids = pc_ids
)
interpro <- interpro[!is.na(interpro$interpro) & interpro$interpro != "", ]

merged <- tx_info
if (nrow(pfam) > 0) merged <- merge(merged, pfam, by = "ensembl_transcript_id", all.x = TRUE)
if (nrow(interpro) > 0) merged <- merge(merged, interpro, by = "ensembl_transcript_id", all.x = TRUE)
merged[merged == ""] <- NA
merged <- unique(merged)

collapsed <- merged %>%
  group_by(ensembl_transcript_id) %>%
  summarise(
    transcript_length = first(transcript_length),
    cds_length = first(cds_length),
    uniprot = paste(unique(na.omit(uniprotswissprot)), collapse = ", "),
    pfam_domains = paste(unique(na.omit(pfam)), collapse = ", "),
    interpro_domains = paste(unique(na.omit(interpro_description)), collapse = "; "),
    .groups = "drop"
  ) %>%
  mutate(
    uniprot = ifelse(uniprot == "", NA_character_, uniprot),
    pfam_domains = ifelse(pfam_domains == "", NA_character_, pfam_domains),
    interpro_domains = ifelse(interpro_domains == "", NA_character_, interpro_domains)
  ) %>%
  rename(transcript_id = ensembl_transcript_id)

fwrite(collapsed, file.path(outdir, "protein_annotations_collapsed.csv"))
message(sprintf("Done: %d transcripts written", nrow(collapsed)))
