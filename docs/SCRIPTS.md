# Pipeline Scripts

Preprocessing scripts are in `scripts/`. They are numbered in the order they run; `<dataset>` is one of `ccle`, `gbm`, `hcc` or `pca`.

| Script | Description |
|---|---|
| `snakemake_processing/` | Snakemake workflow for the patient cohorts: download, STAR genome build and alignment, Salmon quantification |
| `01_tx_counts_<dataset>.R` | Transcript count and TPM matrices from Salmon output |
| `02_tximport_<dataset>.R` | Gene-level counts and TPM with tximport (lengthScaledTPM) |
| `03_transcript_metadata.py` | Transcript annotations and coordinates from the GENCODE v49 basic GTF |
| `04_transcript_exons.py` | Per-exon coordinates from the GENCODE v49 basic GTF |
| `05_pull_protein_annotations.R` | UniProt, Pfam, InterPro from Ensembl BioMart release 115 |
| `06_match_gdsc_response.py` | GDSC1 and GDSC2 drug response matched to CCLE cell lines |
| `07_build_database.R` | Builds one dataset's `quants.db` and summary tables; takes the dataset name as argument |
| `08_deseq2_ccle.R`, `08_deseq2_ccle_transcripts.R` | Pairwise CCLE DESeq2, gene and transcript level |
| `08_deseq2_<dataset>.R` | DESeq2 for each patient cohort with its own design (gbm, hcc, pca) |
| `check_validation.R`, `check_examples.R`, `check_hnf4a_paired.R` | Reproduce the validation values and the README use cases |

## Metadata in `data/`

| File | Description |
|---|---|
| `datasets.csv` | Dataset registry: ID, display name, DESeq2 design |
| `sample_metadata.csv` | CCLE samples, cell lines and tissues |
| `samples_gbm.csv`, `samples_hcc.csv`, `samples_pca.csv` | Patient cohort samples, subjects and groups |
| `tx2gene.csv` | Transcript to gene symbol mapping (GENCODE v49) |

Hugging Face deployment files (Dockerfile, start script) are in `deploy/`.
