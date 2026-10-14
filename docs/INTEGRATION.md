# Adding a Dataset

Each dataset in IsoScope is a folder in `shiny_data/` with its own database. The app scans `shiny_data/` at startup and lists every dataset folder it finds, so adding one needs no changes to `app.R`. The steps below use a hypothetical dataset called `mydata`.

### Step 1: Quantify the samples

Quantify each sample with Salmon against GENCODE v49 basic transcripts, so transcript IDs match the CCLE reference. The Snakemake workflow used for the patient cohorts (STAR alignment, then Salmon in alignment-based mode) is in `scripts/snakemake_processing/`.

### Step 2: Describe the samples

Add `data/samples_mydata.csv` with one row per sample:

| Column | Description | Example |
|---|---|---|
| `sample_id` | Name of the sample's Salmon output folder | `SRR13195827` |
| `sample_label` | Short name shown in plots | `P101_PRE` |
| `subject_id` | Patient or source ID; samples from the same patient share it | `P101` |
| `group_label` | Group used for comparisons | `PRE` |
| `covariate_*` | Optional extra columns, e.g. `covariate_batch` | |

Then add a row to `data/datasets.csv` with the dataset's ID, display name and DESeq2 design (e.g. `~ subject_id + group_label` for paired samples, `~ group_label` otherwise).

### Step 3: Build the dataset

Copy one of the existing patient-cohort scripts (e.g. the `hcc` versions) and change the dataset name and paths:

```bash
Rscript scripts/01_tx_counts_mydata.R      # transcript counts and TPM from Salmon
Rscript scripts/02_tximport_mydata.R       # gene-level counts and TPM (lengthScaledTPM)
Rscript scripts/07_build_database.R mydata # quants.db, group_avg.csv, dominant isoforms, diversity
Rscript scripts/08_deseq2_mydata.R         # optional: precomputed DESeq2 results
```

This produces a folder `mydata/` with the same files as the other datasets (see [Data files and sizes](../README.md#running-locally)).

### Step 4: Launch

Copy the `mydata/` folder into `shiny_data/` and run `shiny::runApp()`. The new dataset appears on the Dataset tab next to the existing cohorts.

In Clustering and Drug Response, the new samples are projected onto the CCLE embedding without changing it. Drug response values exist only for CCLE cell lines, so the new samples show as "No data" there. Without step 3's DESeq2 script, the Differential Expression tab offers only the runtime Wilcoxon test for the new dataset.

> **Note on comparability:** TPM values from different references or library types are not directly comparable. For best results, quantify from raw reads against GENCODE v49 basic, as in step 1.
