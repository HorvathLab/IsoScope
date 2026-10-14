# Adding a Dataset

Each dataset in IsoScope is a folder in `shiny_data/` with its own database. The app scans `shiny_data/` at startup and lists every dataset folder it finds, so adding one needs no changes to `app.R`. The steps below use a hypothetical dataset called `mydata` and are run from the repository root.

### Step 1: Quantify the samples

Align reads to the GRCh38 genome and quantify each sample with Salmon against GENCODE v49 basic transcripts ([release 49 files](https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/); use the GRCh38 files in the main folder, not the `GRCh37_mapping` subfolder). The app matches transcripts by their full versioned ID (e.g. ENST00000316099.10), so data quantified against other GENCODE releases will only partly overlap with the CCLE reference. The reference release is fixed for reproducibility, as in other large resources such as the GDC, which uses GENCODE v36. The Snakemake workflow used for the patient cohorts (STAR alignment, then Salmon in alignment-based mode) is in `scripts/snakemake_processing/`.

### Step 2: Describe the samples

Add `data/samples_mydata.csv` with one row per sample:

| Column | Description | Example |
|---|---|---|
| `sample_id` | Sample name, matching its Salmon output folder | `SRR13195827` |
| `sample_label` | Short name shown in plots | `P101_PRE` |
| `subject_id` | Patient or source ID; samples from the same patient share it | `P101` |
| `group_label` | Group used for comparisons | `PRE` |
| `covariate_*` | Optional extra columns, e.g. `covariate_batch` | |

Then add a row for the dataset to `data/datasets.csv`, following the existing rows:

```
mydata,My dataset,0,data/samples_mydata.csv,out/mydata/salmon_tpm_transcript.csv.gz,~ subject_id + group_label,salmon,gencode,v49,<accession>,,
```

The columns are `dataset_id`, `display_name`, `is_reference` (0 for anything other than CCLE), `samples_file`, `tpm_file`, `design_formula` (`~ subject_id + group_label` for paired samples, `~ group_label` otherwise), `quantifier`, `annotation`, `annotation_version`, `source_accession`, `remapped_from` and `mapping_rate` (the last two can be left empty).

### Step 3: Build the dataset

Copy the `hcc` scripts as templates:

```bash
cp scripts/01_tx_counts_hcc.R scripts/01_tx_counts_mydata.R
cp scripts/02_tximport_hcc.R  scripts/02_tximport_mydata.R
cp scripts/08_deseq2_hcc.R    scripts/08_deseq2_mydata.R
```

Edit the copies before running them:

- **01 and 02:** set `salmon_root` to the folder holding the Salmon output and change `out/hcc` to `out/mydata`. The scripts look for one folder per sample named `<sample_id>_salmon` containing `quant.sf/quant.sf`, and the pattern `"^SRR.*_salmon$"` expects sample IDs starting with `SRR`; adjust the pattern and path to match your layout.
- **08:** set the samples file, output paths and `deseq2_mydata.db`, the group names in `levels` and `contrast`, and the design. Set `EXCLUDE` to a sample to leave out, or to `""` to keep all samples. The pairing step keeps only subjects with a sample in both groups; remove it for an unpaired design.
- The first line of each script adds a library path on our server; it is ignored elsewhere and can be deleted.

Then run:

```bash
Rscript scripts/01_tx_counts_mydata.R      # transcript counts and TPM from Salmon
Rscript scripts/02_tximport_mydata.R       # gene-level counts and TPM (lengthScaledTPM)
Rscript scripts/07_build_database.R mydata # quants.db, group averages, dominant isoforms, diversity
Rscript scripts/08_deseq2_mydata.R         # optional: precomputed DESeq2 results
```

Scripts 01, 07 and 08 need `data.table`, `DBI`, `RSQLite` and `matrixStats`; 02 needs `tximport` and 08 needs `DESeq2`. The results are written to `out/mydata/`: `quants.db`, `group_avg.csv`, `dominant_isoforms.csv`, `diversity_range.csv` and, from step 08, `deseq2_mydata.db`, the same files as the other datasets (see [Data files and sizes](../README.md#running-locally)).

### Step 4: Launch

Copy `out/mydata/` to `shiny_data/mydata/` and run `shiny::runApp()`. The new dataset appears on the Dataset tab next to the existing cohorts.

In Clustering and Drug Response, the new samples are projected onto the CCLE embedding without changing it. Drug response values exist only for CCLE cell lines, so the new samples show as "No data" there. Without the DESeq2 step, the Differential Expression tab offers only the runtime Wilcoxon test for the new dataset.

> **Note on comparability:** TPM values from different references or library types are not directly comparable. For best results, quantify from raw reads against GENCODE v49 basic, as in step 1.
