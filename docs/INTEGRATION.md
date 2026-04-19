# Integrating Your Own Data

IsoScope can incorporate external RNA-seq datasets alongside the CCLE reference. This requires preprocessing your data, merging it into the existing data files, and editing a few lines in `app.R` so the application recognizes your samples. The walkthrough below assumes you have already set up `shiny_data/` as described in the [Running Locally](../README.md#running-locally) section of the main README.

### Step 1: Quantify your samples

Produce Salmon `quant.sf` files for your samples keyed by ENST transcript IDs from a GENCODE release. CCLE was quantified against GENCODE v49; using the same release gives the cleanest join, but any recent GENCODE release will overlap substantially and work in practice.

If your data was already quantified against a different reference (e.g., RefSeq with NM accessions), use Script 09 to remap identifiers to ENST. Edit the paths at the top of the script before running:

```r
# In scripts/09_integrate_external_dataset.R, update these:
input_dir    <- "path/to/your/quant_sf_files"
output_dir   <- "path/to/remapped_output"
mapping_out  <- "shiny_data/enst_to_refseq.csv"
```

### Step 2: Build the TPM matrix

Edit the `salmon_dir` path at the top of `scripts/01_build_count_matrix.py` to point to the directory containing your Salmon output (each subdirectory should contain a `quant.sf` file):

```python
# In scripts/01_build_count_matrix.py, update:
salmon_dir = "path/to/your/salmon_output"
outdir     = "shiny_data"
```

Running the script produces three files:
- `tx2gene.csv`: transcript-to-gene mapping from the GENCODE GTF
- `salmon_counts_transcript.csv`: integer read counts (used by DESeq2)
- `all_quants_all_ENST.csv`: transcript TPM matrix (one column per sample)

The last file is the input for script 07. At this stage it contains only your samples.

### Step 3: Merge with the CCLE TPM matrix

Script 07 expects a TPM matrix that includes both the CCLE reference and your new samples. Download `all_quants_all_ENST.csv` from Zenodo (https://zenodo.org/records/19475049) and append your columns:

```python
import pandas as pd

# Load the CCLE TPM matrix from Zenodo
ccle = pd.read_csv("shiny_data/all_quants_all_ENST.csv", index_col=0).set_index("ENST")

# Load your samples (produced by Step 2)
mine = pd.read_csv("shiny_data/all_quants_mine.csv", index_col=0).set_index("ENST")

# Outer-join on ENST, fill missing with zero, write combined matrix.
# reset_index() preserves ENST as a column (not index), which is what script 07 reads.
combined = ccle.join(mine, how="outer").fillna(0)
combined.index.name = "ENST"
combined.reset_index().to_csv("shiny_data/all_quants_combined.csv")
print(f"Combined: {combined.shape[0]} transcripts x {combined.shape[1]} samples")
```

Save your output from Step 2 as `all_quants_mine.csv` (or rename the file path in the snippet above). Pass `all_quants_combined.csv` as `--quants` in Step 5.

### Step 4: Prepare metadata

Create or update `sample_metadata.csv` to include both CCLE samples and your new samples. The file has three required columns:

| Column | Description | Example |
|--------|-------------|---------|
| `SRR_ID` | Unique sample identifier (does not have to be an SRR accession) | `PATIENT1_PRE` |
| `CellLine` | Short sample name used in plots | `Patient1` |
| `Tissue` | Group label for your samples (must be consistent across replicates) | `MY_PRETREATED` |

For example, if you have 4 tumor samples from 2 patients (pre/post treatment), your rows might look like:

```
SRR_ID,CellLine,Tissue
PATIENT1_PRE,Patient1,PRE_TREATMENT
PATIENT2_PRE,Patient2,PRE_TREATMENT
PATIENT1_POST,Patient1,POST_TREATMENT
PATIENT2_POST,Patient2,POST_TREATMENT
```

Append these rows to the existing `sample_metadata.csv` that contains the 673 CCLE + GBM entries.

### Step 5: Rebuild app data files

Run Script 07 to regenerate the SQLite database and all precomputed tables:

```bash
python scripts/07_preprocess_for_shiny.py \
    --quants   shiny_data/all_quants_combined.csv \
    --metadata shiny_data/sample_metadata.csv \
    --tx2gene  shiny_data/tx2gene.csv \
    --outdir   shiny_data
```

This rebuilds `quants.db`, `tissue_avg.csv`, `dominant_isoforms.csv`, and the diversity and variance tables with your samples included.

If you want precomputed DESeq2 results for comparisons involving your new tissue groups, run Script 08 as well. The script reads pairs from `sample_metadata.csv` and writes to `deseq2.db`:

```bash
Rscript scripts/08_run_deseq2_pairwise.R \
    --counts   shiny_data/salmon_counts_transcript.csv \
    --metadata shiny_data/sample_metadata.csv \
    --outdir   shiny_data
```

### Step 6: Register your tissue groups in app.R

The application needs to know which tissue labels belong to your external dataset. Five sections of `app.R` require editing. The line numbers below reflect the current release and may drift slightly as the code evolves.

**1. Tissue group vectors (around line 300):**

The current code looks like this:

```r
CCLE_TISSUES <- tissue_list[!tissue_list %in% c("PRE_G207", "POST_G207")]
GBM_TISSUES  <- c("PRE_G207", "POST_G207")
```

Add your tissue labels, and exclude them from `CCLE_TISSUES` as well:

```r
CCLE_TISSUES <- tissue_list[!tissue_list %in% c("PRE_G207", "POST_G207", "PRE_TREATMENT", "POST_TREATMENT")]
GBM_TISSUES  <- c("PRE_G207", "POST_G207")
MY_TISSUES   <- c("PRE_TREATMENT", "POST_TREATMENT")
```

**2. Dataset selector UI (around line 1844):**

The current selector offers two modes:

```r
radioButtons(
  "dataset_mode",
  NULL,
  choices = c(
    "CCLE Atlas (667 cell lines, 26 tissues)" = "ccle",
    "CCLE + GBM Patient Samples (673 samples, 28 tissues)" = "both"
  ),
  selected = "ccle"
)
```

Add a new option for your dataset (replace `N` and `M` with your own sample and tissue counts):

```r
radioButtons(
  "dataset_mode",
  NULL,
  choices = c(
    "CCLE Atlas (667 cell lines, 26 tissues)" = "ccle",
    "CCLE + GBM Patient Samples (673 samples, 28 tissues)" = "both",
    "CCLE + My Dataset (N samples, M tissues)" = "my_data"
  ),
  selected = "ccle"
)
```

**3. Reactive expressions for active tissues and cell lines (around lines 1972 and 1983):**

Two reactives filter the active samples and tissues based on `input$dataset_mode`. Add a `my_data` branch to each `switch()`:

```r
active_tissues <- reactive({
  ds <- if (is.null(input$dataset_mode))
    "ccle"
  else
    input$dataset_mode
  switch(ds,
         ccle = CCLE_TISSUES,
         gbm = GBM_TISSUES,
         both = tissue_list,
         my_data = c(CCLE_TISSUES, MY_TISSUES))
})

active_cell_lines <- reactive({
  ds <- if (is.null(input$dataset_mode))
    "ccle"
  else
    input$dataset_mode
  switch(
    ds,
    ccle = sort(unique(sample_meta_df$CellLine[sample_meta_df$Tissue %in% CCLE_TISSUES])),
    gbm = sort(unique(sample_meta_df$CellLine[sample_meta_df$Tissue %in% GBM_TISSUES])),
    both = cell_line_list,
    my_data = sort(unique(sample_meta_df$CellLine[sample_meta_df$Tissue %in% c(CCLE_TISSUES, MY_TISSUES)]))
  )
})
```

**4. Clustering tab projection (around lines 4656 and 4713):**

The Clustering tab projects external samples onto the CCLE-derived embedding. Two conditionals currently check only for `"both"`. Update both to also handle `"my_data"`.

Gene-level warning (line 4656):

```r
# Before
if (!is.null(input$dataset_mode) && input$dataset_mode == "both") {
  showNotification(
    "Gene-level clustering shows CCLE only. GBM samples are available in Isoform (Salmon TPM) mode.",
    type = "warning",
    duration = 6
  )
}

# After
if (!is.null(input$dataset_mode) && input$dataset_mode %in% c("both", "my_data")) {
  showNotification(
    "Gene-level clustering shows CCLE only. External samples are available in Isoform (Salmon TPM) mode.",
    type = "warning",
    duration = 6
  )
}
```

Isoform-level projection (line 4713). The existing block pulls GBM samples and projects them onto the CCLE embedding. Generalize it to pick the right external tissue set based on the mode:

```r
# Before
if (!is.null(input$dataset_mode) &&
    input$dataset_mode == "both") {
  gbm_slabels <- sample_meta_df$sample_label[sample_meta_df$Tissue %in% GBM_TISSUES]
  gbm_cols <- intersect(intersect(sample_cols, gbm_slabels), names(mat_df))
  ...
}

# After
if (!is.null(input$dataset_mode) &&
    input$dataset_mode %in% c("both", "my_data")) {
  ext_tissues <- if (input$dataset_mode == "both") GBM_TISSUES else MY_TISSUES
  gbm_slabels <- sample_meta_df$sample_label[sample_meta_df$Tissue %in% ext_tissues]
  gbm_cols <- intersect(intersect(sample_cols, gbm_slabels), names(mat_df))
  ...
}
```

The internal variable names (`gbm_slabels`, `gbm_cols`, `gbm_tmat`, `gbm_info`) are retained to minimize diff size. They now refer to whichever external cohort is active.

**5. Drug Response tab projection (around line 6316):**

The same pattern applies to the Drug Response tab:

```r
# Before
if (!is.null(input$dataset_mode) &&
    input$dataset_mode == "both") {
  gbm_slabels_drug <- sample_meta_df$sample_label[sample_meta_df$Tissue %in% GBM_TISSUES]
  gbm_cols_drug <- intersect(gbm_slabels_drug, rownames(tmat))
  ...
}

# After
if (!is.null(input$dataset_mode) &&
    input$dataset_mode %in% c("both", "my_data")) {
  ext_tissues <- if (input$dataset_mode == "both") GBM_TISSUES else MY_TISSUES
  gbm_slabels_drug <- sample_meta_df$sample_label[sample_meta_df$Tissue %in% ext_tissues]
  gbm_cols_drug <- intersect(gbm_slabels_drug, rownames(tmat))
  ...
}
```

As with GBM samples, drug response annotations are not available for external samples and their markers will be labeled "No data" in the Drug Response tab.

### Step 7: Launch

With the regenerated files in `shiny_data/` and `app.R` updated, run:

```r
shiny::runApp()
```

In embedding-based analyses (Clustering, Drug Response), your samples are projected onto the CCLE-derived coordinate space and displayed as distinct markers without distorting the reference landscape. In all other tabs, your samples participate directly in the analysis using their own quantification values. Drug response annotations are not available for external samples; they will appear as "No data" in the Drug Response tab.

> **Note on comparability:** TPM values computed against different references (e.g., RefSeq vs GENCODE) are not directly comparable because transcript catalogs and effective lengths differ. For best results, re-quantify from FASTQs against GENCODE v49.
