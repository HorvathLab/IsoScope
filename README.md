# IsoScope

*Interactive visualization of transcript isoform expression across cancer cell lines*

**Live application:** <https://azkasajjad-isoscope.hf.space> ( we need to change this)

## Overview

Most genes produce multiple isoforms, and which one is expressed can matter as much as how much. A textbook example is BCL2L1, where alternative splicing produces BCL-XL, which protects cells from apoptosis, and BCL-XS, which promotes it (Dou et al., 2021). Same gene, opposite roles in cell survival, depending on a single splicing decision. Cancer cells routinely exploit this kind of switch to evade therapy, sustain growth, and adapt to new tissue environments. Despite this, most cancer transcriptomics resources still report a single expression value per gene, collapsing the very signal that distinguishes one isoform from another.

IsoScope was developed to enable systematic exploration of transcript isoform variation. It is an R/Shiny application for interactive exploration of isoform-level expression across **667 cancer cell lines from the Cancer Cell Line Encyclopedia (CCLE)**, covering **26 tissue types**, with optional overlay of **6 glioblastoma patient samples** (28 tissue groups, 674 samples total) in GBM mode. The application links isoform expression to drug sensitivity from GDSC1 and GDSC2, detects isoform switching between tissues, quantifies per-gene isoform diversity, visualizes exon-level structural differences across isoforms, and supports overlay of external datasets for translational context. All analyses are accessible through a point-and-click interface without writing code.

## Features

- **Gene Viewer**: Stacked bar, grouped bar, and heatmap views of isoform expression across tissues or individual cell lines
- **Isoform Switching**: Detect changes in dominant isoform between tissue groups or cell lines, with bootstrap confidence intervals
- **Isoform Diversity**: Normalized Shannon entropy quantifying how evenly isoforms are expressed per gene per tissue
- **Isoform Exon Structure**: Tile-grid visualization of every isoform's exon composition for any gene, with click-to-zoom into a single-transcript detail view showing real genomic proportions, intron sizes, and per-exon coordinates. Available in Diversity, Tissue-Specific, and Switching tabs.
- **Tissue-Specific Isoforms**: Identify transcripts predominantly expressed in a single tissue, with optional gene-name search to query specific genes
- **Clustering**: PCA and UMAP embeddings at isoform or gene level, with external data projection
- **Tissue Correlation**: Pairwise correlation heatmaps across tissues or cell lines
- **Drug Response**: GDSC1/GDSC2 drug sensitivity overlaid on expression embeddings
- **Differential Expression**: Precomputed DESeq2 (325 tissue pairs) and runtime Wilcoxon tests
- **Pathway Enrichment**: KEGG enrichment via clusterProfiler, directional from DE or pooled from switching
- **Circos**: Genome-wide isoform expression landscape with switching event overlay
- **QC & Distribution**: Per-sample quality metrics across tissue types

## Data

- **CCLE RNA-seq:** 667 cell lines, 26 tissues (BioProject [PRJNA523380](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA523380))
- **Alignment:** STAR 2.7.11b, GRCh38, GENCODE v49 basic annotation
- **Quantification:** Salmon 1.10.3 (alignment-based mode)
- **Drug response:** GDSC1 (378 drugs) and GDSC2 (286 drugs) via [Cell Model Passports](https://cellmodelpassports.sanger.ac.uk/downloads)
- **Annotations:** GENCODE v49 (biotype, exons, length, per-exon coordinates), Ensembl BioMart release 115 (UniProt, Pfam, InterPro)
- **External dataset:** 6 glioblastoma patient samples from [GSE162643](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162643) (pre/post G207 oncolytic HSV treatment)
- **Large Dataset:** Data used in running the app locally that are too large for github can be found [here](https://zenodo.org/records/19475049)

## Running locally

The simplest way to use IsoScope is the [live Hugging Face Space](https://azkasajjad-isoform-explorer.hf.space). To run a local copy:

**Prerequisites:**

- R >= 4.2
- Bioconductor packages: `DESeq2`, `clusterProfiler`, `org.Hs.eg.db`
- CRAN packages: `shiny`, `ggplot2`, `plotly`, `DT`, `dplyr`, `tidyr`, `data.table`, `DBI`, `RSQLite`, `uwot`, `bslib`, `circlize`, `RColorBrewer`

```bash
git clone https://github.com/HorvathLab/IsoScope.git
cd IsoScope
```

Place the following data files in `shiny_data/`: ( is this where zenodo link goes?)

| File | Description |
|------|-------------|
| `quants.db` | Transcript TPM matrix (SQLite, 827 MB) |
| `deseq2.db` | Precomputed DESeq2 results (SQLite, distributed as `.zip`) |
| `sample_metadata.csv` | Sample-to-tissue mapping (674 rows) |
| `tissue_avg.csv` | Tissue-averaged expression |
| `dominant_isoforms.csv` | Dominant isoform per gene per tissue with bootstrap CIs |
| `diversity_range.csv` | Per-gene entropy range across tissues |
| `transcript_annotations.csv` | Biotype, exon count, transcript length |
| `transcript_exons.csv.gz` | Per-exon coordinates from GENCODE v49 (powers the exon structure view) |
| `transcript_coords.csv` | Genomic coordinates for Circos |
| `protein_annotations_collapsed.csv` | UniProt, Pfam, InterPro per transcript |
| `enst_to_refseq.csv` | ENST-to-RefSeq NM mapping (53,946 entries) |
| `gdsc1_matched.csv` | GDSC1 drug response matched to CCLE |
| `gdsc2_matched.csv` | GDSC2 drug response matched to CCLE |
| `tpm_matrix.txt.gz` | Gene-level TPM from STAR (optional, for gene-level clustering) |

**Launch:**

```r
shiny::runApp()
```

## Integrating your own data

The application ships with CCLE data as the default reference. Users can overlay their own samples without modifying any application code.

### Adding samples to the existing CCLE reference

1. Quantify your RNA-seq data with Salmon against any GENCODE reference. If your data was quantified against RefSeq, run Script 09 first to remap NM identifiers to ENST.
2. Run Script 01 to build count/TPM matrices from your `quant.sf` files.
3. Run Script 07 with `--metadata` pointing to a CSV containing both CCLE and your new samples (columns: `SRR_ID`, `CellLine`, `Tissue`). Your samples should have unique tissue labels (e.g., `MY_TUMORS`).
4. The script rebuilds `quants.db`, `tissue_avg.csv`, and all derived tables with your samples merged in.
5. Replace the files in `shiny_data/` and relaunch. Your samples appear alongside CCLE in all tabs.

In embedding-based analyses (Clustering, Drug Response), the application projects external samples onto the CCLE-derived coordinate space, displaying them as distinct markers without distorting the reference landscape. In all other tabs, your samples participate directly in the analysis using their own quantification values.

### Running with your own data only (no CCLE)

Follow the same steps above but with only your samples in the metadata file. The application functions as a general-purpose isoform exploration tool. The Dataset selector and CCLE/GBM toggle will not apply.

> **Note:** TPM values from different quantification contexts (RefSeq vs GENCODE) are approximate when compared directly. Re-quantification from FASTQs against the same reference is recommended when overlaying onto CCLE.

## Pipeline scripts

All preprocessing scripts are provided in `scripts/numbered/`:

| # | Script | Description |
|---|--------|-------------|
| 01 | `build_count_matrix.py` | tx2gene mapping + Salmon count/TPM matrix |
| 02 | `build_gene_counts.py` | STAR gene-level count matrix |
| 03 | `calculate_gene_tpm.R` | Gene counts to TPM |
| 04 | `extract_transcript_metadata.py` | Transcript annotations + genomic coordinates from GTF |
| 05 | `pull_protein_annotations.R` | UniProt, Pfam, InterPro from BioMart |
| 06 | `match_gdsc_response.py` | GDSC drug data matching to CCLE |
| 07 | `preprocess_for_shiny.py` | Build SQLite database + all app-ready tables |
| 08 | `run_deseq2_pairwise.R` | Pairwise DESeq2 + SQLite database |
| 09 | `integrate_external_dataset.R` | RefSeq-to-ENST mapping + external data integration |
| 10 | `extract_transcript_exons.py` | Per-exon coordinates from GENCODE GTF (powers the exon structure view) |

## Citation

Azka Sajjad, Siera Martinez, Luke Johnson, Vania Ballesteros Prieto, Hovhannes Arestakesyan, Jewel Dias, Midrar Alhossiny, and Anelia Horvath. *IsoScope: Interactive visualization of transcript isoform expression across cancer cell lines.* George Washington University, 2026.

## License

Released under the MIT License. See `LICENSE` for full text.


## References

1. Dou Z, Zhao D, Chen X, Xu C, Jin X, Zhang X, Wang Y, Xie X, Li Q, Di C, Zhang H. *Aberrant Bcl-x splicing in cancer: from molecular mechanism to therapeutic modulation.* Journal of Experimental & Clinical Cancer Research, 40, 194, 2021. doi:10.1186/s13046-021-02001-w
2. Ghandi M, et al. *Next-generation characterization of the Cancer Cell Line Encyclopedia.* Nature, 569, 503–508, 2019.
3. Garnett MJ, et al. *Systematic identification of genomic markers of drug sensitivity in cancer cells.* Nature, 483, 570–575, 2012.
4. Markert JM, et al. *Phase Ib trial of mutant herpes simplex virus G207 inoculated pre- and post-tumor resection for recurrent GBM.* Molecular Therapy, 17, 199–207, 2009.
5. Dobin A, Davis CA, Schlesinger F, Drenkow J, Zaleski C, Jha S, Batut P, Chaisson M, Gingeras TR. *STAR: ultrafast universal RNA-seq aligner.* Bioinformatics, 29(1), 15–21, 2013. doi:10.1093/bioinformatics/bts635
6. Patro R, Duggal G, Love MI, Irizarry RA, Kingsford C. *Salmon provides fast and bias-aware quantification of transcript expression.* Nature Methods, 14(4), 417–419, 2017. doi:10.1038/nmeth.4197
7. Love MI, Huber W, Anders S. *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biology, 15, 550, 2014.
8. Martin FJ, et al. *Ensembl 2023.* Nucleic Acids Research, 51, D933–D941, 2023.
9. The UniProt Consortium. *UniProt: the Universal Protein Knowledgebase in 2023.* Nucleic Acids Research, 51, D523–D531, 2023.
10. Mistry J, et al. *Pfam: The protein families database in 2021.* Nucleic Acids Research, 49, D412–D419, 2021.
11. Paysan-Lafosse T, et al. *InterPro in 2022.* Nucleic Acids Research, 51, D418–D427, 2023.

## Contact

Horvath Lab, Department of Biochemistry and Molecular Medicine, The George Washington University. [GitHub](https://github.com/HorvathLab)