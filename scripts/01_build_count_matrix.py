"""
Build transcript-to-gene mapping from GENCODE GTF, then construct
a transcript-level count matrix from Salmon quant.sf files.

Input:  GENCODE v49 basic annotation GTF, Salmon quant.sf directories
Output: tx2gene.csv, salmon_counts_transcript.csv
"""
import gzip
import csv
import re
import glob
import pandas as pd
from pathlib import Path

gtf_path   = "data/gencode.v49.basic.annotation.gtf.gz"
salmon_dir = "data/salmon_output"
outdir     = "shiny_data"

RE_TX   = re.compile(r'transcript_id "([^"]+)"')
RE_GENE = re.compile(r'gene_name "([^"]+)"')

# --- Part 1: tx2gene from GTF ---

tx2gene = {}
with gzip.open(gtf_path, "rt") as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        fields = line.split("\t")
        if fields[2] != "transcript":
            continue
        m_tx = RE_TX.search(fields[8])
        m_gene = RE_GENE.search(fields[8])
        if m_tx and m_gene:
            tx2gene[m_tx.group(1)] = m_gene.group(1)

tx2gene_path = f"{outdir}/tx2gene.csv"
with open(tx2gene_path, "w", newline="") as fh:
    writer = csv.writer(fh)
    writer.writerow(["transcript_id", "gene_name"])
    writer.writerows(tx2gene.items())

print(f"tx2gene: {len(tx2gene)} transcripts")

# --- Part 2: Salmon count matrix ---

quant_files = sorted(glob.glob(f"{salmon_dir}/SRR*_salmon/quant.sf/quant.sf"))
if not quant_files:
    quant_files = sorted(glob.glob(f"{salmon_dir}/SRR*/quant.sf"))

counts = {}
for i, qf in enumerate(quant_files):
    parts = Path(qf).parts
    srr = [p for p in parts if p.startswith("SRR")][0].replace("_salmon", "")
    df = pd.read_csv(qf, sep="\t")
    df["Name"] = df["Name"].str.split(".").str[0]
    counts[srr] = df.set_index("Name")["NumReads"]
    if (i + 1) % 100 == 0:
        print(f"  Processed {i + 1}/{len(quant_files)}")

count_mat = pd.DataFrame(counts)
count_mat.index.name = "transcript_id"
count_mat.to_csv(f"{outdir}/salmon_counts_transcript.csv")
print(f"Counts: {count_mat.shape[0]} transcripts x {count_mat.shape[1]} samples")

# --- Part 3: TPM matrix for preprocessing ---

tpms = {}
for i, qf in enumerate(quant_files):
    parts = Path(qf).parts
    srr = [p for p in parts if p.startswith("SRR")][0].replace("_salmon", "")
    df = pd.read_csv(qf, sep="\t")
    df["Name"] = df["Name"].str.split(".").str[0]
    tpms[srr] = df.set_index("Name")["TPM"]

tpm_mat = pd.DataFrame(tpms)
tpm_mat.index.name = "ENST"
tpm_mat.to_csv(f"{outdir}/all_quants_all_ENST.csv")
print(f"TPM: {tpm_mat.shape[0]} transcripts x {tpm_mat.shape[1]} samples")
