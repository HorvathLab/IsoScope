import pandas as pd
import glob
import os

gene_dir = "data/star_output"
files = sorted(glob.glob(f"{gene_dir}/SRR*.ReadsPerGene.out.tab"))
counts = {}
for i, f in enumerate(files):
    srr = os.path.basename(f).replace(".ReadsPerGene.out.tab", "")
    df = pd.read_csv(f, sep="\t", header=None, names=["gene_id", "unstranded", "sense", "antisense"])
    # Skip first 4 summary rows
    df = df.iloc[4:]
    # Strip version numbers
    df["gene_id"] = df["gene_id"].str.split(".").str[0]
    counts[srr] = df.set_index("gene_id")["unstranded"]
    if (i+1) % 100 == 0:
        print(f"  Processed {i+1}/{len(files)}")

count_mat = pd.DataFrame(counts)
count_mat.index.name = "gene_id"
count_mat.to_csv("star_counts_gene.csv")
print(f"Saved: star_counts_gene.csv ({count_mat.shape[0]} genes x {count_mat.shape[1]} samples)")
