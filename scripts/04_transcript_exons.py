import argparse
import gzip
import re
import sqlite3
import pandas as pd

RE_TX = re.compile(r'transcript_id "([^"]+)"')
RE_GENE = re.compile(r'gene_name "([^"]+)"')
RE_EXON = re.compile(r'exon_number (\d+)')

ap = argparse.ArgumentParser()
ap.add_argument("--gtf", required=True)
ap.add_argument("--db", required=True)
ap.add_argument("--out", required=True)
args = ap.parse_args()

con = sqlite3.connect(args.db)
kept = set(pd.read_sql_query("SELECT DISTINCT transcript_id FROM quants", con)["transcript_id"])
con.close()
print(f"{len(kept):,} transcripts in {args.db}")

rows = []
opener = gzip.open if args.gtf.endswith(".gz") else open
with opener(args.gtf, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9 or fields[2] != "exon":
            continue
        attrs = fields[8]
        m_tx = RE_TX.search(attrs)
        if not m_tx or m_tx.group(1) not in kept:
            continue
        m_gene = RE_GENE.search(attrs)
        m_exon = RE_EXON.search(attrs)
        rows.append({
            "gene": m_gene.group(1) if m_gene else "",
            "transcript_id": m_tx.group(1),
            "exon_number": int(m_exon.group(1)) if m_exon else 0,
            "chr": fields[0],
            "start": int(fields[3]),
            "end": int(fields[4]),
            "strand": fields[6],
        })

df = pd.DataFrame(rows, columns=["gene", "transcript_id", "exon_number", "chr", "start", "end", "strand"])
df = df.sort_values(["gene", "transcript_id", "exon_number"]).reset_index(drop=True)
print(f"{len(df):,} exon rows across {df['transcript_id'].nunique():,} transcripts")
df.to_csv(args.out, index=False, compression="gzip" if args.out.endswith(".gz") else None)
print(f"wrote {args.out}")
