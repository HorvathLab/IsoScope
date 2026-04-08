#!/usr/bin/env python
"""
Extract per-exon coordinates from a GENCODE GTF, restricted to transcripts
present in the IsoScope quants.db, and write a compact CSV used by the
Isoform Exon Structure visualization in the Shiny app.

Output schema (transcript_exons.csv.gz):
    gene, transcript_id, exon_number, chr, start, end, strand

transcript_id is version-stripped (ENST00000000123 not ENST00000000123.4)
to match the identifiers used throughout the application.

Usage:
    python 10_extract_transcript_exons.py \\
        --gtf  /path/to/gencode.v49.basic.annotation.gtf.gz \\
        --db   shiny_data/quants.db \\
        --out  shiny_data/transcript_exons.csv.gz
"""

import argparse
import gzip
import re
import sqlite3
import pandas as pd

ATTR_TX = re.compile(r'transcript_id "([^"]+)"')
ATTR_GENE = re.compile(r'gene_name "([^"]+)"')
ATTR_EXON = re.compile(r'exon_number (\d+)')


def load_kept_transcripts(db_path):
    """Return the set of version-stripped transcript_ids present in quants.db."""
    con = sqlite3.connect(db_path)
    tx = pd.read_sql_query("SELECT DISTINCT transcript_id FROM quants", con)
    con.close()
    tx["transcript_id"] = tx["transcript_id"].str.split(".").str[0]
    return set(tx["transcript_id"])


def stream_gtf_exons(gtf_path, kept_tx):
    """Yield per-exon dicts from a GENCODE GTF, filtered to kept transcripts."""
    opener = gzip.open if gtf_path.endswith(".gz") else open
    with opener(gtf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "exon":
                continue
            attrs = fields[8]
            m_tx = ATTR_TX.search(attrs)
            if not m_tx:
                continue
            tx_id = m_tx.group(1).split(".")[0]
            if tx_id not in kept_tx:
                continue
            m_gene = ATTR_GENE.search(attrs)
            m_exon = ATTR_EXON.search(attrs)
            yield {
                "gene": m_gene.group(1) if m_gene else "",
                "transcript_id": tx_id,
                "exon_number": int(m_exon.group(1)) if m_exon else 0,
                "chr": fields[0],
                "start": int(fields[3]),
                "end": int(fields[4]),
                "strand": fields[6],
            }


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gtf", required=True, help="GENCODE GTF (.gtf or .gtf.gz)")
    ap.add_argument("--db", required=True, help="quants.db SQLite file")
    ap.add_argument("--out", required=True, help="Output CSV path (will be gzipped if .gz)")
    args = ap.parse_args()

    print(f"Loading kept transcript ids from {args.db}")
    kept = load_kept_transcripts(args.db)
    print(f"  {len(kept):,} unique transcripts in quants.db")

    print(f"Streaming exons from {args.gtf}")
    rows = list(stream_gtf_exons(args.gtf, kept))
    df = pd.DataFrame(rows, columns=["gene", "transcript_id", "exon_number",
                                     "chr", "start", "end", "strand"])
    df = df.sort_values(["gene", "transcript_id", "exon_number"]).reset_index(drop=True)

    print(f"  {len(df):,} exon rows across {df['transcript_id'].nunique():,} transcripts")
    print(f"Writing {args.out}")
    df.to_csv(args.out, index=False, compression="gzip" if args.out.endswith(".gz") else None)
    print("Done.")


if __name__ == "__main__":
    main()
