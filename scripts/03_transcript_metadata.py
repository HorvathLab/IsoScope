import argparse
import csv
import gzip
import re
from collections import defaultdict

RE_TX = re.compile(r'transcript_id "([^"]+)"')
RE_GENE = re.compile(r'gene_name "([^"]+)"')
RE_GENEID = re.compile(r'gene_id "([^"]+)"')
RE_TYPE = re.compile(r'transcript_type "([^"]+)"')

ap = argparse.ArgumentParser()
ap.add_argument("--gtf", required=True)
ap.add_argument("--outdir", default="out")
args = ap.parse_args()

tx_info = {}
tx_exons = defaultdict(list)
opener = gzip.open if args.gtf.endswith(".gz") else open

with opener(args.gtf, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        fields = line.rstrip("\n").split("\t")
        if len(fields) < 9:
            continue
        attrs = fields[8]
        m_tx = RE_TX.search(attrs)
        if not m_tx:
            continue
        tx_id = m_tx.group(1)
        if fields[2] == "transcript":
            m_gene = RE_GENE.search(attrs)
            m_gid = RE_GENEID.search(attrs)
            m_type = RE_TYPE.search(attrs)
            chrom = fields[0] if fields[0].startswith("chr") else "chr" + fields[0]
            tx_info[tx_id] = {
                "gene": m_gene.group(1) if m_gene else (m_gid.group(1) if m_gid else ""),
                "gene_id": m_gid.group(1) if m_gid else "",
                "biotype": m_type.group(1) if m_type else "unknown",
                "strand": fields[6],
                "chr": chrom,
                "start": int(fields[3]),
                "end": int(fields[4]),
            }
        elif fields[2] == "exon":
            tx_exons[tx_id].append((int(fields[3]), int(fields[4])))

annot = f"{args.outdir.rstrip('/')}/transcript_annotations.csv"
with open(annot, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["transcript_id", "gene", "gene_id", "biotype", "n_exons", "tx_length", "strand", "chr"])
    for tx, info in sorted(tx_info.items()):
        ex = tx_exons.get(tx, [])
        w.writerow([tx, info["gene"], info["gene_id"], info["biotype"],
                    len(ex), sum(e - s + 1 for s, e in ex), info["strand"], info["chr"]])

coords = f"{args.outdir.rstrip('/')}/transcript_coords.csv"
with open(coords, "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["transcript_id", "gene", "chr", "start", "end", "strand"])
    for tx, info in sorted(tx_info.items()):
        w.writerow([tx, info["gene"], info["chr"], info["start"], info["end"], info["strand"]])

print(f"annotations: {len(tx_info)} transcripts -> {annot}")
print(f"coordinates: {len(tx_info)} transcripts -> {coords}")
