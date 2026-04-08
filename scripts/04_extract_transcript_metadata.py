"""
Extract transcript annotations and genomic coordinates from GENCODE GTF.

Input:  GENCODE v49 basic annotation GTF
Output: transcript_annotations.csv (biotype, exons, length)
        transcript_coords.csv (chr, start, end for Circos)
"""
import argparse
import csv
import gzip
import re
from collections import defaultdict


def parse_attributes(attr_str):
    attrs = {}
    for field in attr_str.strip().split(";"):
        field = field.strip()
        if not field:
            continue
        match = re.match(r'(\w+)\s+"([^"]*)"', field)
        if match:
            attrs[match.group(1)] = match.group(2)
    return attrs


def strip_version(eid):
    return eid.split(".")[0]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--outdir", default="shiny_data")
    args = parser.parse_args()

    tx_info = {}
    tx_exons = defaultdict(list)

    opener = gzip.open if args.gtf.endswith(".gz") else open

    with opener(args.gtf, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9:
                continue

            feature_type = fields[2]
            chrom = fields[0]
            strand = fields[6]
            start = int(fields[3])
            end = int(fields[4])
            attrs = parse_attributes(fields[8])

            tx_id_raw = attrs.get("transcript_id", "")
            if not tx_id_raw:
                continue
            tx_id = strip_version(tx_id_raw)

            if feature_type == "transcript":
                gene_name = attrs.get("gene_name", "")
                gene_id = strip_version(attrs.get("gene_id", ""))
                biotype = attrs.get("transcript_type", attrs.get("transcript_biotype", "unknown"))
                gene = gene_name if gene_name else gene_id

                tx_info[tx_id] = {
                    "gene": gene,
                    "gene_id": gene_id,
                    "biotype": biotype,
                    "strand": strand,
                    "chr": chrom if chrom.startswith("chr") else "chr" + chrom,
                    "start": start,
                    "end": end,
                }

            elif feature_type == "exon":
                tx_exons[tx_id].append((start, end))

    outdir = args.outdir.rstrip("/")

    annot_path = f"{outdir}/transcript_annotations.csv"
    with open(annot_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["transcript_id", "gene", "gene_id", "biotype",
                         "n_exons", "tx_length", "strand", "chr"])
        for tx_id, info in sorted(tx_info.items()):
            exons = tx_exons.get(tx_id, [])
            writer.writerow([
                tx_id, info["gene"], info["gene_id"], info["biotype"],
                len(exons), sum(e - s + 1 for s, e in exons),
                info["strand"], info["chr"],
            ])

    coords_path = f"{outdir}/transcript_coords.csv"
    with open(coords_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["transcript_id", "gene", "chr", "start", "end", "strand"])
        for tx_id, info in sorted(tx_info.items()):
            writer.writerow([
                tx_id, info["gene"], info["chr"],
                info["start"], info["end"], info["strand"],
            ])

    print(f"Annotations: {len(tx_info)} transcripts -> {annot_path}")
    print(f"Coordinates: {len(tx_info)} transcripts -> {coords_path}")


if __name__ == "__main__":
    main()
