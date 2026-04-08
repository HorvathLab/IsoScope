"""
Converts Salmon TPM output and sample metadata into the files required by
the Isoform Explorer Shiny app. Transcript IDs are mapped to gene names via
tx2gene.csv, zero-expression transcripts are removed, and all computationally
expensive operations (tissue averages, dominant isoforms with bootstrap CIs,
Shannon diversity, QC stats, isoform variance) are precomputed here so the
app never performs matrix-level operations at runtime.

Usage:
    python 07_preprocess_for_shiny.py \
        --quants   data/all_quants_all_ENST.csv \
        --metadata data/sample_metadata.csv \
        --tx2gene  shiny_data/tx2gene.csv \
        --outdir   shiny_data
"""
import os
import argparse
import sqlite3
import numpy as np
import pandas as pd

N_BOOT = 1000
SEED   = 42


def bootstrap_dominant_prop(gene_tpm, col_idx, n_boot, rng):
    n_samples   = len(col_idx)
    tissue_vals = gene_tpm[:, col_idx]
    observed_avg = tissue_vals.mean(axis=1)
    total = observed_avg.sum()
    if total == 0:
        return None

    observed_props = observed_avg / total
    top_idx        = int(np.argmax(observed_props))
    observed_prop  = float(observed_props[top_idx])

    if n_samples < 3:
        return top_idx, observed_prop, np.nan, np.nan, n_samples

    resample_idx = rng.integers(0, n_samples, size=(n_boot, n_samples))
    boot_avgs    = tissue_vals[:, resample_idx].mean(axis=2)
    boot_totals  = boot_avgs.sum(axis=0)
    safe_totals  = np.where(boot_totals == 0, 1.0, boot_totals)
    boot_props   = boot_avgs[top_idx, :] / safe_totals
    boot_props   = np.where(boot_totals == 0, 0.0, boot_props)

    ci_lower = float(np.percentile(boot_props, 2.5))
    ci_upper = float(np.percentile(boot_props, 97.5))
    return top_idx, observed_prop, ci_lower, ci_upper, n_samples


def main(quants_path, metadata_path, outdir, tx2gene_path):
    os.makedirs(outdir, exist_ok=True)

    print("Loading metadata...")
    meta = pd.read_csv(metadata_path)
    srr_to_cell   = dict(zip(meta['SRR_ID'], meta['CellLine']))
    srr_to_tissue = dict(zip(meta['SRR_ID'], meta['Tissue']))

    print("Loading quant matrix...")
    df = pd.read_csv(quants_path, index_col=0)

    tx2gene_dict = dict(zip(*pd.read_csv(tx2gene_path).values.T))
    enst_col = df['ENST']
    gene_col = enst_col.map(tx2gene_dict)
    mask     = gene_col.notna()
    df       = df[mask].reset_index(drop=True)
    gene_col = gene_col[mask].reset_index(drop=True)
    enst_col = enst_col[mask].reset_index(drop=True)

    sample_cols = [c for c in df.columns if c != 'ENST']
    df_samples  = df[sample_cols].rename(columns={c: c.replace('_salmon', '') for c in sample_cols})
    valid_srrs  = [s for s in df_samples.columns if s in srr_to_cell]
    df_samples  = df_samples[valid_srrs]

    cell_names = [srr_to_cell[s] for s in valid_srrs]
    if len(set(cell_names)) < len(cell_names):
        srr_to_orig = dict(zip(meta['SRR_ID'], meta['OriginalSampleName']))
        col_labels  = [srr_to_orig[s] for s in valid_srrs]
    else:
        col_labels = cell_names
    df_samples.columns = col_labels

    nonzero_mask = df_samples.sum(axis=1) > 0
    df_out = pd.concat([
        pd.DataFrame({'gene': gene_col[nonzero_mask].values,
                       'transcript_id': enst_col[nonzero_mask].values}),
        df_samples[nonzero_mask].reset_index(drop=True)
    ], axis=1)

    mat_np   = df_out.iloc[:, 2:].values
    genes_np = df_out['gene'].values

    # --- Sample metadata ---
    sample_meta_out = pd.DataFrame({
        'sample_label': col_labels,
        'SRR_ID':       valid_srrs,
        'CellLine':     [srr_to_cell[s] for s in valid_srrs],
        'Tissue':       [srr_to_tissue[s] for s in valid_srrs]
    })
    sample_meta_out.to_csv(os.path.join(outdir, 'sample_metadata.csv'), index=False)

    # --- SQLite database ---
    db_out = os.path.join(outdir, 'quants.db')
    if os.path.exists(db_out):
        os.remove(db_out)
    conn = sqlite3.connect(db_out)

    print("Writing quants table...")
    df_out.to_sql('quants', conn, index=False, if_exists='replace', chunksize=200)
    conn.execute('CREATE INDEX idx_gene ON quants (gene)')
    conn.execute('CREATE INDEX idx_transcript ON quants (transcript_id)')
    conn.commit()

    print("Computing QC stats...")
    qc_rows = []
    for i, s in enumerate(col_labels):
        col = mat_np[:, i]
        pos_mask = col > 0
        n_iso  = int(pos_mask.sum())
        n_genes = int(len(np.unique(genes_np[pos_mask]))) if n_iso > 0 else 0
        qc_rows.append({
            'sample_label': s,
            'total_tpm':    float(col.sum()),
            'n_isoforms':   n_iso,
            'n_genes':      n_genes,
            'iso_per_gene': n_iso / max(n_genes, 1)
        })
    pd.DataFrame(qc_rows).to_sql('qc_stats', conn, index=False, if_exists='replace')

    print("Computing isoform variance...")
    pd.DataFrame({
        'gene':          genes_np,
        'transcript_id': df_out['transcript_id'].values,
        'variance':      mat_np.var(axis=1).astype(float)
    }).to_sql('isoform_variance', conn, index=False, if_exists='replace')
    conn.execute('CREATE INDEX idx_var ON isoform_variance (variance DESC)')
    conn.execute('CREATE INDEX idx_gene_var ON isoform_variance (gene)')
    conn.commit()
    conn.close()

    # --- Tissue averages ---
    print("Computing tissue averages...")
    tissue_labels = list(sample_meta_out['Tissue'].unique())
    sample_column_list = list(df_out.columns[2:])

    tissue_sample_map = {
        tissue: [sample_column_list.index(c)
                 for c in sample_meta_out[sample_meta_out['Tissue'] == tissue]['sample_label'].tolist()
                 if c in sample_column_list]
        for tissue in tissue_labels
    }

    tissue_avg = {'gene': genes_np, 'transcript_id': df_out['transcript_id'].values}
    for tissue in tissue_labels:
        t_idx = tissue_sample_map[tissue]
        tissue_avg[tissue] = mat_np[:, t_idx].mean(axis=1) if t_idx else np.zeros(len(df_out))
    pd.DataFrame(tissue_avg).to_csv(os.path.join(outdir, 'tissue_avg.csv'), index=False)

    # --- Dominant isoforms ---
    print(f"Computing dominant isoforms (bootstrap n={N_BOOT})...")
    rng = np.random.default_rng(SEED)
    genes_ta = pd.DataFrame(tissue_avg)['gene'].values
    txids_ta = pd.DataFrame(tissue_avg)['transcript_id'].values
    unique_genes_ta = np.unique(genes_ta)

    dominance_rows = []
    for i, gene in enumerate(unique_genes_ta):
        if i % 5000 == 0 and i > 0:
            print(f"  {i:,} / {len(unique_genes_ta):,} genes")
        gene_mask = genes_ta == gene
        gene_tpm  = mat_np[gene_mask, :]
        txids     = txids_ta[gene_mask]

        for tissue, col_idx in tissue_sample_map.items():
            if not col_idx:
                continue
            result = bootstrap_dominant_prop(gene_tpm, col_idx, N_BOOT, rng)
            if result is None:
                continue
            top_idx, prop, ci_lo, ci_hi, n_samp = result
            dominance_rows.append({
                'gene': gene, 'tissue': tissue,
                'dominant_tx': txids[top_idx], 'dominant_prop': prop,
                'ci_lower': ci_lo, 'ci_upper': ci_hi,
                'n_samples': n_samp,
                'total_tpm': float(gene_tpm[:, col_idx].mean(axis=1).sum())
            })
    pd.DataFrame(dominance_rows).to_csv(os.path.join(outdir, 'dominant_isoforms.csv'), index=False)

    # --- Shannon diversity ---
    print("Computing Shannon diversity...")
    unique_genes = np.unique(genes_np)
    diversity_rows = []
    for i, gene in enumerate(unique_genes):
        if i % 10000 == 0 and i > 0:
            print(f"  {i:,} / {len(unique_genes):,} genes")
        gene_mat = mat_np[genes_np == gene, :]
        entropies = {}
        for tissue, col_idx in tissue_sample_map.items():
            if not col_idx:
                entropies[tissue] = 0.0
                continue
            avg_tpm = gene_mat[:, col_idx].mean(axis=1)
            avg_tpm = avg_tpm[avg_tpm > 0]
            if len(avg_tpm) <= 1:
                entropies[tissue] = 0.0
            else:
                p = avg_tpm / avg_tpm.sum()
                entropies[tissue] = float(-np.sum(p * np.log2(p))) / np.log2(len(avg_tpm))
        vals = list(entropies.values())
        keys = list(entropies.keys())
        diversity_rows.append({
            'gene': gene,
            'entropy_range': max(vals) - min(vals),
            'most_diverse_tissue': keys[vals.index(max(vals))],
            'max_entropy': max(vals),
            'least_diverse_tissue': keys[vals.index(min(vals))],
            'min_entropy': min(vals)
        })
    pd.DataFrame(diversity_rows).sort_values(
        'entropy_range', ascending=False
    ).to_csv(os.path.join(outdir, 'diversity_range.csv'), index=False)

    print(f"\nDone: {df_out['gene'].nunique():,} genes, {len(df_out):,} transcripts")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--quants',   required=True)
    parser.add_argument('--metadata', required=True)
    parser.add_argument('--tx2gene',  required=True)
    parser.add_argument('--outdir',   required=True)
    parser.add_argument('--n_boot',   type=int, default=1000)
    args = parser.parse_args()
    N_BOOT = args.n_boot
    main(os.path.expanduser(args.quants),
         os.path.expanduser(args.metadata),
         os.path.expanduser(args.outdir),
         os.path.expanduser(args.tx2gene))
