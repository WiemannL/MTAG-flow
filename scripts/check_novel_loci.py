import pandas as pd
import sys
import os

# Usage: python check_novel_loci.py <top_signals.txt> <gwas_file> <snp_col> <pval_col> <effect_col> <output_dir>

def main():
    if len(sys.argv) < 7:
        print("Usage: python check_novel_loci.py <top_signals.txt> <gwas_file> <snp_col> <pval_col> <effect_col> <output_dir>")
        sys.exit(1)
    top_signals_file = sys.argv[1]
    gwas_file = sys.argv[2]
    snp_col = sys.argv[3]
    pval_col = sys.argv[4]
    effect_col = sys.argv[5]
    output_dir = sys.argv[6]
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)


    # Load top signals
    top_snps_df = pd.read_csv(top_signals_file)
    top_snps = top_snps_df['Top_SNP'].astype(str).tolist()
    # Load GWAS summary stats
    gwas = pd.read_csv(gwas_file, sep=None, engine='python')
    gwas[snp_col] = gwas[snp_col].astype(str)

    # Filter for top SNPs present in GWAS
    found = gwas[gwas[snp_col].isin(top_snps)].copy()
    found['is_significant'] = found[pval_col] < 0.01

    # Save results
    found.to_csv(os.path.join(output_dir, 'top_signals_in_OGgwas.csv'), index=False)
    not_found = set(top_snps) - set(found[snp_col])
    with open(os.path.join(output_dir, 'top_signals_not_in_OGgwas.txt'), 'w') as f:
        for snp in not_found:
            f.write(f"{snp}\n")

    # --- Novelty check: not within ±500kb of any OG GWAS locus ---
    # Use all genome-wide significant loci from OG GWAS as known loci
    # (P < 5e-8)
    gws_loci = gwas[gwas[pval_col] < 5e-8][['CHR', 'BP']].drop_duplicates()
    # Merge top signals with their positions
    top_snps_pos = found[['SNP', 'CHR', 'BP', 'is_significant']].copy()
    # For each top signal, check if it is within ±500kb of any known locus
    novel_snps = []
    for idx, row in top_snps_pos.iterrows():
        chr1 = row['CHR']
        bp1 = row['BP']
        # Find if any known locus is on same chr and within 500kb
        close = gws_loci[(gws_loci['CHR'] == chr1) & (abs(gws_loci['BP'] - bp1) <= 500000)]
        if close.empty:
            novel_snps.append(row)

    novel_snps_df = pd.DataFrame(novel_snps)
    novel_snps_df.to_csv(os.path.join(output_dir, 'novel_top_signals.csv'), index=False)
    print(f"Checked {len(top_snps)} top signals. {len(found)} found in GWAS, {len(not_found)} not found.")
    print(f"Novel top signals (not within ±500kb of OG GWAS GWS loci): {len(novel_snps_df)}")
    print(f"Results saved in {output_dir}")

if __name__ == "__main__":
    main()
