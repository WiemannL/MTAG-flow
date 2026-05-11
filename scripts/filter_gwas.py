#!/usr/bin/env python3
import pandas as pd
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Filter GWAS/MTAG summary stats file for genome-wide significant hits and create a clumping input file."
    )
    parser.add_argument("input", help="Input GWAS/MTAG summary stats file (_phenotype_formatted.txt)")
    parser.add_argument("--pvalcol", default="P", help="Column name for p-value (default: P)")
    parser.add_argument("--pvalthr", type=float, default=5e-8, help="P-value threshold for significance (default: 5e-8)") ##filter for significant SNPs only 
    parser.add_argument("--out_prefix", help="Custom prefix for output files (default: input file basename without .txt)")

    args = parser.parse_args()
    input_file = args.input
    outprefix = args.out_prefix if args.out_prefix else input_file.rsplit(".", 1)[0]

    # Read data
    df = pd.read_csv(input_file, sep="\t")
    if args.pvalcol not in df.columns:
        raise ValueError(f"P-value column '{args.pvalcol}' not found in file!")

    # Filter
    df_gws = df[df[args.pvalcol] < args.pvalthr].copy()

    # Write filtered full table
    gws_full = f"{outprefix}.postMTAG.txt"
    df_gws.to_csv(gws_full, sep="\t", index=False)

    # Write clumping file
    required_cols = ["SNP", "CHR", "BP", args.pvalcol]
    missing = [c for c in required_cols if c not in df_gws.columns]
    if missing:
        raise ValueError(f"Missing columns in file for clumping: {missing}")
    gws_plink = f"{outprefix}.postMTAG.clump.tsv"
    df_gws[required_cols].to_csv(gws_plink, sep="\t", index=False)

    print(f"Wrote {len(df_gws)} genome-wide significant SNPs to {gws_full} and {gws_plink}")

if __name__ == "__main__":
    main()