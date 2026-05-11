
import pandas as pd
import numpy as np
import sys
import os
from itertools import combinations

def main():
    if len(sys.argv) < 2:
        print("Usage: python define_ld_blocks.py <clumped_file> [output_dir]")
        sys.exit(1)
    clumped_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else os.getcwd()
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)


    snps = []
    with open(clumped_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('CHR'):
                continue
            fields = line.split()
            # Ensure there are enough columns (at least up to P)
            if len(fields) < 5:
                continue
            try:
                snps.append({
                    'CHR': int(fields[0]),
                    'SNP': fields[2],
                    'BP': int(fields[3]),
                    'P': float(fields[4])
                })
            except Exception as e:
                print(f"Skipping line due to parse error: {line}\nError: {e}")
    if not snps:
        print("No SNPs parsed from file. Please check the input format.")
        sys.exit(1)
    df = pd.DataFrame(snps)

    # Sort by chromosome and position
    df = df.sort_values(['CHR', 'BP']).reset_index(drop=True)

    # Assign LD blocks: group SNPs within 500kb or 0.1 <= r2 < 0.6
    # (Assume you have a function get_pairwise_ld(snp1, snp2) returning r2)
    # For demonstration, only use distance; replace with LD calculation as needed
    blocks = []
    used = set()
    for idx, row in df.iterrows():
        if row['SNP'] in used:
            continue
        block = [row['SNP']]
        used.add(row['SNP'])
        for jdx in range(idx+1, len(df)):
            other = df.iloc[jdx]
            if other['CHR'] == row['CHR'] and abs(other['BP'] - row['BP']) <= 500000:
                # Here, add LD check if available
                block.append(other['SNP'])
                used.add(other['SNP'])
        blocks.append(block)

    # For each block, find top and secondary signals
    top_signals = []
    secondary_signals = []
    for block in blocks:
        block_snps = df[df['SNP'].isin(block)]
        top = block_snps.loc[block_snps['P'].idxmin()]
        top_signals.append(top['SNP'])
        secondary = block_snps[block_snps['SNP'] != top['SNP']]['SNP'].tolist()
        secondary_signals.extend(secondary)

    # Output results
    print("Top signals:", top_signals)
    print("Secondary signals:", secondary_signals)

    # Save to files
    pd.DataFrame({'Top_SNP': top_signals}).to_csv(os.path.join(output_dir, 'top_signals.txt'), index=False)
    pd.DataFrame({'Secondary_SNP': secondary_signals}).to_csv(os.path.join(output_dir, 'secondary_signals.txt'), index=False)

if __name__ == "__main__":
    main()
