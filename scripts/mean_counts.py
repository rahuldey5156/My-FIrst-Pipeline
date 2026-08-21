# scripts/mean_counts.py
# Called by Snakemake's mean_counts rule.
# Computes per-gene mean counts across replicates and joins with descriptions.

import pandas as pd

# snakemake object is injected automatically by Snakemake
count_files = snakemake.input.counts
bed_file    = snakemake.input.bed
output_file = snakemake.output.means

# -- Load and average count files --
dfs = []
for f in count_files:
    df = pd.read_csv(f, sep="\t", header=None, names=["gene_id", "count"])
    dfs.append(df.set_index("gene_id")["count"])

means = pd.concat(dfs, axis=1).mean(axis=1).reset_index()
means.columns = ["gene_id", "mean_count"]

# -- Extract gene descriptions from BED file (col 4 = id, col 5 = description) --
bed = pd.read_csv(bed_file, sep="\t", header=None, usecols=[3, 4])
bed.columns = ["gene_id", "description"]
bed = bed.drop_duplicates(subset="gene_id")

# -- Join means with descriptions --
result = means.merge(bed, on="gene_id", how="left")
result["description"] = result["description"].fillna("no_description")

result.to_csv(output_file, sep="\t", index=False)
