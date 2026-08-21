# scripts/fold_changes.py
# Called by Snakemake's fold_changes rule.
# Computes fold-changes between two groups, sorted by absolute FC descending.

import pandas as pd
import numpy as np

fileA       = snakemake.input.fileA
fileB       = snakemake.input.fileB
output_file = snakemake.output.fc
pseudocount = snakemake.params.pseudocount

# Parse group names from the pair wildcard (e.g. "Clone1_induced_T24_vs_WT_uninduced_T0")
pair  = snakemake.wildcards.pair
grpA, grpB = pair.split("_vs_", 1)

# -- Load mean count files --
dfA = pd.read_csv(fileA, sep="\t")
dfB = pd.read_csv(fileB, sep="\t")

# -- Merge on gene_id --
merged = dfA.merge(dfB, on=["gene_id", "description"], suffixes=(f"_{grpA}", f"_{grpB}"))

# -- Compute fold change with pseudocount --
merged["fold_change"] = (
    (merged[f"mean_count_{grpA}"] + pseudocount) /
    (merged[f"mean_count_{grpB}"] + pseudocount)
)

# -- Sort by absolute fold change descending --
merged["abs_fc"] = merged["fold_change"].abs()
merged = merged.sort_values("abs_fc", ascending=False).drop(columns="abs_fc")

merged.to_csv(output_file, sep="\t", index=False)
