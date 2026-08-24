#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 22:07:44 2025

@author: alexandracampbell
"""

import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# --- Load your data ---
asv_df = pd.read_csv("asv_table2.tsv", sep="\t", index_col=0)  # ASVs as rows, samples as columns
metadata_df = pd.read_csv("xesto454_metadata_new.csv")
taxonomy_df = pd.read_csv("taxonomy.tsv", sep="\t", index_col=0)

# 1. Load the ASV table (asv_df) and taxonomy table (taxonomy_df)
asv_df = pd.read_csv("asv.csv", index_col=0)  # ASVs x samples
taxonomy_df = pd.read_csv("cleaned_taxa_data.tsv", sep="\t", index_col=0)  # FeatureID as index

# 2. Transpose ASV table to samples x ASVs
asv_t = asv_df.T

#Test for normal distribution-Shapiro-Wilk and Equal Variance Test
from scipy.stats import shapiro, levene

# Check assumptions for each metric
for metric in ["shannon", "simpson", "chao1", "observed_otus"]:
    print(f"\nChecking assumptions for: {metric}")
    
    # Shapiro-Wilk normality test per group
    for name, group in div_df.groupby("DiseaseState"):
        stat, p = shapiro(group[metric])
        print(f"  {name} - Shapiro p={p:.4f}")

    # Levene's test for homogeneity of variances
    groups = [group[metric].values for name, group in div_df.groupby("DiseaseState")]
    lev_stat, lev_p = levene(*groups)
    print(f"  Levene’s test p={lev_p:.4f}")


# 3. Join with taxonomy (phylum level)
asv_t_with_phylum = asv_t.T.join(taxonomy_df["Phylum"]).T
asv_phylum = asv_t_with_phylum.T.groupby("Phylum").sum().T
asv_phylum.index.name = "SampleID"




##TToubleshoot
# Collapse to phylum level
import os
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

# Assumes `asv_phylum`, `metadata_df` already loaded

# Define your comparison pairs
comparison_pairs = [("Control", "Band"), ("Control", "Diseased"), ("Band", "Diseased"),("HoD", "Diseased"),("Control", "HoD"),("Band", "HoD")]

# Coerce sample ID to index if needed
asv_phylum.index.name = "SampleID"
asv_phylum = asv_phylum.reset_index()

for group1, group2 in comparison_pairs:
    print(f"\n🔍 Comparing {group1} vs {group2}")

    # Merge phylum table with metadata
    merged = asv_phylum.merge(metadata_df, on="SampleID")
    filtered = merged[merged["DiseaseState"].isin([group1, group2])]

    # Melt long
    melted = filtered.melt(id_vars=["SampleID", "DiseaseState"], 
                           var_name="Phylum", value_name="Abundance")
    melted["Abundance"] = pd.to_numeric(melted["Abundance"], errors="coerce")

    # Calculate top phyla based on overall abundance
    top_phyla = melted.groupby("Phylum")["Abundance"].sum().sort_values(ascending=False).head(10).index

    # Create output directory
    outdir = f"boxplots_{group1}_vs_{group2}"
    os.makedirs(outdir, exist_ok=True)

    # Plot each top phylum
    for phylum in top_phyla:
        data = melted[(melted["Phylum"] == phylum) & (melted["Abundance"] > 0)]

        # Ensure there's enough data to make a boxplot
        if data["DiseaseState"].nunique() < 2 or data["Abundance"].sum() == 0:
            print(f"⚠️ Skipping {phylum} — insufficient data to plot.")
            continue

        fig, ax = plt.subplots(figsize=(6, 4))
        sns.boxplot(x="DiseaseState", y="Abundance", data=data, ax=ax)
        ax.set_title(f"{phylum}: {group1} vs {group2}")
        fig.tight_layout()

        # Sanitize filename
        safe_name = phylum.replace(" ", "_").replace("/", "-").replace("__", "_").strip("_")
        outpath = f"{outdir}/boxplot_{safe_name}.png"
        fig.savefig(outpath, dpi=300)
        plt.close(fig)
        print(f"✅ Saved: {outpath}")


    
##ENd troubleshoot



#Now for the Mann-Whitney
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import os

# Load data
asv_df = pd.read_csv("asv_table2.tsv", sep="\t", index_col=0)  # ASVs as rows, samples as columns
metadata_df = pd.read_csv("xesto454_metadata_new.csv")
taxonomy_df = pd.read_csv("cleaned_taxa_data.tsv", sep="\t", index_col=0)  # FeatureID as index
asv_rel = asv_df.div(asv_df.sum(axis=0), axis=1) * 100  # % relative abundance per sample
asv_df=asv_rel

# Transpose ASV table to samples x ASVs
asv_t = asv_df.T
asv_t.index.name = "SampleID"

# Merge ASV data with metadata
merged = asv_t.merge(metadata_df, left_index=True, right_on="SampleID")

# Taxonomic levels to group by
tax_levels = ["Phylum", "Class","Order","Family"]

# Define comparisons
comparison_pairs = [("Control", "Band"), ("Control", "Diseased"), ("Band", "Diseased"),("HoD", "Diseased"),("HoD", "Control"),("HoD", "Band")]

# Loop through each taxonomic level
for level in tax_levels:
    print(f"\n🔁 Aggregating at the {level} level")

    # Join taxonomy to ASV table
    if level not in taxonomy_df.columns:
        print(f"⚠️ Taxonomic level '{level}' not found in taxonomy table. Skipping.")
        continue

    # Add taxonomy to ASV matrix (samples x ASVs → ASVs x samples to join)
    asv_with_tax = asv_t.T.join(taxonomy_df[level])
    grouped = asv_with_tax.groupby(level).sum().T  # Group by phylum or order

    # Merge with metadata again
    grouped["SampleID"] = grouped.index
    grouped = grouped.merge(metadata_df, on="SampleID")

    # Loop through group comparisons
    for group1, group2 in comparison_pairs:
        print(f"🔬 Comparing {group1} vs {group2} at {level} level")

        filtered = grouped[grouped["DiseaseState"].isin([group1, group2])]
        non_taxa_cols = ["SampleID", "DiseaseState", "depth", "latitude", "longitude", "location", "description"]
        tax_data = filtered.drop(columns=[col for col in filtered.columns if col in non_taxa_cols])
        group_labels = filtered["DiseaseState"]

        feature_names = []
        U_stats = []
        p_values = []

        for feature in tax_data.columns:
            g1 = tax_data[feature][group_labels == group1]
            g2 = tax_data[feature][group_labels == group2]

            # Skip feature if both groups are all zero
            if g1.sum() == 0 and g2.sum() == 0:
                continue

            try:
                U, p = mannwhitneyu(g1, g2, alternative='two-sided')
                if np.isfinite(p):
                    feature_names.append(feature)
                    U_stats.append(U)
                    p_values.append(p)
            except Exception as e:
                print(f"⚠️ Skipping {feature} due to error: {e}")
                continue

        if len(p_values) == 0:
            print(f"⚠️ No testable taxa for {group1} vs {group2} at {level} level.")
            continue

        # FDR correction
        reject, pvals_corrected, _, _ = multipletests(p_values, method='fdr_bh')

        # Save results
        result_df = pd.DataFrame({
            level: feature_names,
            "U_statistic": U_stats,
            "p_value": p_values,
            "p_adjusted": pvals_corrected,
            "significant": reject
        }).set_index(level).sort_values("p_adjusted")

        outfile = f"mannwhitney_{group1}_vs_{group2}_{level}.tsv"
        result_df.to_csv(outfile, sep="\t")
        if reject.any():
            sig_only = result_df[result_df["significant"]]
            sig_only["Comparison"] = f"{group1}_vs_{group2}"
            sig_only["Level"] = level
            if 'all_significant' not in locals():
                all_significant = sig_only
            else:
                all_significant = pd.concat([all_significant, sig_only], axis=0)
        print(f"📁 Saved results to: {outfile}")
        
        
if reject.any():
    sig_only = result_df[result_df["significant"]]
    sig_only["Comparison"] = f"{group1}_vs_{group2}"
    sig_only["Level"] = level
    if 'all_significant' not in locals():
        all_significant = sig_only
    else:
        all_significant = pd.concat([all_significant, sig_only], axis=0)

if 'all_significant' in locals():
    all_significant.to_csv("summary_significant_taxa.tsv", sep="\t")
    print("📊 Combined summary of all significant taxa saved to: summary_significant_taxa.tsv")
