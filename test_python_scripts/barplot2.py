#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 29 23:44:12 2025

@author: alexandracampbell
"""
#%%
#This part builds the relative abundance table for generating bar plots
import os


os.getcwd()
wd = "Desktop/Xesto2025/xesto_reads"
os.chdir(wd)
print(f"Working directory changed to: {wd}")


import pandas as pd



df = pd.read_csv("asv_table2.tsv", sep="\t", header=0)  # Ensure first column is sample IDs
df = df.set_index("Feature_ID")

df = df.loc[(df != 0).any(axis=1), (df != 0).any(axis=0)]


#Here is a compositional data transformation to account for enviornmental heterogeneity.
import pandas as pd


taxonomy_df = pd.read_csv("cleaned_taxa_data.tsv", sep="\t", index_col=0)  # FeatureID as index
mapfile = pd.read_csv("xesto454_metadata_new.csv", sep=",", index_col=0)
taxonomy_df.index = taxonomy_df.index.str.strip()

# Convert to relative abundance
df_rel = df.div(df.sum(axis=0), axis=1) 

df_t = df_rel.T
df_t.index.name = "SampleID"
df_with_meta = df_t.join(mapfile, how="inner") 
df_with_meta = df_with_meta.T

merged_table = pd.merge(df_with_meta, taxonomy_df, left_index=True, right_index=True, how='left')
merged_table2 = merged_table.T

valid_states = ["Control", "HoD", "Intermediate", "Diseased"]  # Adjust names to match your metadata
filtered_df_with_meta = merged_table2[merged_table2["DiseaseState"].isin(valid_states)]


#%% Generate barplots
import seaborn as sns
import colorcet as cc
import matplotlib.pyplot as plt

# Define taxonomic levels (only those that exist in your taxonomy_df)
tax_levels = ['Phylum', 'Class', 'Order', 'Family', 'Genus']
available_levels = [lvl for lvl in tax_levels if lvl in taxonomy_df.columns]
print(f"Detected taxonomy levels: {available_levels}")

glasbey_colors = cc.glasbey[:50]

for level in available_levels:
    print(f"Generating {level} plot...")

    # Start fresh each loop
    filtered_df_t = filtered_df_with_meta.T.copy()
    filtered_df_t['Feature_ID'] = filtered_df_t.index

    # Merge taxonomy for this level
    filtered_df_t = filtered_df_t.merge(taxonomy_df[[level]], left_on='Feature_ID', right_index=True)
    filtered_df_t = filtered_df_t.rename(columns={level: 'Taxon'})

    # Melt and merge with metadata
    df_long = filtered_df_t.melt(
        id_vars=['Feature_ID', 'Taxon'],
        var_name='SampleID',
        value_name='Abundance'
    )
    df_long = df_long.merge(mapfile[['DiseaseState']], left_on='SampleID', right_index=True)

    # Group and clean
    grouped_taxa = df_long.groupby(['DiseaseState', 'Taxon'], as_index=False)['Abundance'].sum()
    grouped_taxa['Abundance'] = pd.to_numeric(grouped_taxa['Abundance'], errors='coerce')
    grouped_taxa = grouped_taxa.dropna(subset=['Abundance'])
    
    #Change COntrol to Healthy
    grouped_taxa['DiseaseState'] = grouped_taxa['DiseaseState'].replace({'Control': 'Healthy'})


    # Pick top 50 for this specific level
    top_taxa = grouped_taxa.groupby('Taxon')['Abundance'].sum().nlargest(50).index
    filtered = grouped_taxa[grouped_taxa['Taxon'].isin(top_taxa)]

    # Pivot and normalize within each DiseaseState
    plot_df = filtered.pivot(index='DiseaseState', columns='Taxon', values='Abundance').fillna(0)
    plot_df = plot_df.div(plot_df.sum(axis=1), axis=0)

    # Plot
    plot_df.plot(kind='bar', stacked=True, figsize=(12, 12), color=glasbey_colors)
    plt.ylabel("Relative Abundance")
    plt.title(f"Top 50 Taxa Across Disease States ({level} Level)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left', fontsize=10, title='Taxon')
    plt.tight_layout()
    plt.savefig(f'{level}_barplot_xesto50.png', dpi=300)
    plt.close()
#%% Specific box plots of target taxa
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import kruskal
from statsmodels.stats.multitest import multipletests

# --- Define taxa of interest ---
#Change COntrol to Healthy
filtered['DiseaseState'] = filtered['DiseaseState'].replace({'Control': 'Healthy'})
target_taxa = ["Rhodobacteraceae"]

# --- Filter only taxa that actually exist ---
available_taxa = [t for t in target_taxa if t in filtered_df_t["Taxon"].unique()]
filtered = filtered_df_t[filtered_df_t["Taxon"].isin(available_taxa)].copy()
#Change COntrol to Healthy
filtered['DiseaseState'] = filtered['DiseaseState'].replace({'Control': 'Healthy'})


# --- Kruskal-Wallis + FDR correction ---
results = []
for taxon in available_taxa:
    sub = filtered[filtered["Taxon"] == taxon]
    groups = [group["Relative_Abundance"].values for _, group in sub.groupby("DiseaseState")]
    if len(groups) > 1 and all(len(g) > 1 for g in groups):
        stat, p = kruskal(*groups)
        results.append((taxon, stat, p))
    else:
        results.append((taxon, None, None))

kruskal_df = pd.DataFrame(results, columns=["Taxon", "H_stat", "p_value"])
kruskal_df["p_adj"] = multipletests(
    [p for p in kruskal_df["p_value"] if p is not None],
    method="fdr_bh"
)[1].tolist() + [None] * (len(kruskal_df) - len([p for p in kruskal_df["p_value"] if p is not None]))

filtered = filtered.merge(kruskal_df[["Taxon", "p_adj"]], on="Taxon", how="left")

# --- Faceted Plot ---
sns.set(style="whitegrid", font_scale=1.2)
g = sns.FacetGrid(filtered, col="Taxon", sharey=False, height=4, aspect=0.9)
g.map_dataframe(sns.boxplot, x="DiseaseState", y="Relative_Abundance", palette="Set3")
g.map_dataframe(sns.stripplot, x="DiseaseState", y="Relative_Abundance", color="black", size=3, jitter=True)

# Add p-values as facet titles
for ax, taxon in zip(g.axes.flat, available_taxa):
    pval = filtered.loc[filtered["Taxon"] == taxon, "p_adj"].iloc[0]
    title = f"{taxon}\nFDR p = {pval:.3e}" if pd.notnull(pval) else f"{taxon}\n(p unavailable)"
    ax.set_title(title)
    ax.set_xlabel("")
    ax.set_ylabel("Relative Abundance")

plt.tight_layout()
plt.show()
