#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 22:58:54 2025

@author: alexandracampbell
"""
pip install scikit-bio

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import KMeans
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import QuantileTransformer
import os

os.getcwd()
wd = "Desktop/Xesto2025/xesto_reads"
os.chdir(wd)
print(f"Working directory changed to: {wd}")

# %%



import pandas as pd
import numpy as np
from scipy.spatial.distance import pdist, squareform
from skbio.stats.ordination import pcoa
from skbio.diversity import beta_diversity
from skbio.stats.distance import permanova, DistanceMatrix
from skbio.stats.distance import permdisp
import matplotlib.pyplot as plt
import seaborn as sns
import os

# Load necessary data
asv_df = pd.read_csv("asv_table2.tsv", sep="\t", index_col=0)  # Ensure first column is sample IDs
metadata_df = pd.read_csv("xesto454_metadata_new.csv", sep=",", index_col=0)
asv_rel = asv_df.div(asv_df.sum(axis=0), axis=1) * 100  # Relative abundance
#asv_rel_t = asv_rel.T  # Samples as rows

# Match metadata and ASV table
#asv_rel_t.index.name = "SampleID"
merged_df = asv_rel.join(metadata_df)
merged_metadata = metadata_df.set_index("SampleID").loc[asv_rel.index]
grouping = merged_metadata["DiseaseState"]

# Function to process beta diversity
def beta_analysis(metric_name):
    dist_matrix = beta_diversity(metric=metric_name, counts=asv_rel_t.values, ids=asv_rel_t.index)
    dist_df = pd.DataFrame(dist_matrix.data, index=dist_matrix.ids, columns=dist_matrix.ids)
    dist_df.to_csv(f"/mnt/data/beta_{metric_name}_distance_matrix.tsv", sep="\t")

    # Ordination (PCoA)
    ordination = pcoa(dist_matrix)
    coords = ordination.samples
    coords["DiseaseState"] = grouping.values
    coords.to_csv(f"/mnt/data/beta_{metric_name}_pcoa_coords.tsv", sep="\t")

    # Plot PCoA
    plt.figure(figsize=(8, 6))
    sns.scatterplot(data=coords, x='PC1', y='PC2', hue='DiseaseState', s=100)
    plt.title(f"PCoA - {metric_name}")
    plt.xlabel(f"PC1 ({ordination.proportion_explained[0]*100:.1f}%)")
    plt.ylabel(f"PC2 ({ordination.proportion_explained[1]*100:.1f}%)")
    plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig(f"/mnt/data/beta_{metric_name}_pcoa_plot.png")
    plt.close()

    # PERMANOVA
    permanova_res = permanova(dist_matrix, grouping, permutations=999)
    permanova_df = pd.DataFrame([permanova_res.to_series()])
    permanova_df.to_csv(f"/mnt/data/beta_{metric_name}_permanova.tsv", sep="\t", index=False)

    # PERMDISP
    disp_result = permdisp(dist_matrix, grouping)
    disp_df = pd.DataFrame(disp_result.to_series()).T
    disp_df.to_csv(f"/mnt/data/beta_{metric_name}_permdisp.tsv", sep="\t", index=False)

    return f"{metric_name} analysis complete."

# Run analyses for Bray-Curtis and Jaccard
results = {
    "Bray-Curtis": beta_analysis("braycurtis"),
    "Jaccard": beta_analysis("jaccard")
}

import ace_tools as tools; tools.display_dataframe_to_user(name="Beta Diversity Analysis Results", dataframe=pd.DataFrame(results.items(), columns=["Metric", "Status"]))
