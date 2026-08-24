# -*- coding: utf-8 -*-
"""


"""
# %%
import os
os.chdir("/Users/alexandracampbell/Desktop/Xesto2025/xesto_reads/ML_models")
os.getcwd()
# %%

#%%
import pandas as pd
df = pd.read_csv ("extracted_alpha.tsv", sep = ',')
print(df)

print(df.loc[:,'SampleID'])
#%%

# %%
# Import necessary libraries
import matplotlib.pyplot as plt 
import seaborn as sns



# Basic plotting with Matplotlib
plt.figure(figsize=(8, 6))
plt.plot(df['Type'], df['shannon_entropy'], label='Shannon Entropy', marker='o')
plt.plot(df['Type'], df['ace'], label='ACE Diversity', marker='o')
plt.plot(df['Type'], df['faith_pd'], label="Faith PD", marker='o')
plt.plot(df['Type'], df['observed_features'], label="Observed Features", marker='o')
plt.title('Diversity Measures per Sample')
plt.xlabel('SampleID')
plt.ylabel('Diversity Index')
plt.legend()
plt.xticks(rotation=45)
plt.show()

# Seaborn for a nicer plot
sns.set(style="whitegrid")
plt.figure(figsize=(8, 6))
sns.barplot(x='SampleID', y='shannon_entropy', data=df, color='blue', label='Shannon Entropy')
sns.barplot(x='SampleID', y='ace', data=df, color='red', label='ACE Diversity')
plt.title('Diversity Measures per Sample')
plt.ylabel('Diversity Index')
plt.xticks(rotation=45)
plt.legend()
plt.show()
# %%

# %%
from scipy import stats
# Example of comparing three or more groups
groups = [df[df['Type'] == group]['shannon_entropy'] for group in df['Type'].unique()]    #This is an ANOVA

# Perform One-Way ANOVA
f_stat, p_value = stats.f_oneway(*groups)
print(f"F-statistic: {f_stat}, p-value: {p_value}")
# %%

# %%
from scipy.stats import linregress

# Example: Predict Shannon diversity based on Observed_Features
slope, intercept, r_value, p_value, std_err = linregress(df['faith_pd'], df['shannon_entropy'])

print(f"Slope: {slope}, Intercept: {intercept}")
print(f"R-squared: {r_value**2}, p-value: {p_value}")
# %%

# %%
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import kruskal

# Boxplot comparing Shannon diversity across different types
#sns.boxplot(x='Type', y='shannon_entropy', data=df)+ plt.title('Shannon Entropy by Type')
#plt.show()

# Boxplot comparing Shannon diversity across different types
#sns.boxplot(x='Type', y='faith_pd', data=df)+plt.title('Faith PD by Type')
#plt.show()

# %%
# Boxplot comparing Ace diversity across different types
sns.set_style("white")
sns.boxplot(x='Type', y='ace', data=df), plt.title('ACE Diversity by Type'), plt.yscale("log"),sns.color_palette("colorblind")
plt.show()

# %%
plt.figure(figsize=(8,6))
ax = sns.boxplot(data=df, x="Type", y="ace",
                 width=0.5,  # Controls the width of the boxes
                 boxprops={'facecolor':'white', 'edgecolor':'black', 'linewidth':1.5},  # White fill, black border
                 whiskerprops={'color':'black', 'linewidth':1.5},  # Black whiskers
                 capprops={'color':'black', 'linewidth':1.5},  # Black caps
                 medianprops={'color':'red', 'linewidth':2})  # Red median line

# %%
from scipy.stats import shapiro, ks_2samp

# Example: Checking normality for 'ace' column
stat, p = shapiro(df["ace"])  # Shapiro-Wilk Test
print(f"Shapiro-Wilk Test: p-value = {p}") #p-value = 3.4767473462583266e-10, not normal

# If p < 0.05, the data is **not normally distributed**

# %%
from scipy.stats import kruskal

stat, p = kruskal(*groups)
print(f"Kruskal-Wallis Test: p-value = {p}")

# %%
import seaborn as sns
import matplotlib.pyplot as plt

sns.boxplot(data=df, x="Type", y="ace")
plt.show()

print(df.dtypes)
# %%

from sklearn.utils import resample #This resamples data to hopefully increase statistical power
df_balanced = df.groupby("Type", group_keys=False).apply(lambda x: resample(x, n_samples=min(df["Type"].value_counts()), random_state=42))

# %%
from skbio.stats.distance import DistanceMatrix
from scipy.spatial.distance import cdist
from scipy.spatial.distance import pdist
df2 = df.apply(pd.to_numeric, errors='coerce')

euclidean_dist = pdist(df2, metric='euclidean')
distance_matrix = DistanceMatrix(euclidean_dist, ids=df2.index)
print(distance_matrix)

# %%
df2 = df2.dropna()

print(df2.dtypes)

euclidean_dist = pdist(df2, metric='euclidean')
distance_matrix = DistanceMatrix(euclidean_dist, ids=df2.index)

if df2.empty:
    print("The DataFrame is empty.")
