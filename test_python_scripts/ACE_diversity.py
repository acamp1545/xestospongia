#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 22 23:50:14 2025

@author: alexandracampbell
"""

import numpy as np
import pandas as pd

def calculate_ace(sample_counts):
    sample_counts = np.array(sample_counts)
    sample_counts = sample_counts[sample_counts > 0]

    if len(sample_counts) == 0:
        return np.nan

    rare = sample_counts[sample_counts <= 10]
    abund = sample_counts[sample_counts > 10]

    S_rare = len(rare)
    S_abund = len(abund)
    N_rare = rare.sum()

    F1 = np.sum(rare == 1)

    if N_rare == 0 or S_rare == 0:
        return S_abund  # fallback if no rare species

    C_ace = 1 - (F1 / N_rare)

    if C_ace == 0:
        return np.nan

    freqs = pd.Series(rare).value_counts()
    i_vals = freqs.index
    fi = freqs.values
    sum_fi_i = np.sum(i_vals * fi)
    sum_fi_i_i1 = np.sum(i_vals * (i_vals - 1) * fi)

    gamma_sq = max((S_rare * sum_fi_i_i1 / (sum_fi_i**2)) - 1, 0)

    ace = S_abund + (S_rare / C_ace) + (F1 / C_ace) * gamma_sq
    return ace

ace_scores = counts.apply(calculate_ace, axis=1)
div_df["ACE"] = ace_scores
