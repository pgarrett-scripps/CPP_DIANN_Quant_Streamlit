import pandas as pd
import numpy as np
import math
import peptacular as pt


def ratio_inf(row):
    light_value = row['Ms1.Area.L']
    heavy_value = row['Ms1.Area.H']

    if pd.isna(heavy_value):  # Proper NaN check
        return np.inf
    elif pd.isna(light_value):  # Proper NaN check
        return 0
    else:
        return light_value / heavy_value


def safe_log2(x):
    if np.isinf(x):
        return np.inf
    elif x == 0:
        return -np.inf
    elif x > 0:
        return math.log2(x)
    else:
        st.warning(f"Invalid value for log2: {x}")
        return np.nan


def calculate_accessibility(row):
    log2ratio = row['Light/Heavy.Log2Ratio']

    if np.isneginf(log2ratio):
        return 0

    if np.isinf(log2ratio):
        return 100

    return 100 * 2 ** log2ratio / (2 ** log2ratio + 1)


def cached_find_peptide_indexes(protein_sequence, stripped_sequence):
    return pt.find_subsequence_indices(protein_sequence, stripped_sequence)


def cached_find_site_indexes(stripped_sequence, site):
    return pt.find_subsequence_indices(stripped_sequence, site)
