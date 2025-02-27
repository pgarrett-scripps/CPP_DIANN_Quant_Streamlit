import math

import streamlit as st
import pandas as pd
import numpy as np
import plotly.express as px
from fastaframes import to_df
from peptacular.protein import find_peptide_indexes
from scipy.stats import norm


# Define a function to safely calculate log2 to handle division by zero or negative numbers


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


with st.sidebar:
    st.title("CPP Analysis Tool")

    results_file = st.file_uploader("Upload results file", type=[".tsv", ".zip", ".parquet"],
                                    accept_multiple_files=True)
    fasta_file = st.file_uploader("Upload fasta file", type=[".fasta"])
    channel_q_value_filter = st.number_input("Channel Q-value filter", value=0.10, min_value=0.0, max_value=1.0,
                                             step=0.01,
                                             help='Channel.Q.Value reflects the confidence that the precursor is indeed'
                                                  ' present in the respective channel')

    remove_zeros = st.checkbox("Remove missing values", value=True,
                               help='Remove data which have either a missing L or H channel value')
    fill_inf_with = st.number_input("Fill inf with", value=100)

    c1, c2 = st.columns(2)
    min_lysine_count = c1.number_input("Min lysine count", value=1, min_value=0, max_value=10, step=1,
                                       help='Minimum number of lysines in a peptide')
    max_lysine_count = c2.number_input("Max lysine count", value=1, min_value=0, max_value=10, step=1,
                                       help='Maximum number of lysines in a peptide')

    min_ms1_area = st.number_input("Min MS1 area", value=10_000, min_value=0, help='Minimum MS1 area')

    min_evidence = st.number_input("Min Evidence", value=0.0, min_value=0.0, max_value=1.0,
                                   help='Minimum Evidence')
    min_channel_evidence = st.number_input("Min Channel Evidence", value=0.0, min_value=0.0, max_value=1.0,
                                           help='Minimum Channel Evidence')
    max_channel_q_value = st.number_input("Max Channel Q-value", value=0.05, min_value=0.0, max_value=1.0,
                                          help='Maximum Channel Q-value')

    filter_unique_sites = st.checkbox("Filter for unique sites", value=False, help='Filter for unique sites')

    # peptide filter:
    should_filter_peptides = st.checkbox("Filter peptides", value=False)

    filter_peptides = set()
    if should_filter_peptides:
        filter_peptides = st.text_input("Stripped peptides to keep (Comma seperated)", value="")
        filter_peptides = set(filter_peptides.split(','))

    should_filter_sites = st.checkbox("Filter sites", value=False)

    filter_sites = set()
    if should_filter_sites:
        filter_sites = st.text_input("Sites to keep (Comma seperated)", value="")
        filter_sites = set(filter_sites.split(','))

    replace_missing_channel = None
    if st.checkbox("Replace Missing Channel with Value"):
        replace_missing_channel = st.number_input("Replace Missing Channel with Value", value=10)

    run_btn = st.button('Run', use_container_width=True, type='primary')

if not run_btn:
    st.warning('Click the "Run" button to start the analysis.')
    st.stop()

# Check if the file is uploaded
if results_file is None or len(results_file) == 0:
    st.warning("Please upload a file to proceed.")
    st.stop()

# load file
dfs = []
for results_file in results_file:
    if results_file.name.endswith('.zip'):
        df = pd.read_csv(results_file, sep='\t', compression='zip')
    else:
        try:
            df = pd.read_csv(results_file, sep='\t')
        except:
            df = pd.read_parquet(results_file)

    dfs.append(df)

df = pd.concat(dfs)

# drop rows that dont have a Channel value in [light_channel_label, heavy_channel_label]
df = df[df['Channel'].isin(['L', 'H'])]

# NEW: Pivot data if using new long-format with a 'channel' column
index_cols = ['Run', 'Precursor.Id', 'Modified.Sequence', 'Stripped.Sequence', 'Precursor.Charge',
              'Precursor.Lib.Index', 'Decoy', 'Proteotypic', 'Protein.Ids', 'Protein.Group',
              'Protein.Names', 'Genes']

ratio_df = df.pivot_table(
    index=index_cols,
    columns=['Channel'],
    values=['Ms1.Area', 'Evidence', 'Mass.Evidence', 'Channel.Evidence', 'Q.Value', 'Global.Q.Value',
            'Channel.Q.Value'],
    aggfunc=lambda x: x if len(x) == 1 else (_ for _ in ()).throw(ValueError(f"Duplicate entries found for {x.index}"))
).reset_index()

# Collapse multi-index column names, but keep index columns unchanged
ratio_df.columns = [
    col[0] if col[1] == '' else f"{col[0]}.{col[1]}"
    for col in ratio_df.columns.to_flat_index()
]

# Remove rows with both channels missing
ratio_df.dropna(subset=['Ms1.Area.L', 'Ms1.Area.H'], how='all', inplace=True)

# filter based on the number of lysine's in the peptide
ratio_df['K_count'] = ratio_df['Stripped.Sequence'].str.count('K')
ratio_df = ratio_df[(ratio_df['K_count'] <= max_lysine_count) & (ratio_df['K_count'] >= min_lysine_count)]

if len(filter_peptides) > 0:
    ratio_df = ratio_df[ratio_df['Stripped.Sequence'].isin(filter_peptides)]

# replace 0 with NaN in Ms1.Area.L and Ms1.Area.H
ratio_df['Ms1.Area.L'] = ratio_df['Ms1.Area.L'].replace(0, np.nan)
ratio_df['Ms1.Area.H'] = ratio_df['Ms1.Area.H'].replace(0, np.nan)

# Remove rows with zero values in either Channel.L or Channel.H if required
if remove_zeros:
    ratio_df.dropna(subset=['Ms1.Area.L', 'Ms1.Area.H'], inplace=True)

# Ensure NaN values are replaced with 0 before filtering
ratio_df = ratio_df[(ratio_df['Ms1.Area.L'].fillna(0) >= min_ms1_area) |
                    (ratio_df['Ms1.Area.H'].fillna(0) >= min_ms1_area)]

# Replace missing values in Ms1.Area.L and Ms1.Area.H with the specified value
if replace_missing_channel is not None:
    ratio_df['Ms1.Area.L'].fillna(replace_missing_channel, inplace=True)
    ratio_df['Ms1.Area.H'].fillna(replace_missing_channel, inplace=True)

# Calculate the light/heavy ratio
ratio_df['Light/Heavy.Ratio'] = ratio_df.apply(ratio_inf, axis=1)  # Light/Heavy.Ratio in range of [0 - inf]

# Use the safe_log2 function for the 'Light/Heavy.Log2Ratio' column
ratio_df['Light/Heavy.Log2Ratio'] = ratio_df['Light/Heavy.Ratio'].apply(
    safe_log2)  # Light/Heavy.Log2Ratio in range of [-inf - inf]

# Calculate Accessibility from Log2Ratio, but override with 0 or 100 if
ratio_df['Accessibility'] = ratio_df.apply(calculate_accessibility, axis=1)

# replace inf values with fill_inf_with
ratio_df['Light/Heavy.Log2Ratio'] = ratio_df['Light/Heavy.Log2Ratio'].replace([np.inf], fill_inf_with)
ratio_df['Light/Heavy.Log2Ratio'] = ratio_df['Light/Heavy.Log2Ratio'].replace([-np.inf], -fill_inf_with)

# check if Evidence.L or Evidence.H is greater than min_evidence
ratio_df = ratio_df[(ratio_df['Evidence.L'] >= min_evidence) | (ratio_df['Evidence.H'] >= min_evidence)]

# check if Channel.Evidence.L or H is greater than min_evidence
ratio_df = ratio_df[
    (ratio_df['Channel.Evidence.L'] >= min_channel_evidence) | (ratio_df['Channel.Evidence.H'] >= min_channel_evidence)]

# check if Q.Value.L or H is less than max_q_value
ratio_df = ratio_df[(ratio_df['Q.Value.L'] <= max_channel_q_value) | (ratio_df['Q.Value.H'] <= max_channel_q_value)]

# check if Ms1.Area.L or H is greater than min_ms1_area
ratio_df = ratio_df[(ratio_df['Ms1.Area.L'] >= min_ms1_area) | (ratio_df['Ms1.Area.H'] >= min_ms1_area)]


def cached_find_peptide_indexes(protein_sequence, stripped_sequence):
    return find_peptide_indexes(protein_sequence, stripped_sequence)


# Cache the site index lookup function

def cached_find_site_indexes(stripped_sequence, site):
    return find_peptide_indexes(stripped_sequence, site)


with st.expander("Warnings", expanded=True):
    if fasta_file is not None:

        fasta_df = to_df(fasta_file)

        protein_name_to_sequence = {}
        for index, row in fasta_df.iterrows():
            protein_name_to_sequence[row['unique_identifier']] = row['protein_sequence']

        peptide_indexes_list = []
        site_index_list = []
        peptide_strings = []
        protein_strings = []
        for index, row in ratio_df.iterrows():
            protein_names = row['Protein.Ids'].split(';')
            stripped_sequence = row['Stripped.Sequence']
            indexes_by_protein = []
            sites_by_protein = []
            for protein_name in protein_names:
                # find index of sequence in protein sequence
                protein_sequence = protein_name_to_sequence[protein_name]
                if stripped_sequence not in protein_sequence:
                    st.warning(f"Sequence {stripped_sequence} not found in protein {protein_name}")

                peptide_indexes = cached_find_peptide_indexes(protein_sequence, stripped_sequence)
                indexes_by_protein.append([i + 1 for i in peptide_indexes])

                site_indexes_by_peptide = []
                for peptide_index in peptide_indexes:

                    if 'K' not in stripped_sequence:
                        st.warning(f"Sequence {stripped_sequence} does not contain any lysine's")

                    site_indexes = cached_find_site_indexes(stripped_sequence, 'K')

                    for site_index in site_indexes:
                        site_indexes_by_peptide.append(peptide_index + site_index + 1)

                sites_by_protein.append(site_indexes_by_peptide)

            peptide_index_str = ';'.join([';'.join(map(str, x)) for x in indexes_by_protein])
            site_index_str = ';'.join([';'.join(map(str, x)) for x in sites_by_protein])

            peptide_strs = [[stripped_sequence + '@' + str(index) for index in indexes] for indexes in
                            indexes_by_protein]
            protein_strs = [[protein_name + '@' + str(index) for index in indexes] for protein_name, indexes in
                            zip(protein_names, sites_by_protein)]

            peptide_indexes_list.append(peptide_index_str)
            site_index_list.append(site_index_str)
            peptide_strings.append(';'.join([';'.join(x) for x in peptide_strs]))
            protein_strings.append(';'.join([';'.join(x) for x in protein_strs]))

        ratio_df['Peptide.Indexes'] = peptide_indexes_list
        ratio_df['Site.Indexes'] = site_index_list
        ratio_df['Peptide.Index.Strings'] = peptide_strings
        ratio_df['Protein.Site.Strings'] = protein_strings

        if should_filter_sites:
            ratio_df = ratio_df[ratio_df['Protein.Site.Strings'].isin(filter_sites)]

        # if it contains ;
        if filter_unique_sites:
            ratio_df = ratio_df[~ratio_df['Protein.Site.Strings'].str.contains(';')]

# Group by 'Stripped.Sequence' and calculate the required statistics
peptide_stats_df = \
    ratio_df.groupby(['Run', 'Stripped.Sequence', 'Protein.Ids', 'Protein.Group', 'Protein.Names', 'Genes'])[
        'Light/Heavy.Log2Ratio'].agg(
        ['mean', 'std', 'count', 'sem', 'median', 'min', 'max'])
peptide_stats_df.reset_index(inplace=True)

# Rename the columns for clarity
peptide_stats_df.columns = [
    'Run',
    'Stripped.Sequence',
    'Protein.Ids',
    'Protein.Group',
    'Protein.Names',
    'Genes',
    'Log2Ratio.Mean',
    'Log2Ratio.Std',
    'Log2Ratio.Count',
    'Log2Ratio.SEM',
    'Log2Ratio.Median',
    'Log2Ratio.Min',
    'Log2Ratio.Max'
]

peptide_stats_df['Accessibility'] = 100 * 2 ** peptide_stats_df['Log2Ratio.Mean'] / (
        2 ** peptide_stats_df['Log2Ratio.Mean'] + 1)
peptide_stats_df['Accessibility.Min'] = 100 * 2 ** peptide_stats_df['Log2Ratio.Min'] / (
        2 ** peptide_stats_df['Log2Ratio.Min'] + 1)
peptide_stats_df['Accessibility.Max'] = 100 * 2 ** peptide_stats_df['Log2Ratio.Max'] / (
        2 ** peptide_stats_df['Log2Ratio.Max'] + 1)

# Create the scatter plot with mean on the x-axis and standard deviation on the y-axis
fig = px.scatter(peptide_stats_df, x='Log2Ratio.Mean', y='Log2Ratio.SEM', hover_name='Stripped.Sequence',
                 hover_data=['Run', 'Protein.Names', 'Log2Ratio.Count', 'Log2Ratio.SEM', 'Log2Ratio.Median',
                             'Log2Ratio.Min', 'Log2Ratio.Max'])

# Enhance the plot with titles and labels
fig.update_traces(textposition='top center')
fig.update_layout(
    title='Mean vs. SEM (per Sequence)',
    xaxis_title='Mean of Log2 Ratios',
    yaxis_title='SEM of Log2 Ratios',
    showlegend=False
)

#st.plotly_chart(fig)


# filter based on run and protein site string
site_stats_df = \
    ratio_df.groupby(['Run', 'Protein.Site.Strings', 'Protein.Ids', 'Protein.Group', 'Protein.Names', 'Genes'])[
        'Light/Heavy.Log2Ratio'].agg(
        ['mean', 'std', 'count', 'sem', 'median', 'min', 'max'])
site_stats_df.reset_index(inplace=True)

# Rename the columns for clarity
site_stats_df.columns = [
    'Run',
    'Protein.Site.Strings',
    'Protein.Ids',
    'Protein.Group',
    'Protein.Names',
    'Genes',
    'Log2Ratio.Mean',
    'Log2Ratio.Std',
    'Log2Ratio.Count',
    'Log2Ratio.SEM',
    'Log2Ratio.Median',
    'Log2Ratio.Min',
    'Log2Ratio.Max'
]

site_stats_df['Accessibility'] = 100 * 2 ** site_stats_df['Log2Ratio.Mean'] / (2 ** site_stats_df['Log2Ratio.Mean'] + 1)
site_stats_df['Accessibility.Min'] = 100 * 2 ** site_stats_df['Log2Ratio.Min'] / (
        2 ** site_stats_df['Log2Ratio.Min'] + 1)
site_stats_df['Accessibility.Max'] = 100 * 2 ** site_stats_df['Log2Ratio.Max'] / (
        2 ** site_stats_df['Log2Ratio.Max'] + 1)

st.subheader("Filtered Data")
st.metric(label="Number of peptides", value=ratio_df.shape[0])

st.subheader("Raw Data")
st.dataframe(ratio_df)

# download button
st.download_button(label="Download Raw Data",
                   data=ratio_df.to_csv(index=False).encode('utf-8'),
                   file_name='data.csv',
                   mime='text/csv',
                   use_container_width=True,
                   type='primary')

st.subheader("Peptide Data")
st.dataframe(peptide_stats_df)

st.download_button(label="Download Peptide Data",
                   data=peptide_stats_df.to_csv(index=False).encode('utf-8'),
                   file_name='peptide_stats.csv',
                   mime='text/csv',
                   use_container_width=True,
                   type='primary')

st.subheader("Site Data")
st.dataframe(site_stats_df)

st.download_button(label="Download Site Data",
                   data=site_stats_df.to_csv(index=False).encode('utf-8'),
                   file_name='site_stats.csv',
                   mime='text/csv',
                   use_container_width=True,
                   type='primary')
