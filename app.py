import math

import requests
import streamlit as st
import pandas as pd
import numpy as np  # It's good practice to import numpy for numerical operations
import plotly.express as px
from fastaframes import to_df
from peptacular.protein import find_peptide_indexes


# Define a function to safely calculate log2 to handle division by zero or negative numbers
def safe_log2(x):

    if x > 0:
        return math.log2(x)
    elif x == 0:
        return float('-inf')  # Return NaN if the ratio is zero or negative
    else:
        return np.nan


def ratio_inf(row):
    light_value = row['Channel.L']
    heavy_value = row['Channel.H']

    if heavy_value == 0:
        return float('inf')
    else:
        return light_value / heavy_value

def calculate_accessibility(row):
    if row['Light/Heavy.Log2Ratio'] == float('-inf'):
        return 0

    if row['Light/Heavy.Log2Ratio'] == float('inf'):
        return 100

    return 100 * 2 ** row['Light/Heavy.Log2Ratio'] / (2 ** row['Light/Heavy.Log2Ratio'] + 1)

st.title("CPP Analysis Tool")

results_file = st.file_uploader("Upload results file", type=[".tsv"])
fasta_file = st.file_uploader("Upload fasta file", type=[".fasta"])
channel_q_value_filter = st.number_input("Channel Q-value filter", value=0.10, min_value=0.0, max_value=1.0, step=0.01,
                                         help='Channel.Q.Value reflects the confidence that the precursor is indeed present in the respective channel')
remove_zeros = st.checkbox("Remove zeros", value=True, help='Remove zero values from the data')
fill_inf_with = st.number_input("Fill inf with", value=100)
merge_identical_ratios = st.checkbox("Merge identical ratios", value=True, help='Merge identicle ratios')
merge_tolerance = float(st.text_input("Merge tolerance", value='0.0001', help='Merge tolerance'))

c1, c2 = st.columns(2)
min_lysine_count = c1.number_input("Min lysine count", value=1, min_value=0, max_value=10, step=1,
                                   help='Minimum number of lysines in a peptide')
max_lysine_count = c2.number_input("Max lysine count", value=1, min_value=0, max_value=10, step=1,
                                   help='Maximum number of lysines in a peptide')

min_ms1_area = st.number_input("Min MS1 area", value=10_000, min_value=0, help='Minimum MS1 area')

c1, c2 = st.columns(2)
min_evidence_ms1 = c2.number_input("Min evidence MS1", value=0.0, min_value=0.0, max_value=1.0, help='Minimum evidence MS1')
min_evidence_ms2 = c1.number_input("Min evidence MS2", value=0.0, min_value=0.0, help='Minimum evidence MS2')

# Check if the file is uploaded
if results_file is None:
    st.warning("Please upload a file to proceed.")
    st.stop()

# load file
df = pd.read_csv(results_file, sep='\t')

# Filter based on the number of lysine's in the peptide
df = df[(df['Stripped.Sequence'].str.count('K') <= max_lysine_count) &
        (df['Stripped.Sequence'].str.count('K') >= min_lysine_count)]

# remove rows with both channels missing
df = df[(df['Channel.L'] != 0) | (df['Channel.H'] != 0)]

# Remove rows with zero values in either Channel.L or Channel.H
if remove_zeros:
    df = df[(df['Channel.L'] != 0) & (df['Channel.H'] != 0)]

# Remove rows which don't have a Channel.L or Channel.H value above the min_ms1_area
df = df[(df['Channel.L'] >= min_ms1_area) | (df['Channel.H'] >= min_ms1_area)]

# Apply channel.evidence filter
df = df[(df['Channel.Evidence.Ms1'] >= min_evidence_ms1)]
df = df[(df['Channel.Evidence.Ms2'] >= min_evidence_ms2)]

# Filter the data based on the Q-value
df = df[df['Channel.Q.Value'] <= channel_q_value_filter]

# Calculate the light/heavy ratio
df['Light/Heavy.Ratio'] = df.apply(ratio_inf, axis=1)  # Light/Heavy.Ratio in range of [0 - inf]

# Use the safe_log2 function for the 'Light/Heavy.Log2Ratio' column
df['Light/Heavy.Log2Ratio'] = df['Light/Heavy.Ratio'].apply(safe_log2)  # Light/Heavy.Log2Ratio in range of [-inf - inf]

# Calculate Accessibility from Log2Ratio, but override with 0 or 100 if
df['Accessibility'] = df.apply(calculate_accessibility, axis=1)

# replace inf values with fill_inf_with
df['Light/Heavy.Log2Ratio'] = df['Light/Heavy.Log2Ratio'].replace([np.inf], fill_inf_with)
df['Light/Heavy.Log2Ratio'] = df['Light/Heavy.Log2Ratio'].replace([-np.inf], -fill_inf_with)


if merge_identical_ratios:
    # In the dataframe df, some Stripped.Sequence's have Light/Heavy.Log2Ratio's that are nearly identical, but not exactly the same
    # For these instances, we must keep only the first occurrence of each row

    # Sort by 'Stripped.Sequence' and 'Light/Heavy.Log2Ratio' to ensure duplicates are ordered
    df.sort_values(by=['Stripped.Sequence', 'Light/Heavy.Log2Ratio'], inplace=True)

    # Use 'duplicated' to mark rows that have an identical sequence and a very close ratio as duplicates
    df['is_duplicate'] = df.duplicated(subset=['Stripped.Sequence'], keep='first') & \
                         (df.groupby('Stripped.Sequence')['Light/Heavy.Log2Ratio'].diff().abs().fillna(
                             0) < merge_tolerance)

    # Keep rows where 'is_duplicate' is False and Light/Heavy.Ratio is not 0 or inf
    df = df[(~df['is_duplicate'])].drop(columns='is_duplicate')


if fasta_file is not None:

    fasta_df = to_df(fasta_file)

    protein_name_to_sequence = {}
    for index, row in fasta_df.iterrows():
        protein_name_to_sequence[row['unique_identifier']] = row['protein_sequence']

    peptide_indexes_list = []
    site_index_list = []
    peptide_strings = []
    protein_strings = []
    for index, row in df.iterrows():
        protein_names = row['Protein.Ids'].split(';')
        stripped_sequence = row['Stripped.Sequence']
        indexes_by_protein = []
        sites_by_protein = []
        for protein_name in protein_names:
            # find index of sequence in protein sequence
            protein_sequence = protein_name_to_sequence[protein_name]
            if stripped_sequence not in protein_sequence:
                st.warning(f"Sequence {stripped_sequence} not found in protein {protein_name}")

            peptide_indexes = find_peptide_indexes(protein_sequence, stripped_sequence)
            indexes_by_protein.append([i + 1 for i in peptide_indexes])

            site_indexes_by_peptide = []
            for peptide_index in peptide_indexes:

                if 'K' not in stripped_sequence:
                    st.warning(f"Sequence {stripped_sequence} does not contain any lysine's")

                site_indexes = find_peptide_indexes(stripped_sequence, 'K')

                for site_index in site_indexes:
                    site_indexes_by_peptide.append(peptide_index + site_index + 1)

            sites_by_protein.append(site_indexes_by_peptide)

        peptide_index_str = ';'.join([','.join(map(str, x)) for x in indexes_by_protein])
        site_index_str = ';'.join([','.join(map(str, x)) for x in sites_by_protein])

        peptide_strs = [[stripped_sequence + '@' + str(index) for index in indexes] for indexes in indexes_by_protein]
        protein_strs = [[protein_name + '@' + str(index) for index in indexes] for protein_name, indexes in zip(protein_names, sites_by_protein)]

        peptide_indexes_list.append(peptide_index_str)
        site_index_list.append(site_index_str)
        peptide_strings.append(';'.join([','.join(x) for x in peptide_strs]))
        protein_strings.append(';'.join([','.join(x) for x in protein_strs]))

    df['Peptide.Indexes'] = peptide_indexes_list
    df['Site.Indexes'] = site_index_list
    df['Peptide.Index.Strings'] = peptide_strings
    df['Protein.Site.Strings'] = protein_strings

st.subheader("Filtered Data")
st.metric(label="Number of peptides", value=df.shape[0])
st.dataframe(df)

# Group by 'Stripped.Sequence' and calculate the required statistics
stats_df = df.groupby(['Stripped.Sequence', 'Protein.Ids', 'Protein.Group', 'Protein.Names', 'Genes'])['Light/Heavy.Log2Ratio'].agg(
    ['mean', 'std', 'count', 'sem', 'median', 'min', 'max'])
stats_df.reset_index(inplace=True)

# Rename the columns for clarity
stats_df.columns = [
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

st.subheader("Peptide level Statistics")
stats_df['Accessibility'] = 100 * 2 ** stats_df['Log2Ratio.Mean'] / (2 ** stats_df['Log2Ratio.Mean'] + 1)

st.dataframe(stats_df)

# Create the scatter plot with mean on the x-axis and standard deviation on the y-axis
fig = px.scatter(stats_df, x='Log2Ratio.Mean', y='Log2Ratio.SEM', hover_name='Stripped.Sequence',
                 hover_data=['Protein.Names', 'Log2Ratio.Count', 'Log2Ratio.SEM', 'Log2Ratio.Median',
                             'Log2Ratio.Min', 'Log2Ratio.Max'])

# Enhance the plot with titles and labels
fig.update_traces(textposition='top center')
fig.update_layout(
    title='Mean vs. SEM (per Sequence)',
    xaxis_title='Mean of Log2 Ratios',
    yaxis_title='SEM of Log2 Ratios',
    showlegend=False
)

st.plotly_chart(fig)

if fasta_file is not None:
    # split by protein_str
    data = []
    for index, row in df.iterrows():
        protein_site_strs = []
        unique = len(row['Protein.Site.Strings'].split(';')) == 1
        for protein_strs in row['Protein.Site.Strings'].split(';'):
            for protein_str in protein_strs.split(','):
                protein_site_strs.append(protein_str)

        for protein_site_str in protein_site_strs:
            protein_name, site_index = protein_site_str.split('@')
            data.append([protein_name, int(site_index), row['Light/Heavy.Log2Ratio'], int(unique)])

    protein_site_df = pd.DataFrame(data, columns=['Protein.Id', 'Site.Index', 'Log2Ratio', 'Unique'])

    #st.subheader("Site level Data")
    #st.dataframe(protein_site_df)

    # Group by 'Stripped.Sequence' and calculate the required statistics
    protein_site_df = protein_site_df.groupby(['Protein.Id', 'Site.Index'])[
        'Log2Ratio'].agg(['mean', 'std', 'count', 'sem', 'median', 'min', 'max'])
    protein_site_df.reset_index(inplace=True)

    # Rename the columns for clarity
    protein_site_df.columns = [
        'Protein.Id',
        'Site.Index',
        'Log2Ratio.Mean',
        'Log2Ratio.Std',
        'Log2Ratio.Count',
        'Log2Ratio.SEM',
        'Log2Ratio.Median',
        'Log2Ratio.Min',
        'Log2Ratio.Max'
    ]

    protein_site_df['Protein.Site.String'] = protein_site_df['Protein.Id'] + '@' + protein_site_df['Site.Index'].astype(str)

    protein_site_df['Accessibility'] = 100 * 2 ** protein_site_df['Log2Ratio.Mean'] / (2 ** protein_site_df['Log2Ratio.Mean'] + 1)

    st.subheader("Site level Statistics")
    st.dataframe(protein_site_df)

    # Create the scatter plot with mean on the x-axis and standard deviation on the y-axis
    fig = px.scatter(protein_site_df, x='Log2Ratio.Mean', y='Log2Ratio.SEM', hover_name='Protein.Site.String',
                     hover_data=['Protein.Id', 'Site.Index', 'Log2Ratio.Count', 'Log2Ratio.SEM', 'Log2Ratio.Median',
                                 'Log2Ratio.Min', 'Log2Ratio.Max'])

    # Enhance the plot with titles and labels
    fig.update_traces(textposition='top center')
    fig.update_layout(
        title='Mean vs. SEM (per Site)',
        xaxis_title='Mean of Log2 Ratios',
        yaxis_title='SEM of Log2 Ratios',
        showlegend=False
    )

    st.plotly_chart(fig)
