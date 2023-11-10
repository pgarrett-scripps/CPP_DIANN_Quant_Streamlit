import math

import requests
import streamlit as st
import pandas as pd
import numpy as np  # It's good practice to import numpy for numerical operations
import plotly.express as px
#from fastaframes import to_df
#from peptacular.protein import find_peptide_indexes


# Define a function to safely calculate log2 to handle division by zero or negative numbers
def safe_log2(x):
    if x > 0:
        return math.log2(x)
    else:
        return np.nan  # Return NaN if the ratio is zero or negative

@st.cache_data
def fetch_sequence_from_uniprot(accession_number):
    url = f"https://www.uniprot.org/uniprot/{accession_number}.fasta"
    response = requests.get(url)
    if response.status_code != 200:
        st.error(f"Error fetching sequence from UniProt: {response.status_code}")
        st.stop()
        return None
    return ''.join(response.text.split('\n')[1:])  # Remove the header line


st.title("CPP Analysis Tool")

results_file = st.file_uploader("Upload results file", type=[".tsv"])
#fasta_file = st.file_uploader("Upload fasta file", type=[".fasta"])
channel_q_value_filter = st.number_input("Channel Q-value filter", value=0.10, min_value=0.0, max_value=1.0, step=0.01, help='Channel.Q.Value reflects the confidence that the precursor is indeed present in the respective channel')
remove_zeros = st.checkbox("Remove zeros", value=True, help='Remove zero values from the data')
merge_identical_ratios = st.checkbox("Merge identicle ratios", value=True, help='Merge identicle ratios')
merge_tolerance = float(st.text_input("Merge tolerance", value='0.0001', help='Merge tolerance'))
min_lysine_count = st.number_input("Min lysine count", value=1, min_value=0, max_value=10, step=1, help='Minimum number of lysines in a peptide')
max_lysine_count = st.number_input("Max lysine count", value=1, min_value=0, max_value=10, step=1, help='Maximum number of lysines in a peptide')

# Check if the file is uploaded
if results_file is not None:
    df = pd.read_csv(results_file, sep='\t')

    # Ensure the 'Stripped.Sequence' column exists to avoid KeyError
    if 'Stripped.Sequence' in df.columns:
        df = df[(df['Stripped.Sequence'].str.count('K') <= max_lysine_count) & (df['Stripped.Sequence'].str.count('K') >= min_lysine_count)]

        # Ensure the 'Channel.L' and 'Channel.H' columns exist to avoid KeyError
        if all(x in df.columns for x in ['Channel.L', 'Channel.H']):
            df['Light/Heavy.Ratio'] = df['Channel.L'] / df['Channel.H']
            # Use the safe_log2 function for the 'Light/Heavy.Log2Ratio' column
            df['Light/Heavy.Log2Ratio'] = df['Light/Heavy.Ratio'].apply(safe_log2)

            # Filter the data based on the Q-value
            df = df[df['Channel.Q.Value'] <= channel_q_value_filter]

            # Remove rows with zero values in either Channel.L or Channel.H
            if remove_zeros:
                df = df[(df['Channel.L'] != 0) & (df['Channel.H'] != 0)]

            if merge_identical_ratios:
                # In the dataframe df, some Stripped.Sequence's have Light/Heavy.Log2Ratio's that are nearly identical, but not exactly the same
                # For these instances, we must keep only the first occurrence of each row

                # Sort by 'Stripped.Sequence' and 'Light/Heavy.Log2Ratio' to ensure duplicates are ordered
                df.sort_values(by=['Stripped.Sequence', 'Light/Heavy.Log2Ratio'], inplace=True)

                # Use 'duplicated' to mark rows that have an identical sequence and a very close ratio as duplicates
                df['is_duplicate'] = df.duplicated(subset=['Stripped.Sequence'], keep='first') & \
                                     (df.groupby('Stripped.Sequence')['Light/Heavy.Log2Ratio'].diff().abs().fillna(
                                         0) < merge_tolerance)

                # Keep rows where 'is_duplicate' is False
                df = df[~df['is_duplicate']].drop(columns='is_duplicate')

            st.subheader("Data")
            st.dataframe(df)

            # Group by 'Stripped.Sequence' and calculate the required statistics
            stats_df = df.groupby(['Stripped.Sequence', 'Protein.Names'])['Light/Heavy.Log2Ratio'].agg(['mean', 'std', 'count', 'sem', 'median', 'min', 'max'])
            stats_df.reset_index(inplace=True)

            # Rename the columns for clarity
            stats_df.columns = [
                'Stripped.Sequence',
                'Protein.Names',
                'Log2Ratio.Mean',
                'Log2Ratio.Std',
                'Log2Ratio.Count',
                'Log2Ratio.SEM',
                'Log2Ratio.Median',
                'Log2Ratio.Min',
                'Log2Ratio.Max'
            ]

            st.subheader("Statistics")
            st.dataframe(stats_df)

            # Create the scatter plot with mean on the x-axis and standard deviation on the y-axis
            fig = px.scatter(stats_df, x='Log2Ratio.Mean', y='Log2Ratio.SEM', hover_name='Stripped.Sequence', hover_data=['Protein.Names', 'Log2Ratio.Count', 'Log2Ratio.SEM', 'Log2Ratio.Median', 'Log2Ratio.Min', 'Log2Ratio.Max'])
            # Enhance the plot with titles and labels
            fig.update_traces(textposition='top center')
            fig.update_layout(
                title='Volcano Plot: Mean vs. SEM',
                xaxis_title='Mean of Log2 Ratios',
                yaxis_title='SEM of Log2 Ratios',
                showlegend=False
            )

            st.plotly_chart(fig)



        else:
            st.error("The uploaded file does not contain 'Channel.L' or 'Channel.H' columns.")
    else:
        st.error("The uploaded file does not contain 'Stripped.Sequence' column.")
else:
    st.warning("Please upload a file to proceed.")
