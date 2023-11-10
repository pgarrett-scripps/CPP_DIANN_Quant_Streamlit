import math
import streamlit as st
import pandas as pd
import numpy as np  # It's good practice to import numpy for numerical operations
import plotly.express as px

# Define a function to safely calculate log2 to handle division by zero or negative numbers
def safe_log2(x):
    if x > 0:
        return math.log2(x)
    else:
        return np.nan  # Return NaN if the ratio is zero or negative

st.title("CPP Analysis Tool")

results_file = st.file_uploader("Upload results file", type=["tsv"])

# Check if the file is uploaded
if results_file is not None:
    df = pd.read_csv(results_file, sep='\t')

    # Ensure the 'Stripped.Sequence' column exists to avoid KeyError
    if 'Stripped.Sequence' in df.columns:
        df = df[df['Stripped.Sequence'].str.count('K') == 1]

        # Ensure the 'Channel.L' and 'Channel.H' columns exist to avoid KeyError
        if all(x in df.columns for x in ['Channel.L', 'Channel.H']):
            df['Light/Heavy.Ratio'] = df['Channel.L'] / df['Channel.H']
            # Use the safe_log2 function for the 'Light/Heavy.Log2Ratio' column
            df['Light/Heavy.Log2Ratio'] = df['Light/Heavy.Ratio'].apply(safe_log2)

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
            fig = px.scatter(stats_df, x='Log2Ratio.Mean', y='Log2Ratio.Std', hover_name='Stripped.Sequence', hover_data=['Protein.Names', 'Log2Ratio.Count', 'Log2Ratio.SEM', 'Log2Ratio.Median', 'Log2Ratio.Min', 'Log2Ratio.Max'])
            # Enhance the plot with titles and labels
            fig.update_traces(textposition='top center')
            fig.update_layout(
                title='Volcano Plot: Mean vs. Standard Deviation',
                xaxis_title='Mean of Log2 Ratios',
                yaxis_title='Standard Deviation of Log2 Ratios',
                showlegend=False
            )

            st.plotly_chart(fig)

        else:
            st.error("The uploaded file does not contain 'Channel.L' or 'Channel.H' columns.")
    else:
        st.error("The uploaded file does not contain 'Stripped.Sequence' column.")
else:
    st.warning("Please upload a file to proceed.")
