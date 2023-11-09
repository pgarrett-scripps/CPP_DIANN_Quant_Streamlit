import math

import streamlit as st
import pandas as pd

light_file = st.file_uploader("Upload light file", type=["tsv"])
heavy_file = st.file_uploader("Upload heavy file", type=["tsv"])

nterm_modification = st.text_input('N-term Modification', value="(UniMod:330)")
light_modification = st.text_input('Light Modification', value="(UniMod:36)")
heavy_modification = st.text_input('Heavy Modification', value="(UniMod:330)")

df_light = pd.read_csv(light_file, sep='\t')
df_heavy = pd.read_csv(heavy_file, sep='\t')

# Remove peptides with more than one K
df_light = df_light[df_light['Stripped.Sequence'].str.count('K') == 1]
df_heavy = df_heavy[df_heavy['Stripped.Sequence'].str.count('K') == 1]

# remove unused columns:
df_light = df_light[['Modified.Sequence', 'Precursor.Charge', 'Ms1.Area']]
df_heavy = df_heavy[['Modified.Sequence', 'Precursor.Charge', 'Ms1.Area']]

# merge duplicate rows
df_light = df_light.groupby(['Modified.Sequence', 'Precursor.Charge']).mean().reset_index()
df_heavy = df_heavy.groupby(['Modified.Sequence', 'Precursor.Charge']).mean().reset_index()


sequence_to_area_light = {}
for index, row in df_light.iterrows():
    modified_sequence = row['Modified.Sequence'].replace(light_modification, heavy_modification)
    precursor_charge = row['Precursor.Charge']
    sequence_to_area_light[(modified_sequence, precursor_charge)] = row['Ms1.Area']

sequence_to_area_heavy = {}
for index, row in df_heavy.iterrows():
    modified_sequence = row['Modified.Sequence']
    precursor_charge = row['Precursor.Charge']
    sequence_to_area_heavy[(modified_sequence, precursor_charge)] = row['Ms1.Area']

ratio_data = {}
for (modified_sequence, precursor_charge), area_light in sequence_to_area_light.items():
    area_heavy = sequence_to_area_heavy.get((modified_sequence, precursor_charge), 0)
    if area_heavy > 0:
        ratio_data[(modified_sequence, 2)] = area_light / area_heavy

df = pd.DataFrame.from_dict(ratio_data, orient='index', columns=['Ratio'])

# Fix the index
df.index = pd.MultiIndex.from_tuples(df.index, names=['Sequence', 'Precursor.Charge'])

# Add log2 column
df['log2'] = df['Ratio'].apply(lambda x: math.log2(x))

st.dataframe(df)


