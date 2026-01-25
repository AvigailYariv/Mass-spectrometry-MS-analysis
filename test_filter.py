#!/usr/bin/env python
"""Test script to verify ribosomal protein filtering"""

import pandas as pd
import random
from filter import filter_ribosomal_proteins, is_ribosomal_protein, get_filter_statistics
import logging

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Load sample 1 - parse correctly
df_raw = pd.read_excel(r'MS data\sample1.xlsx', sheet_name='Proteins', header=None)
headers = df_raw.iloc[1].tolist()

# Filter for protein rows only
protein_rows = []
for i in range(2, len(df_raw)):
    row = df_raw.iloc[i]
    if pd.notna(row.iloc[0]) and row.iloc[0] in [False, True]:
        protein_rows.append(i)

data = df_raw.iloc[protein_rows].reset_index(drop=True)
data.columns = headers

print('=== Sample 1 - Filtering Results ===\n')
print(f'Total proteins before filtering: {len(data)}')

# Apply filter
filtered = filter_ribosomal_proteins(data, 'Description')

print(f'Total ribosomal proteins after filtering: {len(filtered)}')
print(f'Percentage retained: {100*len(filtered)/len(data):.1f}%')

# Show some examples of filtered proteins
print('\n=== Examples of Filtered Ribosomal Proteins ===\n')
if len(filtered) > 0:
    sample_indices = random.sample(range(len(filtered)), min(15, len(filtered)))
    
    for idx, row_idx in enumerate(sample_indices, 1):
        row = filtered.iloc[row_idx]
        gene = row['Gene Symbol']
        desc = str(row['Description'])[:120]
        coverage = row['Coverage [%]']
        print(f'{idx:2d}. {gene:15s} | Coverage: {coverage:6.0f} | {desc}')

print('\n=== Examples of Excluded Proteins ===\n')
# Show some examples of proteins that were excluded
excluded = data[~data.index.isin(filtered.index)]
if len(excluded) > 0:
    sample_indices = random.sample(range(len(excluded)), min(10, len(excluded)))
    
    for idx, row_idx in enumerate(sample_indices, 1):
        row = excluded.iloc[row_idx]
        gene = row['Gene Symbol']
        desc = str(row['Description'])[:120]
        coverage = row['Coverage [%]']
        print(f'{idx:2d}. {gene:15s} | Coverage: {coverage:6.0f} | {desc}')
