#!/usr/bin/env python
"""Script to display random protein descriptions from sample 1 for verification"""

import pandas as pd
import random

# Load without header to understand structure
df_raw = pd.read_excel(r'MS data\sample1.xlsx', sheet_name='Proteins', header=None)

# Row 1 has headers
headers = df_raw.iloc[1].tolist()

# Filter for protein rows (first column is False/True for Checked, second column is "Specific" for proteins)
# Skip rows that start with NaN or have "Confidence" in column 1 (these are peptide sub-rows)
protein_rows = []
for i in range(2, len(df_raw)):
    row = df_raw.iloc[i]
    # Check if this looks like a protein row (not empty first column, contains Accession data)
    if pd.notna(row.iloc[0]) and row.iloc[0] in [False, True]:
        protein_rows.append(i)

# Create dataframe with only protein rows
data = df_raw.iloc[protein_rows].reset_index(drop=True)
data.columns = headers

print('=== Sample 1 - Protein Data Overview ===\n')
print(f'Total proteins in sample: {len(data)}')
print(f'DataFrame shape: {data.shape}')

# Get 30 random indices
random_indices = random.sample(range(len(data)), min(30, len(data)))

print(f'\n=== 30 Random Protein Descriptions from Sample 1 ===\n')

for i, idx in enumerate(random_indices, 1):
    row = data.iloc[idx]
    gene_symbol = str(row['Gene Symbol']) if pd.notna(row['Gene Symbol']) else 'N/A'
    description = str(row['Description']) if pd.notna(row['Description']) else 'N/A'
    coverage = str(row['Coverage [%]']) if pd.notna(row['Coverage [%]']) else 'N/A'
    
    # Truncate long descriptions
    if len(description) > 130:
        description = description[:130] + "..."
    
    print(f'{i:2d}. Gene: {gene_symbol:15s} | Coverage: {coverage:7s} | Description:')
    print(f'    {description}\n')
