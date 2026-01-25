#!/usr/bin/env python
"""Test script for the analysis module - aggregating ribosomal proteins from all samples"""

import pandas as pd
import logging
from pathlib import Path
from config import DATA_DIR

logging.basicConfig(level=logging.INFO, format='%(name)s - %(levelname)s - %(message)s')

print("\n" + "="*70)
print("STEP 1: Loading All Samples")
print("="*70)

# Load all samples with correct parsing
samples = {}
sample_files = sorted(DATA_DIR.glob("*.xlsx"))

for excel_file in sample_files:
    sample_name = excel_file.stem
    
    try:
        # Load without header to parse correctly
        df_raw = pd.read_excel(excel_file, sheet_name='Proteins', header=None)
        
        # Get headers from row 1
        headers = df_raw.iloc[1].tolist()
        
        # Filter for protein rows only
        protein_rows = []
        for i in range(2, len(df_raw)):
            row = df_raw.iloc[i]
            if pd.notna(row.iloc[0]) and row.iloc[0] in [False, True]:
                protein_rows.append(i)
        
        # Create DataFrame with only protein rows
        data = df_raw.iloc[protein_rows].reset_index(drop=True)
        data.columns = headers
        
        samples[sample_name] = data
        print(f"  {sample_name:10s}: {len(data):5d} proteins loaded")
    except Exception as e:
        print(f"  ERROR loading {sample_name}: {str(e)}")

print(f"\nTotal samples loaded: {len(samples)}")

print("\n" + "="*70)
print("STEP 2: Filtering for Ribosomal Proteins")
print("="*70)

from filter import filter_all_samples

filtered_samples = filter_all_samples(samples, 'Description')

print("\nFiltered ribosomal proteins per sample:")
for sample_name in sorted(filtered_samples.keys()):
    df = filtered_samples[sample_name]
    print(f"  {sample_name:10s}: {len(df):3d} ribosomal proteins")

print("\n" + "="*70)
print("STEP 3: Aggregating Across All Samples")
print("="*70)

from analysis import aggregate_ribosomal_proteins, get_aggregation_summary, print_aggregation_summary, prepare_long_format

aggregated = aggregate_ribosomal_proteins(filtered_samples, 'Gene Symbol', 'Coverage [%]', 'Description')

print(f"\nAggregated DataFrame shape: {aggregated.shape}")
print(f"Columns (first 12): {aggregated.columns.tolist()[:12]}")

# Get summary statistics
summary = get_aggregation_summary(aggregated)
print_aggregation_summary(summary)

print("\n" + "="*70)
print("STEP 4: Preparing Long Format for Visualization")
print("="*70)

long_df = prepare_long_format(aggregated)

print(f"\nLong-format DataFrame shape: {long_df.shape}")
print("\nFirst 10 rows:")
print(long_df.head(10).to_string())

print("\n" + "="*70)
print("STEP 5: Sample Data Preview")
print("="*70)

print("\nTop 15 proteins by mean coverage:")
top_proteins = aggregated.nlargest(15, 'Mean_Coverage')[['Gene_Symbol', 'Mean_Coverage', 'Detected_Samples', 'Description']]
for idx, row in top_proteins.iterrows():
    desc = row['Description'][:60] if pd.notna(row['Description']) else 'N/A'
    print(f"  {row['Gene_Symbol']:15s}: {row['Mean_Coverage']:6.1f}% (in {int(row['Detected_Samples']):2d}/13 samples) - {desc}")

print("\n" + "="*70)
