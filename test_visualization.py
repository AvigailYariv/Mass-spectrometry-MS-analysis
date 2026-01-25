#!/usr/bin/env python
"""Test script for the visualization module"""

import pandas as pd
import logging
from pathlib import Path
from config import DATA_DIR, OUTPUT_DIR

logging.basicConfig(level=logging.INFO, format='%(name)s - %(levelname)s - %(message)s')

print("\n" + "="*70)
print("VISUALIZATION TEST - Loading and Aggregating Data")
print("="*70 + "\n")

# Load and parse all samples
samples = {}
sample_files = sorted(DATA_DIR.glob("*.xlsx"))

for excel_file in sample_files[:3]:  # Test with first 3 samples for speed
    sample_name = excel_file.stem
    
    try:
        df_raw = pd.read_excel(excel_file, sheet_name='Proteins', header=None)
        headers = df_raw.iloc[1].tolist()
        
        protein_rows = []
        for i in range(2, len(df_raw)):
            row = df_raw.iloc[i]
            if pd.notna(row.iloc[0]) and row.iloc[0] in [False, True]:
                protein_rows.append(i)
        
        data = df_raw.iloc[protein_rows].reset_index(drop=True)
        data.columns = headers
        
        samples[sample_name] = data
        print(f"Loaded {sample_name}: {len(data)} proteins")
    except Exception as e:
        print(f"Error loading {sample_name}: {str(e)}")

# Filter for ribosomal proteins
from filter import filter_all_samples

filtered_samples = filter_all_samples(samples, 'Description')

print("\nFiltered ribosomal proteins:")
for sample_name, df in filtered_samples.items():
    print(f"  {sample_name}: {len(df)} ribosomal proteins")

# Aggregate
from analysis import aggregate_ribosomal_proteins, prepare_long_format

aggregated = aggregate_ribosomal_proteins(filtered_samples, 'Gene Symbol', 'Coverage [%]', 'Description')
long_df = prepare_long_format(aggregated)

print(f"\nAggregated: {aggregated.shape[0]} unique proteins")
print(f"Long format: {long_df.shape[0]} measurements")

# Create visualizations
print("\n" + "="*70)
print("Creating Visualizations")
print("="*70 + "\n")

from visualization import (create_scatter_plot, create_boxplot, 
                          create_heatmap, create_faceted_scatter_plot)

# Scatter plot
print("Creating scatter plot...")
scatter_path = create_scatter_plot(long_df, output_file="scatter_coverage_all")
if scatter_path:
    print(f"✓ Saved: {scatter_path}")

# Boxplot
print("Creating boxplot...")
boxplot_path = create_boxplot(long_df, output_file="boxplot_coverage")
if boxplot_path:
    print(f"✓ Saved: {boxplot_path}")

# Heatmap
print("Creating heatmap...")
heatmap_path = create_heatmap(aggregated, output_file="heatmap_coverage")
if heatmap_path:
    print(f"✓ Saved: {heatmap_path}")

print("\n" + "="*70)
print("Visualization Test Complete!")
print(f"Output directory: {OUTPUT_DIR}")
print("="*70 + "\n")
