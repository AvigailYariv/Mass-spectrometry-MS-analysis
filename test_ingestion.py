#!/usr/bin/env python
"""Test script for data ingestion module"""

from data_ingestion import load_all_samples, identify_key_columns
import logging

logging.basicConfig(level=logging.INFO)

# Load all samples
print("\n=== Loading Samples ===")
samples = load_all_samples()

print(f"\n=== Test Results ===")
print(f"Total samples: {len(samples)}")
print(f"Sample names: {list(samples.keys())}")

# Check first sample structure
if samples:
    first_sample = next(iter(samples.values()))
    print(f"\nFirst sample shape: {first_sample.shape}")
    print(f"Columns (first 15): {first_sample.columns.tolist()[:15]}")
    
    # Identify key columns
    key_cols = identify_key_columns(first_sample)
    print(f"\nKey columns identified: {key_cols}")
    
    # Show sample data
    print("\nFirst protein entry:")
    print(first_sample.iloc[0])
