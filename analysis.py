"""
Analysis Module for Multi-Sample Ribosomal Protein Coverage Aggregation

This module aggregates filtered ribosomal protein data from multiple samples
into a unified data structure for comparative analysis.
"""

import logging
import pandas as pd
import numpy as np
from typing import Dict, Tuple, List

logger = logging.getLogger(__name__)


def aggregate_ribosomal_proteins(filtered_samples: Dict[str, pd.DataFrame],
                                  gene_col: str = 'Gene Symbol',
                                  coverage_col: str = 'Coverage [%]',
                                  description_col: str = 'Description') -> pd.DataFrame:
    """
    Aggregate ribosomal protein coverage data from all samples into a single DataFrame.
    
    Creates a unified view where:
    - Each row represents a unique ribosomal protein (by Gene Symbol)
    - Each column represents coverage values from a sample
    - Missing values (proteins not detected in a sample) are NaN
    
    Args:
        filtered_samples: Dictionary mapping sample names to filtered DataFrames
        gene_col: Name of the gene symbol column
        coverage_col: Name of the coverage percentage column
        description_col: Name of the description column
        
    Returns:
        DataFrame with proteins as rows and samples as columns, plus metadata columns
    """
    
    if not filtered_samples:
        logger.error("No filtered samples provided")
        return pd.DataFrame()
    
    # Dictionary to store aggregated data
    protein_data = {}  # {gene_symbol: {sample_name: coverage_value, ...}}
    protein_metadata = {}  # {gene_symbol: description}
    
    logger.info(f"Aggregating data from {len(filtered_samples)} samples")
    
    # Iterate through each sample
    for sample_name, df in filtered_samples.items():
        logger.debug(f"Processing sample '{sample_name}' with {len(df)} ribosomal proteins")
        
        for idx, row in df.iterrows():
            gene_symbol = str(row[gene_col]) if pd.notna(row[gene_col]) else f"Unknown_{idx}"
            coverage = row[coverage_col] if pd.notna(row[coverage_col]) else np.nan
            description = str(row[description_col]) if pd.notna(row[description_col]) else "N/A"
            
            # Initialize protein entry if not exists
            if gene_symbol not in protein_data:
                protein_data[gene_symbol] = {}
                protein_metadata[gene_symbol] = description
            
            # Store coverage value for this sample
            protein_data[gene_symbol][sample_name] = coverage
    
    # Convert to DataFrame
    aggregated_df = pd.DataFrame.from_dict(protein_data, orient='index')
    
    # Add metadata columns
    aggregated_df.insert(0, 'Description', aggregated_df.index.map(protein_metadata))
    aggregated_df.insert(0, 'Gene_Symbol', aggregated_df.index)
    
    # Calculate statistics
    sample_names = sorted([s for s in aggregated_df.columns if s not in ['Gene_Symbol', 'Description']])
    
    # Coverage statistics across samples
    aggregated_df['Mean_Coverage'] = aggregated_df[sample_names].mean(axis=1)
    aggregated_df['Std_Coverage'] = aggregated_df[sample_names].std(axis=1)
    aggregated_df['Min_Coverage'] = aggregated_df[sample_names].min(axis=1)
    aggregated_df['Max_Coverage'] = aggregated_df[sample_names].max(axis=1)
    aggregated_df['Detected_Samples'] = aggregated_df[sample_names].notna().sum(axis=1)
    aggregated_df['Detection_Rate'] = 100 * aggregated_df['Detected_Samples'] / len(sample_names)
    
    logger.info(f"Aggregated {len(aggregated_df)} unique ribosomal proteins")
    logger.info(f"Sample columns: {sample_names}")
    
    return aggregated_df


def get_aggregation_summary(aggregated_df: pd.DataFrame) -> Dict:
    """
    Generate summary statistics about the aggregated data.
    
    Args:
        aggregated_df: The aggregated protein coverage DataFrame
        
    Returns:
        Dictionary with summary statistics
    """
    
    # Find sample columns (exclude metadata columns)
    metadata_cols = {'Gene_Symbol', 'Description', 'Mean_Coverage', 'Std_Coverage',
                     'Min_Coverage', 'Max_Coverage', 'Detected_Samples', 'Detection_Rate'}
    sample_cols = [col for col in aggregated_df.columns if col not in metadata_cols]
    
    summary = {
        'total_unique_proteins': len(aggregated_df),
        'num_samples': len(sample_cols),
        'sample_names': sample_cols,
        'proteins_detected_all_samples': len(aggregated_df[aggregated_df['Detection_Rate'] == 100]),
        'mean_detection_rate': aggregated_df['Detection_Rate'].mean(),
        'mean_coverage': aggregated_df['Mean_Coverage'].mean(),
        'median_coverage': aggregated_df['Mean_Coverage'].median(),
        'min_coverage': aggregated_df['Mean_Coverage'].min(),
        'max_coverage': aggregated_df['Mean_Coverage'].max(),
    }
    
    # Per-sample statistics
    summary['per_sample_stats'] = {}
    for sample in sample_cols:
        non_null = aggregated_df[sample].notna().sum()
        summary['per_sample_stats'][sample] = {
            'proteins_detected': non_null,
            'mean_coverage': aggregated_df[sample].mean(),
            'median_coverage': aggregated_df[sample].median(),
            'min_coverage': aggregated_df[sample].min(),
            'max_coverage': aggregated_df[sample].max(),
        }
    
    return summary


def prepare_long_format(aggregated_df: pd.DataFrame) -> pd.DataFrame:
    """
    Convert the wide-format aggregated DataFrame to long format for visualization.
    
    Long format has columns: Gene_Symbol, Sample, Coverage
    This format is ideal for scatter plots and other visualizations.
    
    Args:
        aggregated_df: The wide-format aggregated DataFrame
        
    Returns:
        Long-format DataFrame suitable for plotting
    """
    
    # Find sample columns (exclude metadata)
    metadata_cols = {'Gene_Symbol', 'Description', 'Mean_Coverage', 'Std_Coverage',
                     'Min_Coverage', 'Max_Coverage', 'Detected_Samples', 'Detection_Rate'}
    sample_cols = [col for col in aggregated_df.columns if col not in metadata_cols]
    
    # Melt the DataFrame
    long_df = aggregated_df[['Gene_Symbol', 'Description'] + sample_cols].melt(
        id_vars=['Gene_Symbol', 'Description'],
        var_name='Sample',
        value_name='Coverage'
    )
    
    # Remove NaN values (proteins not detected in a sample)
    long_df = long_df.dropna(subset=['Coverage'])
    
    logger.info(f"Created long-format DataFrame with {len(long_df)} coverage measurements")
    
    return long_df


def print_aggregation_summary(summary: Dict):
    """
    Print a formatted summary of the aggregation.
    
    Args:
        summary: Dictionary from get_aggregation_summary()
    """
    print("\n" + "="*70)
    print("RIBOSOMAL PROTEIN AGGREGATION SUMMARY")
    print("="*70)
    print(f"\nTotal unique ribosomal proteins detected: {summary['total_unique_proteins']}")
    print(f"Number of samples analyzed: {summary['num_samples']}")
    print(f"Proteins detected in all samples: {summary['proteins_detected_all_samples']}")
    print(f"Average detection rate across proteins: {summary['mean_detection_rate']:.1f}%")
    print(f"\nCoverage Statistics (across all proteins):")
    print(f"  Mean coverage: {summary['mean_coverage']:.1f}%")
    print(f"  Median coverage: {summary['median_coverage']:.1f}%")
    print(f"  Min coverage: {summary['min_coverage']:.1f}%")
    print(f"  Max coverage: {summary['max_coverage']:.1f}%")
    
    print(f"\nPer-Sample Statistics:")
    print("-" * 70)
    for sample in summary['sample_names']:
        stats = summary['per_sample_stats'][sample]
        print(f"{sample:15s}: {stats['proteins_detected']:3d} proteins | "
              f"Cov: {stats['mean_coverage']:5.1f}% ± {stats['max_coverage']-stats['min_coverage']:5.1f}%")
    print("="*70 + "\n")
