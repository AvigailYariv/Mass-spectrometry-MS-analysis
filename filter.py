"""
Filtering Module for Ribosomal Protein Selection

This module filters MS protein data to identify ribosomal and ribosome-binding proteins
based on description keywords, excluding mitochondrial ribosomes.
"""

import logging
import pandas as pd
from typing import List, Dict
from config import RIBOSOMAL_KEYWORDS, EXCLUSION_KEYWORDS

logger = logging.getLogger(__name__)


def is_ribosomal_protein(description: str) -> bool:
    """
    Determine if a protein is a ribosomal or ribosome-binding protein.
    
    A protein is considered ribosomal if its description contains:
    - Any of the RIBOSOMAL_KEYWORDS
    - AND does NOT contain any EXCLUSION_KEYWORDS (e.g., mitochondrial)
    
    Args:
        description: The protein description string to evaluate
        
    Returns:
        True if the protein is a ribosomal protein (non-mitochondrial), False otherwise
    """
    if not isinstance(description, str) or pd.isna(description):
        return False
    
    description_lower = description.lower()
    
    # Check if description contains any exclusion keywords (mitochondrial, etc.)
    for keyword in EXCLUSION_KEYWORDS:
        if keyword.lower() in description_lower:
            logger.debug(f"Excluding protein due to keyword '{keyword}': {description[:50]}")
            return False
    
    # Check if description contains any ribosomal keywords
    for keyword in RIBOSOMAL_KEYWORDS:
        if keyword.lower() in description_lower:
            logger.debug(f"Including protein due to keyword '{keyword}': {description[:50]}")
            return True
    
    return False


def filter_ribosomal_proteins(df: pd.DataFrame, description_col: str = 'Description') -> pd.DataFrame:
    """
    Filter a dataframe to retain only ribosomal proteins.
    
    Args:
        df: DataFrame containing protein data (typically from a single sample)
        description_col: Name of the column containing protein descriptions
        
    Returns:
        Filtered DataFrame containing only ribosomal proteins
    """
    if description_col not in df.columns:
        logger.error(f"Description column '{description_col}' not found in dataframe")
        logger.error(f"Available columns: {df.columns.tolist()}")
        return pd.DataFrame()
    
    # Apply filter to each row
    mask = df[description_col].apply(is_ribosomal_protein)
    filtered_df = df[mask].copy()
    
    logger.info(f"Filtered {len(df)} proteins to {len(filtered_df)} ribosomal proteins "
                f"({100*len(filtered_df)/len(df):.1f}%)")
    
    return filtered_df


def filter_all_samples(samples: Dict[str, pd.DataFrame], 
                       description_col: str = 'Description') -> Dict[str, pd.DataFrame]:
    """
    Filter ribosomal proteins from all samples.
    
    Args:
        samples: Dictionary mapping sample names to DataFrames
        description_col: Name of the column containing protein descriptions
        
    Returns:
        Dictionary mapping sample names to filtered DataFrames (ribosomal proteins only)
    """
    filtered_samples = {}
    
    for sample_name, df in samples.items():
        logger.info(f"Processing sample: {sample_name}")
        filtered_df = filter_ribosomal_proteins(df, description_col)
        filtered_samples[sample_name] = filtered_df
    
    return filtered_samples


def get_filter_statistics(original_samples: Dict[str, pd.DataFrame],
                          filtered_samples: Dict[str, pd.DataFrame]) -> Dict:
    """
    Generate statistics about the filtering process.
    
    Args:
        original_samples: Dictionary of original DataFrames
        filtered_samples: Dictionary of filtered DataFrames
        
    Returns:
        Dictionary with filtering statistics
    """
    stats = {
        'total_original_proteins': sum(len(df) for df in original_samples.values()),
        'total_ribosomal_proteins': sum(len(df) for df in filtered_samples.values()),
        'sample_stats': {}
    }
    
    for sample_name in original_samples.keys():
        orig_count = len(original_samples.get(sample_name, pd.DataFrame()))
        filt_count = len(filtered_samples.get(sample_name, pd.DataFrame()))
        
        stats['sample_stats'][sample_name] = {
            'original': orig_count,
            'filtered': filt_count,
            'percentage': (100 * filt_count / orig_count) if orig_count > 0 else 0
        }
    
    return stats
