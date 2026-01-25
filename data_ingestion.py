"""
Data Ingestion Module for Mass Spectrometry Analysis

This module handles reading and parsing MS data from Excel files.
It identifies relevant sheets, extracts protein data, and normalizes column names.
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple
import pandas as pd
from config import DATA_DIR, PROTEIN_NAME_KEYWORDS, COVERAGE_KEYWORDS, DESCRIPTION_KEYWORDS

logger = logging.getLogger(__name__)


def find_sample_files(data_dir: Path = DATA_DIR) -> List[Path]:
    """
    Find all Excel sample files in the data directory.
    
    Args:
        data_dir: Path to the data directory
        
    Returns:
        List of sorted Path objects for all .xlsx files
    """
    if not data_dir.exists():
        logger.error(f"Data directory not found: {data_dir}")
        return []
    
    sample_files = sorted(data_dir.glob("*.xlsx"))
    logger.info(f"Found {len(sample_files)} Excel files in {data_dir}")
    
    return sample_files


def get_sheet_names(excel_file: Path) -> List[str]:
    """
    Get all sheet names from an Excel file.
    
    Args:
        excel_file: Path to the Excel file
        
    Returns:
        List of sheet names
    """
    try:
        xl_file = pd.ExcelFile(excel_file)
        return xl_file.sheet_names
    except Exception as e:
        logger.error(f"Error reading sheet names from {excel_file.name}: {str(e)}")
        return []


def find_data_sheet(excel_file: Path) -> str:
    """
    Find the sheet containing protein data.
    Prioritizes sheets named "Proteins" or "Sheet1".
    
    Args:
        excel_file: Path to the Excel file
        
    Returns:
        Sheet name to use for data extraction, or None if not found
    """
    sheet_names = get_sheet_names(excel_file)
    
    # Check for common protein data sheet names
    priority_sheets = ["Proteins", "Sheet1", "Data", "Results"]
    
    for sheet in priority_sheets:
        if sheet in sheet_names:
            logger.debug(f"Using sheet '{sheet}' from {excel_file.name}")
            return sheet
    
    # If no priority sheet found, use the first one
    if sheet_names:
        logger.debug(f"Using first sheet '{sheet_names[0]}' from {excel_file.name}")
        return sheet_names[0]
    
    logger.warning(f"No sheets found in {excel_file.name}")
    return None


def identify_key_columns(df: pd.DataFrame) -> Dict[str, str]:
    """
    Identify key columns in the dataframe by searching for keywords.
    Searches both column names and first data row.
    
    Args:
        df: DataFrame to search
        
    Returns:
        Dictionary mapping data types to column names:
        {
            'protein_id': 'column_name',
            'gene_symbol': 'column_name',
            'description': 'column_name',
            'coverage': 'column_name'
        }
    """
    columns_found = {}
    
    # Search through columns for matches
    column_list = df.columns.tolist()
    
    for col in column_list:
        col_lower = str(col).lower()
        
        # Check for protein ID column
        if 'accession' in col_lower or 'protein' in col_lower:
            if 'protein_id' not in columns_found:
                columns_found['protein_id'] = col
        
        # Check for gene symbol column
        if 'gene' in col_lower and 'symbol' in col_lower:
            if 'gene_symbol' not in columns_found:
                columns_found['gene_symbol'] = col
        
        # Check for description column
        if 'description' in col_lower:
            if 'description' not in columns_found:
                columns_found['description'] = col
        
        # Check for coverage column
        if 'coverage' in col_lower:
            if 'coverage' not in columns_found:
                columns_found['coverage'] = col
    
    logger.debug(f"Identified columns: {columns_found}")
    return columns_found


def read_sample_file(excel_file: Path) -> Tuple[pd.DataFrame, str]:
    """
    Read and parse a single Excel sample file.
    
    Args:
        excel_file: Path to the Excel file
        
    Returns:
        Tuple of (DataFrame with protein data, sample name)
    """
    try:
        # Find the correct sheet
        sheet_name = find_data_sheet(excel_file)
        if not sheet_name:
            logger.warning(f"Could not find data sheet in {excel_file.name}")
            return None, excel_file.stem
        
        # Read the Excel file
        # The header is in the first data row, so we skip it and set it as header
        df = pd.read_excel(excel_file, sheet_name=sheet_name, header=0)
        
        logger.info(f"Loaded {len(df)} proteins from {excel_file.name}")
        
        # Extract sample name from filename
        sample_name = excel_file.stem
        
        return df, sample_name
        
    except Exception as e:
        logger.error(f"Error reading {excel_file.name}: {str(e)}")
        return None, excel_file.stem


def load_all_samples() -> Dict[str, pd.DataFrame]:
    """
    Load all sample files from the data directory.
    
    Returns:
        Dictionary mapping sample names to DataFrames
    """
    sample_files = find_sample_files()
    samples = {}
    
    for excel_file in sample_files:
        df, sample_name = read_sample_file(excel_file)
        if df is not None:
            samples[sample_name] = df
        else:
            logger.warning(f"Skipped {sample_name} due to read errors")
    
    logger.info(f"Successfully loaded {len(samples)} samples")
    return samples
