"""
Configuration module for Mass Spectrometry Ribosomal Protein Analysis Pipeline
"""

import os
from pathlib import Path

# Project root directory
PROJECT_ROOT = Path(__file__).parent
DATA_DIR = PROJECT_ROOT / "MS data"
OUTPUT_DIR = PROJECT_ROOT / "output"
SAMPLES_DIR = DATA_DIR / "my_samples"

# Ensure output directory exists
OUTPUT_DIR.mkdir(exist_ok=True)

# Column name mappings (case-insensitive keywords to search for in Excel files)
PROTEIN_NAME_KEYWORDS = ["protein", "id", "gene"]
COVERAGE_KEYWORDS = ["coverage", "cover"]
DESCRIPTION_KEYWORDS = ["description", "name", "gene"]

# Ribosomal protein filtering keywords for DESCRIPTION field
# Include: ribosomal proteins and ribosome-binding proteins
# Exclude: mitochondrial ribosomes
RIBOSOMAL_KEYWORDS = [
    "ribosomal",
    "ribosome",
    "ribosome-binding",
    "rpl",  # Ribosomal Protein Large
    "rps",  # Ribosomal Protein Small
    "eif3",  # Eukaryotic translation initiation factor 3 (ribosome-related)
]

# Keywords to EXCLUDE (e.g., mitochondrial proteins)
EXCLUSION_KEYWORDS = [
    "mitochondrial",
    "mitochondrion",
    "mrps",  # Mitochondrial RPS
    "mrpl",  # Mitochondrial RPL
]

# Visualization settings
PLOT_DPI = 300
PLOT_FORMAT = "png"  # Can be 'png' or 'pdf'
FIGURE_SIZE = (14, 8)
SCATTER_SIZE = 100
