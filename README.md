# Mass Spectrometry Ribosomal Protein Analysis Pipeline
A comprehensive Python-based bioinformatics tool for comparative proteomics analysis that automates the discovery, filtering, and visualization of ribosomal protein coverage across multiple mass spectrometry samples.

---

## Overview
This tool is designed to automate the comparative analysis of Mass Spectrometry (MS) data with a specific focus on the **ribosomal proteome**. 

While standard MS tables contain thousands of proteins, this tool automatically filters for ribosomal proteins and ribosome-binding factors, aggregates data across multiple samples, and generates comprehensive visualizations. It utilizes the **Coverage [%]** metric to compare the abundance and presence of these specific proteins across samples (e.g., different *E. coli* cell cultures or experimental conditions).

## Key Features
* **Automated Data Ingestion:** Reads and parses multiple Excel MS files (`.xlsx`) with complex multi-row headers
* **Targeted Filtering:** Intelligently isolates ribosomal proteins (40S/60S subunits) and ribosome-binding factors from large MS datasets
* **Mitochondrial Exclusion:** Automatically excludes mitochondrial ribosomes (MRPS, MRPL) to focus on cytoplasmic ribosomes
* **Multi-Sample Aggregation:** Combines data from 13+ samples into a unified analytical framework
* **Comprehensive Statistics:** Calculates mean, median, min/max coverage, detection rates, and per-sample metrics
* **Multiple Visualizations:** Generates scatter plots, boxplots, and heatmaps for intuitive data interpretation
* **Production-Ready:** Includes logging, error handling, and CSV export of results

## Input Data
The program expects mass spectrometry output tables in **Excel format** (`.xlsx`).

**Sample Source:** 
- Validated on *E. coli* lysates
- Applicable to any organism/cell type

**Required Format:** 
Each Excel file must contain a "Proteins" sheet with:
- **Accession:** Protein ID/accession number
- **Gene Symbol:** Gene name (e.g., RPS19, RPL12)
- **Description:** Full protein description (used for ribosomal protein filtering)
- **Coverage [%]:** Coverage percentage (the metric used for comparison)
- Additional columns are acceptable and will be preserved in aggregated output

**Example Data Structure:**
```
Input: 13 Excel files (sample1.xlsx through sample13.xlsx)
  ├── sample1.xlsx: 899 proteins → 77 ribosomal proteins
  ├── sample2.xlsx: 1,123 proteins → 90 ribosomal proteins
  └── ... (total 162 unique ribosomal proteins across all samples)
```

## Methodology & Pipeline

The analysis pipeline proceeds through **6 main steps**:

### Step 1: Data Ingestion
- Locates all Excel files in the `MS data/` directory
- Parses complex multi-row headers and identifies the protein data sheet
- Handles variable column naming across different MS processing pipelines
- Loads protein data with accession, gene symbol, description, and coverage fields

### Step 2: Ribosomal Protein Filtering
- Scans protein descriptions for ribosomal keywords:
  - `"ribosomal"`, `"ribosome"`, `"ribosome-binding"`, `"rpl"`, `"rps"`, `"eif3"`
- **Excludes** mitochondrial ribosomes (keywords: `"mitochondrial"`, `"mrps"`, `"mrpl"`)
- Retains only cytoplasmic ribosomal proteins and ribosome-binding factors

### Step 3: Data Aggregation
- Combines filtered results from all 13 samples into a unified dataset
- Creates protein-to-sample mapping with coverage values
- Handles missing proteins gracefully (NaN for proteins not detected in specific samples)
- Calculates aggregate statistics: mean coverage, detection rates, min/max ranges

### Step 4: Comparative Analysis
- Generates comprehensive summary statistics
- Computes per-protein and per-sample metrics
- Identifies proteins detected across all samples vs. sample-specific proteins

### Step 5: Visualization
- **Scatter Plot:** All proteins with coverage values colored by sample (shows heterogeneity)
- **Faceted Scatter Plots:** Splits proteins into manageable groups (25 proteins per plot)
- **Boxplot:** Coverage distribution for each protein across all samples
- **Heatmap:** Full protein × sample coverage matrix

### Step 6: Results Export
- Saves aggregated data as CSV for downstream analysis
- Generates summary report with statistics
- Stores all visualizations as high-resolution PNG files (300 DPI)

## Visualization & Output

### Generated Outputs

**Visualization Files (PNG, 300 DPI):**
1. `ribosomal_coverage_scatter_all.png` - All proteins scatter plot
2. `ribosomal_coverage_scatter_part[1-N].png` - Faceted scatter plots (25 proteins each)
3. `ribosomal_coverage_boxplot.png` - Coverage distribution boxplot
4. `ribosomal_coverage_heatmap.png` - Full coverage heatmap

**Data Files:**
- `ribosomal_proteins_aggregated.csv` - Complete aggregated dataset with all proteins and samples
- `summary_report.txt` - Comprehensive statistics and per-sample analysis

### Plot Descriptions

**Scatter Plot:**
- **X-Axis:** Ribosomal Protein Gene Symbols (sorted alphabetically)
- **Y-Axis:** Coverage [%] (0-100)
- **Points:** Each dot represents a specific sample (different colors)
- **Purpose:** Visualize coverage heterogeneity between samples for each protein

**Boxplot:**
- **X-Axis:** Ribosomal Protein Gene Symbols
- **Y-Axis:** Coverage [%]
- **Boxes:** Show coverage distribution (quartiles) for each protein across all samples
- **Purpose:** Identify proteins with high/low/variable coverage

**Heatmap:**
- **Rows:** Individual ribosomal proteins
- **Columns:** Samples
- **Colors:** Intensity represents coverage [%] (yellow to red gradient)
- **Purpose:** Quick overview of detection patterns across the entire dataset

*These visualizations enable users to quickly identify ribosomal proteins with differential coverage across samples, potentially indicating changes in ribosome composition or assembly.*

---

## Installation & Setup

### Prerequisites
- **Python 3.8+**
- **pip** (Python package manager)

### Installation Steps

1. **Clone or Download the Repository:**
   ```bash
   cd Mass-spectrometry-MS-analysis
   ```

2. **Install Dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

   This installs:
   - `pandas` - Data manipulation and analysis
   - `openpyxl` - Excel file reading
   - `matplotlib` - Plotting library
   - `seaborn` - Statistical data visualization
   - `numpy` - Numerical computing

3. **Verify Installation:**
   ```bash
   python -c "import pandas, openpyxl, matplotlib, seaborn, numpy; print('All dependencies installed successfully!')"
   ```

---

## Usage

### Quick Start

1. **Place your MS data files** in the `MS data/` directory (Excel files named `sample1.xlsx`, `sample2.xlsx`, etc.)

2. **Run the pipeline:**
   ```bash
   python main.py
   ```

3. **Check the output:**
   Results will be saved to the `output/` directory including visualizations and CSV data

### Input Directory Structure
```
Mass-spectrometry-MS-analysis/
└── MS data/
    ├── sample1.xlsx
    ├── sample2.xlsx
    ├── sample3.xlsx
    └── ... (up to 13+ samples)
```

### Output Files
```
output/
├── ribosomal_coverage_scatter_all.png          # Main scatter plot
├── ribosomal_coverage_scatter_part[1-N].png   # Faceted plots (if >30 proteins)
├── ribosomal_coverage_boxplot.png              # Box plot
├── ribosomal_coverage_heatmap.png              # Heatmap
├── ribosomal_proteins_aggregated.csv           # Aggregated data (all proteins × samples)
└── summary_report.txt                          # Statistical summary
```

---

## Project Structure

```
Mass-spectrometry-MS-analysis/
├── main.py                              # Main orchestration script
├── config.py                            # Configuration, paths, and keywords
├── data_ingestion.py                    # Excel file parsing module
├── filter.py                            # Ribosomal protein filtering logic
├── analysis.py                          # Data aggregation and statistics
├── visualization.py                     # Plot generation module
├── requirements.txt                     # Python dependencies
├── README.md                            # This file
├── MS data/                             # Input directory
│   ├── sample1.xlsx
│   ├── sample2.xlsx
│   └── ... (your MS data files)
└── output/                              # Generated results
    ├── *.png                            # Visualization plots
    ├── *.csv                            # Aggregated data
    └── *.txt                            # Summary reports
```

---

## Module Documentation

### `config.py`
Configuration file containing:
- Directory paths (input, output)
- Ribosomal protein filtering keywords
- Exclusion keywords (mitochondrial proteins)
- Visualization settings (DPI, figure size, colors)

### `data_ingestion.py`
Handles:
- Finding and locating Excel sample files
- Parsing complex multi-row headers
- Identifying correct data sheets
- Extracting protein information

### `filter.py`
Implements ribosomal protein filtering:
- Keyword matching in protein descriptions
- Mitochondrial ribosome exclusion
- Per-sample and multi-sample filtering

### `analysis.py`
Aggregates and analyzes filtered data:
- Combines data from multiple samples
- Calculates coverage statistics
- Prepares data for visualization

### `visualization.py`
Generates publication-quality plots:
- Scatter plots with sample differentiation
- Boxplots for distribution analysis
- Heatmaps for overview visualization
- Faceted plots for large datasets

### `main.py`
Orchestrates the complete pipeline:
- Coordinates all analysis steps
- Implements error handling
- Logs pipeline progress
- Exports results

---

## Customization

### Modifying Filter Keywords

Edit `config.py` to change which proteins are considered "ribosomal":

```python
RIBOSOMAL_KEYWORDS = [
    "ribosomal",
    "ribosome",
    # Add more keywords as needed
]

EXCLUSION_KEYWORDS = [
    "mitochondrial",
    # Add proteins to exclude
]
```

### Changing Output Format

Modify visualization settings in `config.py`:

```python
PLOT_FORMAT = "pdf"    # Change from "png" to "pdf"
PLOT_DPI = 600         # Increase resolution
FIGURE_SIZE = (16, 10) # Larger plots
```

### Adjusting Plot Parameters

Edit visualization settings in `visualization.py` or pass custom parameters when calling plot functions.

---

## Troubleshooting

### Issue: "ModuleNotFoundError: No module named 'seaborn'"
**Solution:** Install missing packages:
```bash
pip install seaborn
```

---

