# Mass Spectrometry Ribosomal Protein Analysis Pipeline - Project Plan

## Project Overview

Build a Python tool to read 13 MS Excel sample files, filter for ribosomal proteins across all samples, and analyze coverage distribution. The pipeline ingests multi-sheet Excel files, applies protein name filtering, aggregates coverage data, and produces comparative scatter plot visualizations showing protein abundance heterogeneity between samples.

## Implementation Steps

### 1. Set up project structure and dependencies
Create `main.py`, `config.py`, and a `requirements.txt` with pandas, openpyxl, matplotlib, and seaborn for data handling and visualization.

### 2. Build Excel data ingestion module
Write a function to read each `.xlsx` file, identify sheets containing protein data, extract columns for protein names, descriptions, and coverage percentages, and handle variable column naming across samples.

### 3. Implement ribosomal protein filtering logic
Create a filter function using keyword matching (e.g., "ribosomal", "L-protein", "S-protein", "rpl", "rps") on protein names/descriptions to isolate ribosome-related entries.
Consider multiple naming conventions, including:
    * `"ribosomal"`
    * `"ribosome-binding"`

### 4. Aggregate multi-sample coverage data
Combine filtered results from all 13 samples into a unified data structure mapping each protein to coverage values across samples, handling missing proteins gracefully.

### 5. Generate comparative scatter plot visualization
Create plots with protein names on X-axis and coverage [%] on Y-axis, with each sample represented by distinct colors/markers, saved as high-quality PNG/PDF output.

### 6. Create main execution script
Build a driver script that orchestrates the pipeline: loads all samples → filters → aggregates → plots, with progress logging and error handling.

## Context & Data

- **Total Samples:** 13 Excel files (`sample1.xlsx` through `sample13.xlsx`)
- **Location:** `MS data/` directory
- **File Format:** `.xlsx` (Excel workbook)
- **Key Columns Expected:** Protein Names/IDs, Description/Gene Names, Coverage [%]
- **Target Organism:** E. coli (validated use case, extensible to others)

## Key Design Considerations

1. **Excel sheet structure assumption** 
   - Do you know which sheet(s) contain the protein data in each Excel file? 
   - Some samples might have multiple sheets or use different sheet names
   - Should the tool auto-detect the relevant sheet or use a configuration file?

2. **Ribosomal protein definition**
   - Which protein naming conventions does the dataset use (UniProt IDs, gene symbols like "RPL23", "RPS19", or full descriptions)?
   - Should we validate against a known ribosomal protein database?

3. **Output and reporting**
   - Besides scatter plots, would you want a summary CSV/table listing all filtered proteins and their coverage across samples?
   - Should we include statistical comparisons between samples?

## Expected Deliverables

- `main.py` — Main execution script orchestrating the pipeline
- `config.py` — Configuration file with paths, filtering keywords, plot settings
- `data_ingestion.py` — Module for reading Excel files
- `filter.py` — Module for protein filtering logic
- `analysis.py` — Module for data aggregation and analysis
- `visualization.py` — Module for scatter plot generation
- `requirements.txt` — Python package dependencies
- Output plots (PNG/PDF) comparing protein coverage across samples
- Optional: Summary CSV with filtered protein results
