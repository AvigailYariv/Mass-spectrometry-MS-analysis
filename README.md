# Mass-spectrometry-MS-analysis
A targeted bioinformatics tool for comparative proteomics that quantifies and analyzes changes in ribosomal protein composition across different mass spectrometry samples.
---

## Overview
This program is a bioinformatics tool designed to automate the comparative analysis of Mass Spectrometry (MS) data, with a specific focus on the **ribosomal proteome**.

While standard MS tables contain thousands of proteins with various features, this tool specifically filters for ribosomal proteins and associated factors (e.g., assembly factors, chaperones). It utilizes the **Coverage [%]** metric to compare the abundance and presence of these specific proteins across multiple samples (e.g., different *E. coli* cell cultures or conditions).

## Key Features
* **Targeted Filtering:** Automatically isolates ribosome-related entries from large, complex MS datasets.
* **Multi-Sample Comparison:** Capable of ingesting multiple MS tables corresponding to different samples/cells simultaneously.
* **Coverage Analysis:** Focuses on the `Coverage [%]` feature to assess the relative detection of proteins across samples.
* **Visual Output:** Generates a comparative scatter plot for intuitive data interpretation.

## Input Data
The program expects raw Mass Spectrometry output tables (CSV, TSV, or Excel).
* **Sample Source:** Validated on *E. coli* lysates, but applicable to any organism/cell type.
* **Required Format:** The input tables must contain:
    1.  **Protein Names/IDs**
    2.  **Description/Gene Names** (for filtering purposes)
    3.  **Coverage [%]** (The specific metric used for comparison)

## Methodology
The analysis pipeline proceeds in three steps:

1.  **Data Ingestion:** The program reads multiple MS tables representing different samples.
2.  **Ribosome Filtering:** It parses the protein list to retain only relevant entries, such as:
    * Ribosomal Subunits (e.g., L, S proteins)
    * Ribosome Assembly Factors
3.  **Comparative Plotting:** It aligns the filtered data to visualize the differences in Coverage [%] between the samples.

## Visualization / Output
The primary output of the software is a comparative plot designed to highlight heterogeneity between samples.

**Plot Breakdown:**
* **X-Axis:** Ribosomal Protein Names.
* **Y-Axis:** Coverage [%] (Metric of abundance/detection).
* **Data Points:** Each dot represents a specific sample.

*This visualization allows users to quickly identify if specific ribosomal proteins have higher or lower coverage in certain cells compared to others, potentially indicating changes in ribosome stoichiometry or assembly.*

