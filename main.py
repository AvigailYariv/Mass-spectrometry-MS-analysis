"""
Main execution script for Mass Spectrometry Ribosomal Protein Analysis Pipeline

This script orchestrates the complete pipeline:
1. Load all MS Excel sample files (13 samples)
2. Filter for ribosomal proteins
3. Aggregate coverage data across samples
4. Generate comparative visualizations
5. Generate summary statistics and reports
"""

import logging
import pandas as pd
from pathlib import Path
from config import PROJECT_ROOT, DATA_DIR, OUTPUT_DIR
from data_ingestion import find_sample_files
from filter import filter_all_samples, get_filter_statistics
from analysis import (aggregate_ribosomal_proteins, get_aggregation_summary, 
                     print_aggregation_summary, prepare_long_format)
from visualization import (create_scatter_plot, create_boxplot, 
                          create_heatmap, create_faceted_scatter_plot)

# Configure logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)


def load_all_samples_raw() -> dict:
    """
    Load all sample files with proper parsing.
    
    Returns:
        Dictionary mapping sample names to DataFrames
    """
    logger.info("="*70)
    logger.info("STEP 1: Loading All Sample Files")
    logger.info("="*70)
    
    samples = {}
    sample_files = find_sample_files()
    
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
            logger.info(f"  {sample_name:12s}: {len(data):5d} proteins loaded")
        except Exception as e:
            logger.error(f"Error loading {sample_name}: {str(e)}")
    
    logger.info(f"\nTotal samples loaded: {len(samples)}")
    return samples


def filter_ribosomal_proteins(samples: dict) -> dict:
    """
    Filter all samples for ribosomal proteins.
    
    Returns:
        Dictionary mapping sample names to filtered DataFrames
    """
    logger.info("\n" + "="*70)
    logger.info("STEP 2: Filtering for Ribosomal Proteins")
    logger.info("="*70 + "\n")
    
    filtered_samples = filter_all_samples(samples, 'Description')
    
    logger.info("Filtered ribosomal proteins per sample:")
    for sample_name in sorted(filtered_samples.keys()):
        df = filtered_samples[sample_name]
        orig_count = len(samples.get(sample_name, pd.DataFrame()))
        pct = 100 * len(df) / orig_count if orig_count > 0 else 0
        logger.info(f"  {sample_name:12s}: {len(df):3d} / {orig_count:5d} proteins ({pct:5.1f}%)")
    
    # Get filter statistics
    stats = get_filter_statistics(samples, filtered_samples)
    logger.info(f"\nTotal ribosomal proteins across all samples: {stats['total_ribosomal_proteins']}")
    logger.info(f"Total original proteins: {stats['total_original_proteins']}")
    
    return filtered_samples


def aggregate_data(filtered_samples: dict) -> tuple:
    """
    Aggregate filtered data across all samples.
    
    Returns:
        Tuple of (aggregated_df, long_df, summary_dict)
    """
    logger.info("\n" + "="*70)
    logger.info("STEP 3: Aggregating Data Across All Samples")
    logger.info("="*70 + "\n")
    
    aggregated = aggregate_ribosomal_proteins(filtered_samples, 'Gene Symbol', 'Coverage [%]', 'Description')
    
    logger.info(f"Aggregated DataFrame shape: {aggregated.shape}")
    logger.info(f"Unique proteins: {len(aggregated)}")
    
    # Prepare long format and get summary
    long_df = prepare_long_format(aggregated)
    summary = get_aggregation_summary(aggregated)
    
    return aggregated, long_df, summary


def generate_visualizations(aggregated_df: pd.DataFrame, long_df: pd.DataFrame) -> list:
    """
    Generate all visualization plots.
    
    Returns:
        List of paths to generated plots
    """
    logger.info("\n" + "="*70)
    logger.info("STEP 4: Generating Visualizations")
    logger.info("="*70 + "\n")
    
    output_files = []
    
    # Scatter plot
    logger.info("Creating scatter plot (all proteins)...")
    scatter_path = create_scatter_plot(long_df, output_file="ribosomal_coverage_scatter_all")
    if scatter_path:
        output_files.append(scatter_path)
        logger.info(f"  ✓ Saved: {scatter_path.name}")
    
    # Faceted scatter plots (if many proteins)
    num_proteins = long_df['Gene_Symbol'].nunique()
    if num_proteins > 30:
        logger.info(f"Creating faceted scatter plots ({num_proteins} proteins)...")
        faceted_paths = create_faceted_scatter_plot(long_df, proteins_per_plot=25, 
                                                    output_prefix="ribosomal_coverage_scatter")
        output_files.extend(faceted_paths)
        logger.info(f"  ✓ Saved {len(faceted_paths)} plot parts")
    
    # Boxplot
    logger.info("Creating boxplot...")
    boxplot_path = create_boxplot(long_df, output_file="ribosomal_coverage_boxplot")
    if boxplot_path:
        output_files.append(boxplot_path)
        logger.info(f"  ✓ Saved: {boxplot_path.name}")
    
    # Heatmap
    logger.info("Creating heatmap...")
    heatmap_path = create_heatmap(aggregated_df, output_file="ribosomal_coverage_heatmap")
    if heatmap_path:
        output_files.append(heatmap_path)
        logger.info(f"  ✓ Saved: {heatmap_path.name}")
    
    return output_files


def save_results(aggregated_df: pd.DataFrame, summary: dict):
    """
    Save aggregated results to CSV and generate summary report.
    """
    logger.info("\n" + "="*70)
    logger.info("STEP 5: Saving Results and Reports")
    logger.info("="*70 + "\n")
    
    # Save aggregated data to CSV
    csv_path = OUTPUT_DIR / "ribosomal_proteins_aggregated.csv"
    aggregated_df.to_csv(csv_path, index=False)
    logger.info(f"Saved aggregated data to: {csv_path.name}")
    
    # Generate and save summary report
    report_path = OUTPUT_DIR / "summary_report.txt"
    with open(report_path, 'w') as f:
        f.write("="*70 + "\n")
        f.write("RIBOSOMAL PROTEIN COVERAGE ANALYSIS - SUMMARY REPORT\n")
        f.write("="*70 + "\n\n")
        
        f.write(f"Total unique ribosomal proteins detected: {summary['total_unique_proteins']}\n")
        f.write(f"Number of samples analyzed: {summary['num_samples']}\n")
        f.write(f"Samples: {', '.join(summary['sample_names'])}\n\n")
        
        f.write(f"Proteins detected in all samples: {summary['proteins_detected_all_samples']}\n")
        f.write(f"Average detection rate across proteins: {summary['mean_detection_rate']:.1f}%\n\n")
        
        f.write("Coverage Statistics (across all proteins):\n")
        f.write(f"  Mean coverage: {summary['mean_coverage']:.1f}%\n")
        f.write(f"  Median coverage: {summary['median_coverage']:.1f}%\n")
        f.write(f"  Min coverage: {summary['min_coverage']:.1f}%\n")
        f.write(f"  Max coverage: {summary['max_coverage']:.1f}%\n\n")
        
        f.write("Per-Sample Statistics:\n")
        f.write("-" * 70 + "\n")
        for sample in summary['sample_names']:
            stats = summary['per_sample_stats'][sample]
            f.write(f"{sample:15s}: {stats['proteins_detected']:3d} proteins | ")
            f.write(f"Cov: {stats['mean_coverage']:5.1f}% ± {stats['max_coverage']-stats['min_coverage']:5.1f}%\n")
        f.write("=" * 70 + "\n")
    
    logger.info(f"Saved summary report to: {report_path.name}")


def main():
    """
    Main execution function for the MS analysis pipeline.
    Orchestrates all steps: load → filter → aggregate → visualize → report
    """
    logger.info("\n" + "#"*70)
    logger.info("# MASS SPECTROMETRY RIBOSOMAL PROTEIN ANALYSIS PIPELINE")
    logger.info("#"*70 + "\n")
    
    try:
        # Verify data directory exists
        if not DATA_DIR.exists():
            logger.error(f"Data directory not found: {DATA_DIR}")
            return False
        
        # Step 1: Load all samples
        samples = load_all_samples_raw()
        if not samples:
            logger.error("No samples loaded. Pipeline aborted.")
            return False
        
        # Step 2: Filter for ribosomal proteins
        filtered_samples = filter_ribosomal_proteins(samples)
        if not filtered_samples:
            logger.error("No ribosomal proteins found. Pipeline aborted.")
            return False
        
        # Step 3: Aggregate data
        aggregated_df, long_df, summary = aggregate_data(filtered_samples)
        if aggregated_df.empty:
            logger.error("Aggregation failed. Pipeline aborted.")
            return False
        
        # Step 4: Generate visualizations
        plot_files = generate_visualizations(aggregated_df, long_df)
        logger.info(f"Generated {len(plot_files)} visualization files")
        
        # Step 5: Save results
        save_results(aggregated_df, summary)
        
        # Print summary
        logger.info("\n")
        print_aggregation_summary(summary)
        
        logger.info("\n" + "#"*70)
        logger.info("# PIPELINE COMPLETED SUCCESSFULLY")
        logger.info("#"*70)
        logger.info(f"Output directory: {OUTPUT_DIR}\n")
        
        return True
        
    except Exception as e:
        logger.error(f"Pipeline execution failed: {str(e)}", exc_info=True)
        return False


if __name__ == "__main__":
    success = main()
    exit(0 if success else 1)
