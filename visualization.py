"""
Visualization Module for Ribosomal Protein Coverage Analysis

This module creates comparative scatter plots showing coverage of ribosomal proteins
across multiple samples, highlighting the heterogeneity between samples.
"""

import logging
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
from typing import Optional, List
from config import OUTPUT_DIR, PLOT_DPI, PLOT_FORMAT, FIGURE_SIZE, SCATTER_SIZE

logger = logging.getLogger(__name__)

# Set seaborn style
sns.set_style("whitegrid")
plt.rcParams['figure.dpi'] = PLOT_DPI


def create_scatter_plot(long_df: pd.DataFrame,
                       title: str = "Ribosomal Protein Coverage Across Samples",
                       output_file: Optional[str] = None,
                       figsize: tuple = FIGURE_SIZE) -> Optional[Path]:
    """
    Create a scatter plot showing protein coverage across samples.
    
    Each protein gets a position on the X-axis, and coverage values from different
    samples are shown as colored points, revealing heterogeneity between samples.
    
    Args:
        long_df: Long-format DataFrame with columns: Gene_Symbol, Sample, Coverage
        title: Title for the plot
        output_file: Optional filename (without extension) to save the plot
        figsize: Figure size as (width, height)
        
    Returns:
        Path to the saved figure, or None if not saved
    """
    
    if long_df.empty:
        logger.error("Cannot create plot from empty DataFrame")
        return None
    
    try:
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Get unique samples and create color palette
        samples = sorted(long_df['Sample'].unique())
        colors = sns.color_palette("husl", len(samples))
        sample_colors = dict(zip(samples, colors))
        
        logger.info(f"Creating scatter plot with {len(samples)} samples and "
                   f"{long_df['Gene_Symbol'].nunique()} proteins")
        
        # Get unique proteins in order (sorted)
        proteins = sorted(long_df['Gene_Symbol'].unique())
        protein_positions = {protein: i for i, protein in enumerate(proteins)}
        
        # Plot each sample with different color
        for sample in samples:
            sample_data = long_df[long_df['Sample'] == sample]
            
            # Convert protein names to numeric positions with slight jitter
            x_positions = sample_data['Gene_Symbol'].map(protein_positions).values
            x_positions_jittered = x_positions + np.random.normal(0, 0.08, size=len(x_positions))
            
            ax.scatter(x_positions_jittered, 
                      sample_data['Coverage'].values,
                      c=[sample_colors[sample]],
                      label=sample,
                      s=SCATTER_SIZE,
                      alpha=0.7,
                      edgecolors='black',
                      linewidth=0.5)
        
        # Formatting
        ax.set_xlabel('Ribosomal Proteins', fontsize=12, fontweight='bold')
        ax.set_ylabel('Coverage [%]', fontsize=12, fontweight='bold')
        ax.set_title(title, fontsize=14, fontweight='bold', pad=20)
        
        # Set X-axis ticks at protein positions
        ax.set_xticks(range(len(proteins)))
        ax.set_xticklabels(proteins, rotation=90, fontsize=9)
        
        # Set Y-axis limits
        ax.set_ylim(-5, 105)
        ax.set_xlim(-1, len(proteins))
        
        # Add legend
        ax.legend(title='Samples', bbox_to_anchor=(1.05, 1), loc='upper left', 
                 frameon=True, fontsize=9, title_fontsize=10)
        
        # Add grid
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        
        # Save figure
        if output_file:
            output_path = OUTPUT_DIR / f"{output_file}.{PLOT_FORMAT}"
            plt.savefig(output_path, dpi=PLOT_DPI, format=PLOT_FORMAT, bbox_inches='tight')
            logger.info(f"Saved plot to {output_path}")
            plt.close()
            return output_path
        else:
            plt.show()
            return None
            
    except Exception as e:
        logger.error(f"Error creating scatter plot: {str(e)}", exc_info=True)
        plt.close()
        return None


def create_faceted_scatter_plot(long_df: pd.DataFrame,
                               proteins_per_plot: int = 20,
                               output_prefix: str = "coverage_plot",
                               figsize: tuple = (14, 8)) -> List[Path]:
    """
    Create multiple faceted scatter plots if there are many proteins.
    
    Splits proteins into groups and creates one plot per group to avoid
    overcrowding when there are many proteins.
    
    Args:
        long_df: Long-format DataFrame with columns: Gene_Symbol, Sample, Coverage
        proteins_per_plot: Number of proteins to display per plot
        output_prefix: Prefix for output filenames
        figsize: Figure size as (width, height)
        
    Returns:
        List of paths to saved figures
    """
    
    proteins = sorted(long_df['Gene_Symbol'].unique())
    samples = sorted(long_df['Sample'].unique())
    
    output_files = []
    
    # Create plots for each group of proteins
    num_plots = (len(proteins) + proteins_per_plot - 1) // proteins_per_plot
    
    logger.info(f"Creating {num_plots} faceted plots with ~{proteins_per_plot} proteins each")
    
    for plot_num in range(num_plots):
        start_idx = plot_num * proteins_per_plot
        end_idx = min((plot_num + 1) * proteins_per_plot, len(proteins))
        
        proteins_subset = proteins[start_idx:end_idx]
        
        # Filter data for this subset
        mask = long_df['Gene_Symbol'].isin(proteins_subset)
        subset_df = long_df[mask].copy()
        
        # Create figure
        fig, ax = plt.subplots(figsize=figsize)
        
        # Get color palette
        colors = sns.color_palette("husl", len(samples))
        sample_colors = dict(zip(samples, colors))
        
        # Map proteins to positions
        protein_positions = {protein: i for i, protein in enumerate(proteins_subset)}
        
        # Plot each sample
        for sample in samples:
            sample_data = subset_df[subset_df['Sample'] == sample]
            
            x_positions = sample_data['Gene_Symbol'].map(protein_positions).values
            x_positions_jittered = x_positions + np.random.normal(0, 0.08, size=len(x_positions))
            
            ax.scatter(x_positions_jittered,
                      sample_data['Coverage'].values,
                      c=[sample_colors[sample]],
                      label=sample,
                      s=SCATTER_SIZE,
                      alpha=0.7,
                      edgecolors='black',
                      linewidth=0.5)
        
        # Formatting
        ax.set_xlabel('Ribosomal Proteins', fontsize=12, fontweight='bold')
        ax.set_ylabel('Coverage [%]', fontsize=12, fontweight='bold')
        
        plot_title = f"Ribosomal Protein Coverage - Plot {plot_num + 1}/{num_plots} " \
                    f"({proteins_subset[0]} to {proteins_subset[-1]})"
        ax.set_title(plot_title, fontsize=13, fontweight='bold', pad=15)
        
        ax.set_xticks(range(len(proteins_subset)))
        ax.set_xticklabels(proteins_subset, rotation=45, ha='right', fontsize=10)
        
        ax.set_ylim(-5, 105)
        ax.set_xlim(-1, len(proteins_subset))
        
        ax.legend(title='Samples', bbox_to_anchor=(1.05, 1), loc='upper left',
                 frameon=True, fontsize=9, title_fontsize=10)
        
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        
        # Save
        output_file = f"{output_prefix}_part{plot_num + 1}"
        output_path = OUTPUT_DIR / f"{output_file}.{PLOT_FORMAT}"
        plt.savefig(output_path, dpi=PLOT_DPI, format=PLOT_FORMAT, bbox_inches='tight')
        logger.info(f"Saved plot part {plot_num + 1} to {output_path}")
        plt.close()
        
        output_files.append(output_path)
    
    return output_files


def create_boxplot(long_df: pd.DataFrame,
                  output_file: Optional[str] = None,
                  figsize: tuple = (14, 6)) -> Optional[Path]:
    """
    Create a boxplot showing coverage distribution for each protein across samples.
    
    Args:
        long_df: Long-format DataFrame with columns: Gene_Symbol, Sample, Coverage
        output_file: Optional filename (without extension) to save the plot
        figsize: Figure size as (width, height)
        
    Returns:
        Path to the saved figure, or None if not saved
    """
    
    if long_df.empty:
        logger.error("Cannot create boxplot from empty DataFrame")
        return None
    
    try:
        fig, ax = plt.subplots(figsize=figsize)
        
        proteins = sorted(long_df['Gene_Symbol'].unique())
        
        logger.info(f"Creating boxplot for {len(proteins)} proteins")
        
        # Prepare data for boxplot
        boxplot_data = [long_df[long_df['Gene_Symbol'] == protein]['Coverage'].values 
                       for protein in proteins]
        
        # Create boxplot
        bp = ax.boxplot(boxplot_data, labels=proteins, patch_artist=True)
        
        # Color the boxes
        colors = sns.color_palette("pastel", 1)
        for patch in bp['boxes']:
            patch.set_facecolor(colors[0])
        
        ax.set_xlabel('Ribosomal Proteins', fontsize=12, fontweight='bold')
        ax.set_ylabel('Coverage [%]', fontsize=12, fontweight='bold')
        ax.set_title('Coverage Distribution of Ribosomal Proteins Across Samples',
                    fontsize=13, fontweight='bold', pad=15)
        
        ax.set_xticklabels(proteins, rotation=90, fontsize=9)
        ax.set_ylim(-5, 105)
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        
        # Save
        if output_file:
            output_path = OUTPUT_DIR / f"{output_file}.{PLOT_FORMAT}"
            plt.savefig(output_path, dpi=PLOT_DPI, format=PLOT_FORMAT, bbox_inches='tight')
            logger.info(f"Saved boxplot to {output_path}")
            plt.close()
            return output_path
        else:
            plt.show()
            return None
            
    except Exception as e:
        logger.error(f"Error creating boxplot: {str(e)}", exc_info=True)
        plt.close()
        return None


def create_heatmap(aggregated_df: pd.DataFrame,
                   output_file: Optional[str] = None,
                   figsize: tuple = (14, 10)) -> Optional[Path]:
    """
    Create a heatmap showing coverage values for all proteins and samples.
    
    Args:
        aggregated_df: Wide-format aggregated DataFrame (from analysis.aggregate_ribosomal_proteins)
        output_file: Optional filename (without extension) to save the plot
        figsize: Figure size as (width, height)
        
    Returns:
        Path to the saved figure, or None if not saved
    """
    
    if aggregated_df.empty:
        logger.error("Cannot create heatmap from empty DataFrame")
        return None
    
    try:
        # Find sample columns (exclude metadata)
        metadata_cols = {'Gene_Symbol', 'Description', 'Mean_Coverage', 'Std_Coverage',
                        'Min_Coverage', 'Max_Coverage', 'Detected_Samples', 'Detection_Rate'}
        sample_cols = [col for col in aggregated_df.columns if col not in metadata_cols]
        sample_cols = sorted(sample_cols)
        
        # Prepare data for heatmap
        heatmap_data = aggregated_df[sample_cols].copy()
        heatmap_data.index = aggregated_df['Gene_Symbol']
        
        logger.info(f"Creating heatmap with {len(heatmap_data)} proteins and {len(sample_cols)} samples")
        
        fig, ax = plt.subplots(figsize=figsize)
        
        # Create heatmap
        sns.heatmap(heatmap_data, cmap='YlOrRd', cbar_kws={'label': 'Coverage [%]'},
                   ax=ax, linewidths=0.5, linecolor='gray')
        
        ax.set_xlabel('Samples', fontsize=12, fontweight='bold')
        ax.set_ylabel('Ribosomal Proteins', fontsize=12, fontweight='bold')
        ax.set_title('Coverage Heatmap: Ribosomal Proteins Across Samples',
                    fontsize=13, fontweight='bold', pad=15)
        
        plt.tight_layout()
        
        # Save
        if output_file:
            output_path = OUTPUT_DIR / f"{output_file}.{PLOT_FORMAT}"
            plt.savefig(output_path, dpi=PLOT_DPI, format=PLOT_FORMAT, bbox_inches='tight')
            logger.info(f"Saved heatmap to {output_path}")
            plt.close()
            return output_path
        else:
            plt.show()
            return None
            
    except Exception as e:
        logger.error(f"Error creating heatmap: {str(e)}", exc_info=True)
        plt.close()
        return None
