import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib.gridspec as gridspec
from scipy.sparse import issparse
import os

def scatter_gene_expression(fea, Y: np.ndarray, marker_gene_indices: list, marker_gene_name: list, 
                            figsize=(15, 9), subplot_size=None, save_fig=False, save_path=None):
    """
    Create scatter plots of gene expression.

    :param fea: Feature matrix for gene expression.
    :param Y: 2D array of coordinates for scatter plot.
    :param marker_gene_indices: List of indices for marker genes.
    :param marker_gene_name: List of marker gene names.
    :param figsize: Size of the figure.
    :param subplot_size: Number of rows and columns for subplots.
    :param save_fig: Boolean indicating whether to save the figure.
    :param save_path: File path to save the figure.
    :raises ValueError: If save_fig is True and save_path is None.
    """
    if save_fig and not save_path:
        raise ValueError("Please provide a valid save_path to save the figure.")
    if issparse(fea):
        fea = fea.toarray()
    n_cols = 5 if subplot_size is None else subplot_size[1]
    n_rows = int(np.ceil((len(marker_gene_indices) + 1) / n_cols))

    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    total_plots = len(marker_gene_indices)
    
    for ii, gene_idx in enumerate(marker_gene_indices):
        color = fea[:, gene_idx]
        ax = axes.flatten()[ii]
        scatter = ax.scatter(Y[:, 0], Y[:, 1], c=color, s=0.1)
        ax.set_title(marker_gene_name[ii])
        fig.colorbar(scatter, ax=ax)
        ax.axis('off')
    
    # Remove extra axes if any
    for extra_ax in axes.flatten()[total_plots:]:
        fig.delaxes(extra_ax)

    plt.tight_layout()

    if save_fig:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        fig.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved at {save_path}")
    
    plt.show()

def plot_mean_expression_heatmap(fea: np.ndarray, marker_gene_indices: list, marker_gene_name: list, 
                                 cell_type: np.ndarray, save_fig=False, save_path=None):
    """
    Plot mean expression heatmap for marker genes across cell types.

    :param fea: 2D array of gene expression features.
    :param marker_gene_indices: List of indices for marker genes.
    :param marker_gene_name: List of marker gene names.
    :param cell_type: Array of cell type labels.
    :param save_fig: Boolean indicating whether to save the figure.
    :param save_path: File path to save the figure.
    :raises ValueError: If save_fig is True and save_path is None.
    """
    if save_fig and not save_path:
        raise ValueError("Please provide a valid save_path to save the figure.")
    if issparse(fea):
        fea = fea.toarray()
    unique_cell_type = np.unique(cell_type)
    mean_expression = np.array([
        [np.mean(fea[cell_type == ct, idx]) for idx in marker_gene_indices]
        for ct in unique_cell_type
    ])
    scaled_expression = (mean_expression - mean_expression.min()) / (mean_expression.max() - mean_expression.min())

    plt.figure(figsize=(10, 8))
    ax = sns.heatmap(scaled_expression, annot=False, cmap='Blues', 
                     cbar_kws={'label': 'Scaled Expression'})
    ax.set_xticks(np.arange(len(marker_gene_name)) + 0.5)
    ax.set_xticklabels(marker_gene_name, rotation=90, ha='center')
    ax.set_yticks(np.arange(len(unique_cell_type)) + 0.5)
    ax.set_yticklabels(unique_cell_type.astype(int))
    plt.xlabel('Gene')
    plt.ylabel('Cluster')

    if save_fig:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved at {save_path}")

    plt.show()

def plot_clustered_heatmap(fea: np.ndarray, marker_gene_indices: list, marker_gene_name: list, 
                           cell_type: np.ndarray, save_fig=False, save_path=None):
    """
    Plot a clustered heatmap of gene expression.

    :param fea: 2D array of gene expression features.
    :param marker_gene_indices: List of indices for marker genes.
    :param marker_gene_name: List of marker gene names.
    :param cell_type: Array of cell type labels.
    :param save_fig: Boolean indicating whether to save the figure.
    :param save_path: File path to save the figure.
    :raises ValueError: If save_fig is True and save_path is None.
    """
    if issparse(fea):
        fea = fea.toarray()
    if save_fig and not save_path:
        raise ValueError("Please provide a valid save_path to save the figure.")
    cell_type = cell_type.reshape(-1)
    indices = np.argsort(cell_type)
    cluster_sorted = cell_type[indices]
    gene_expression = np.vstack([fea[indices, idx] for idx in marker_gene_indices]).T

    unique_clusters, counts = np.unique(cluster_sorted, return_counts=True)
    cluster_positions = np.insert(np.cumsum(counts), 0, 0)
    cluster_centers = cluster_positions[:-1] + np.diff(cluster_positions) // 2

    fig = plt.figure(figsize=(5, 10))
    gs = gridspec.GridSpec(1, 2, width_ratios=[1, 30])
    
    ax0 = plt.subplot(gs[0])
    ax0.imshow(cluster_sorted.reshape(-1, 1), aspect='auto', interpolation='none', cmap='viridis')
    ax0.set_xticks([])
    ax0.set_yticks(cluster_centers)
    ax0.set_yticklabels(unique_clusters)
    ax0.set_ylabel('Cluster')
    
    ax1 = plt.subplot(gs[1])
    sns.heatmap(gene_expression, annot=False, cmap='plasma', ax=ax1, cbar_kws={"shrink": 0.5})
    ax1.set_xticks(np.arange(len(marker_gene_name)) + 0.5)
    ax1.set_xticklabels(marker_gene_name, rotation=90, ha='center')
    ax1.set_yticks([])
    ax1.set_xlabel('Gene')
    
    plt.tight_layout()

    if save_fig:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved at {save_path}")
    
    plt.show()


def compare_violin_plot_between_classes(fea: np.ndarray, marker_gene_indices: list, marker_gene_name: list, 
                                        cluster1: np.ndarray, cluster2: np.ndarray, class_name=['1', '2'], 
                                        figsize=(10, 3), save_fig=False, save_path=None):
    """
    Compare gene expression distributions between two classes using violin plots.

    :param fea: 2D array of gene expression features.
    :param marker_gene_indices: List of indices for marker genes.
    :param marker_gene_name: List of marker gene names.
    :param cluster1: Array of indices for the first cluster.
    :param cluster2: Array of indices for the second cluster.
    :param class_name: Names of the classes for labeling.
    :param figsize: Size of the figure.
    :param save_fig: Boolean indicating whether to save the figure.
    :param save_path: File path to save the figure.
    :raises ValueError: If save_fig is True and save_path is None.
    """
    if issparse(fea):
        fea = fea.toarray()
    if save_fig and not save_path:
        raise ValueError("Please provide a valid save_path to save the figure.")
    
    celltype = np.concatenate((cluster1, cluster2))
    celltype_name = np.array([class_name[0]] * len(cluster1) + [class_name[1]] * len(cluster2))
    palette = 'Set2'
    
    fig, axes = plt.subplots(1, len(marker_gene_indices), figsize=figsize, squeeze=False)
    
    for ii, gene_idx in enumerate(marker_gene_indices):
        gene_expression = fea[celltype, gene_idx].flatten()
        data = pd.DataFrame({'Cell Type': celltype_name, 'Gene Expression': gene_expression})
        
        sns.violinplot(
            x='Gene Expression', y='Cell Type', data=data, orient='h', hue='Cell Type', 
            palette=palette, ax=axes[0, ii], legend=False
        )
        axes[0, ii].set_xlabel(marker_gene_name[ii])
        
        if ii == 0:
            axes[0, ii].set_ylabel('Cluster Name')
        else:
            axes[0, ii].set_ylabel('')
            axes[0, ii].set_yticks([])
        
        for spine in ['top', 'right', 'bottom', 'left']:
            axes[0, ii].spines[spine].set_visible(False)
        axes[0, ii].tick_params(top=False, bottom=False, left=False, right=False)
        axes[0, ii].set_xticks([])
    
    plt.tight_layout()

    if save_fig:
        os.makedirs(os.path.dirname(save_path), exist_ok=True)
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        print(f"Figure saved at {save_path}")
    
    plt.show()