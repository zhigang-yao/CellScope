from scipy.stats import wasserstein_distance
from joblib import Parallel, delayed
import numpy as np
import plotly.graph_objects as go
import pandas as pd
def Gene_Analysis(fea, layer):
    """
    Perform gene analysis using Wasserstein distance.

    :param fea: Feature matrix (2D NumPy array).
    :param layer: List or array indicating clusters.
    :return: Wasserstein distances, labels, label strings, and flow labels.
    :raises ValueError: If input parameters or data are invalid.
    """
    # Check if fea is a 2-dimensional NumPy array
    if not isinstance(fea, np.ndarray) or fea.ndim != 2:
        raise ValueError("Input 'fea' must be a 2-dimensional NumPy array.")
    
    # Check if len(layer) is even and greater than 0
    if len(layer) == 0 or len(layer) % 2 != 0:
        raise ValueError("Input 'layer' must have an even number of elements and must not be empty.")

    len_flow = int(len(layer) / 2)
    num_Gene = fea.shape[1]
    Res = np.zeros((num_Gene, len_flow))

    def calculate_wasserstein(ii, jj):
        cluster1 = layer[ii * 2]
        cluster2 = layer[ii * 2 + 1]
        Gene1 = fea[cluster1, jj]
        Gene2 = fea[cluster2, jj]
        return wasserstein_distance(Gene1, Gene2)

    # Calculate Wasserstein distances in parallel
    for ii in range(len_flow):
        Res[:, ii] = Parallel(n_jobs=-1)(delayed(calculate_wasserstein)(ii, jj) for jj in range(num_Gene))

    # Define gene labels based on Wasserstein distances
    label = np.zeros_like(Res).astype(int)
    label[Res < 0.5] = 0 
    label[(Res >= 0.5) & (Res <= 1)] = 1 
    label[Res > 1] = 2 

    # Define label strings
    label_str = np.zeros_like(label).astype(object)
    label_str[label == 0] = 'Housekeeper Gene'
    label_str[label == 1] = 'Type-Related Gene'
    label_str[label == 2] = 'Type-Determining Gene'

    # Create flow labels
    unique_samples, counts = np.unique(label, axis=0, return_counts=True)
    flow_labels = np.zeros(label.shape[0], dtype=object) 
    for i, sample in enumerate(unique_samples):
        flow_label = f"Flow {i}"
        indices = np.where((label == sample).all(axis=1))[0]
        flow_labels[indices] = flow_label

    return Res, label, label_str, flow_labels

def plot_sankey(label_str, width=1100, height=600,save_fig=False, save_path='./sankey_diagram.png'):
    """
    Plot a Sankey diagram based on the gene labels across columns in label_str.

    :param label_str: A 2D NumPy array where each column contains gene labels 
                      ('Housekeeper Gene', 'Type-Related Gene', 'Type-Determining Gene').
    """
    # Define all unique labels
    unique_labels = ['Housekeeper Gene', 'Type-Related Gene', 'Type-Determining Gene']
    
    # Initialize source, target, and value lists
    sources = []
    targets = []
    values = []
    
    num_columns = label_str.shape[1]
    gene_counts = pd.DataFrame(columns=[f'Layer {i+1}' for i in range(num_columns)],
                            index=unique_labels)

    for i in range(num_columns):
        col = label_str[:, i]
        for label in unique_labels:
            gene_counts.iloc[unique_labels.index(label), i] = np.sum(col == label)
    # Create mappings from column to column for the Sankey diagram
    for i in range(num_columns - 1):
        col1 = label_str[:, i]
        col2 = label_str[:, i + 1]

        # Iterate through each pair of labels
        for label_from in unique_labels:
            for label_to in unique_labels:
                count = np.sum((col1 == label_from) & (col2 == label_to))
                if count > 0:
                    sources.append(unique_labels.index(label_from) + i * len(unique_labels))
                    targets.append(unique_labels.index(label_to) + (i + 1) * len(unique_labels))
                    values.append(count)

    # Define node labels for the Sankey diagram
    node_labels = []
    for i in range(num_columns):
        node_labels.extend([f'{label} (Layer {i + 1})' for label in unique_labels])
    
    # Create Sankey diagram
    fig = go.Figure(go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=node_labels
        ),
        link=dict(
            source=sources,
            target=targets,
            value=values
        )
    ))

    fig.update_layout(width = int(width), height = int(height))
    if save_fig:
        # Check if the save_path has a valid file extension
        if not save_path.endswith(('.png', '.pdf', '.jpeg')):
            raise ValueError("Save path must end with '.png', '.pdf', or '.jpeg'")
        
        # Save the figure using kaleido
        fig.write_image(save_path, engine='kaleido')
        print(f"Figure saved to {save_path}")
    fig.show()
    return gene_counts