import matplotlib.pyplot as plt
import numpy as np
import umap
from sklearn.metrics import adjusted_rand_score
import os
from scipy.sparse import issparse

def find_map(label1: np.ndarray, label2: np.ndarray) -> np.ndarray:
    """
    Create a mapping matrix between two sets of labels.

    :param label1: First set of labels.
    :param label2: Second set of labels.
    :return: Mapping matrix indicating overlaps between classes.
    """
    unique_label_1 = np.unique(label1)
    unique_label_2 = np.unique(label2)
    n_class_1 = len(unique_label_1)
    n_class_2 = len(unique_label_2)
    mapping = np.zeros((n_class_1, n_class_2))
    
    for ii in range(n_class_1):
        for jj in range(n_class_2):
            temp1 = np.where(label1 == unique_label_1[ii])[0]
            temp2 = np.where(label2 == unique_label_2[jj])[0]
            temp = np.intersect1d(temp1, temp2)
            mapping[ii, jj] = len(temp)
    
    return mapping != 0

def generate_tree_structured(fea, T, step0 = None, step1 = None, cell_type = None) -> tuple:
    """
    Generate the tree-structured representation of the data.

    :return: Tuple containing various representations and mappings.
    """
    # Check if fea is a 2-dimensional NumPy array
    if not isinstance(fea, np.ndarray) or fea.ndim != 2:
        raise ValueError("Input 'fea' must be a 2-dimensional NumPy array.")
    
    # Check if T is a 2-dimensional NumPy array
    if not isinstance(T, np.ndarray) or T.ndim != 2:
        raise ValueError("Input 'T' must be a 2-dimensional NumPy array.")
    
    # Check if step0 is None or a non-negative integer
    if step0 is not None and (not isinstance(step0, int) or step0 < 0):
        raise ValueError("'step0' must be None or a non-negative integer.")
    
    # Check if step1 is None or a non-negative integer
    if step1 is not None and (not isinstance(step1, int) or step1 < 0):
        raise ValueError("'step1' must be None or a non-negative integer.")
    
    # Check if step0 < step1 when both are not None
    if step0 is not None and step1 is not None and step0 >= step1:
        raise ValueError("'step0' must be less than 'step1' when both are provided.")
    
    # Check if cell_type is provided when step1 is None
    if step1 is None and cell_type is None:
        raise ValueError("'cell_type' must be provided when 'step1' is None.")
    
    # Check if cell_type is a 1-dimensional NumPy array when provided
    if cell_type is not None and (not isinstance(cell_type, np.ndarray) or cell_type.ndim != 1):
        raise ValueError("'cell_type' must be a 1-dimensional NumPy array.")
    umap_model = umap.UMAP(n_neighbors=min(15,fea.shape[0]-1), min_dist=0.1, n_components=2)
    Y_initial = umap_model.fit_transform(fea)
    Y_1, Title_1, index_1 = [], [], []
    
    len_T = T.shape[1]
    if step0 == None:
        for ii in range(len_T):
            label = T[:, ii]
            if len(np.unique(label)) != 1:
                break
        step0 = ii

    label_step0 = T[:, step0]
    initial_dict = {ii: f'Cluster{ii}' for ii in np.unique(label_step0)}
    
    for ii in np.unique(label_step0):
        umap_model = umap.UMAP(n_neighbors=min(15,len(np.where(label_step0==ii)[0])-1), min_dist=0.1, n_components=2)
        Y = umap_model.fit_transform(fea[label_step0 == ii, :])
        Y_1.append(Y)
        Title_1.append(initial_dict[ii])
        index_1.append(np.where(label_step0 == ii)[0])
    
    dict_1 = initial_dict
    Y_all, Title_all, index_all = [], [], []
    
    if step1 == None:
        ARI = np.zeros(T.shape[1])
        for ii in range(T.shape[1]):
            ARI[ii] = adjusted_rand_score(T[:, ii], cell_type)
        step1 = np.argmax(ARI)

    for ii in range(step0 + 1, step1 + 1):
        dict_2 = {}
        label_1 = T[:, ii - 1]
        label_2 = T[:, ii]
        unique_label_2 = np.unique(label_2)
        mapping = find_map(label_2, label_1)
        sum_map = np.sum(mapping, axis=0)
        num_rep = 1
        
        for jj in range(len(unique_label_2)):
            temp = np.where(mapping[jj, :] == 1)[0]
            temp = temp[0]
            
            if sum_map[temp] == 1:
                dict_2[jj] = dict_1[temp]
            else:
                umap_model = umap.UMAP(n_neighbors=min(15,len(np.where(label_2==jj)[0])-1), min_dist=0.1, n_components=2)
                Y = umap_model.fit_transform(fea[label_2 == jj, :])
                dict_2[jj] = dict_1[temp] + str(num_rep)
                Y_all.append(Y)
                index_all.append(np.where(label_2 == jj)[0])
                Title_all.append(dict_2[jj])
                num_rep += 1
        
        dict_1 = dict_2
    
    return Y_initial, label_step0, Y_1, Title_1, Y_all, Title_all, index_1, index_all, step0, step1

def visualize_tree_structured(Y_initial, label_step0, Y_1, Title_1, Y_all, Title_all, index_1, index_all, step0, step1, T, save_fig: bool = False, save_path: str = ''):
    """
    Visualize the tree-structured representation.

    :param figsize: Figure size (default is None, which sets a default size).
    :return: Various representations and mappings.
    """
    if save_fig and not isinstance(save_path, str):
        raise ValueError("'save_path' must be a valid string when 'save_fig' is True.")
    np.random.seed(0)
    celltype = T[:, step1].astype(int)
    unique_celltypes = np.unique(celltype)
    colors = np.random.rand(len(unique_celltypes), 3)
    color_map = {label: colors[i] for i, label in enumerate(unique_celltypes)}
    colors_array = np.array([color_map[label] for label in celltype])
    
    if save_fig:
        os.makedirs(save_path, exist_ok=True)
        plt.figure(figsize=(6, 6))
        plt.scatter(Y_initial[:, 0], Y_initial[:, 1], c=colors_array, s=1)
        plt.title('Initial')
        plt.axis('off')
        plt.savefig(os.path.join(save_path, 'Initial.pdf'), dpi=300, bbox_inches='tight', format='pdf')
        plt.close()
        for ii in range(len(Y_1)):
            Y = Y_1[ii]
            n_points = len(Y)
            point_size = 10 if n_points < 200 else (5 if n_points < 500 else 1)
            plt.figure(figsize=(6, 6))
            plt.scatter(Y[:, 0], Y[:, 1], c=colors_array[index_1[ii]], s=point_size)
            plt.title(Title_1[ii])
            plt.axis('off')
            plt.savefig(os.path.join(save_path, Title_1[ii]+'.pdf'), dpi=300, bbox_inches='tight', format='pdf')
            plt.close()
        for ii in range(len(Y_all)):
            Y = Y_all[ii]
            n_points = len(Y)
            point_size = 10 if n_points < 200 else (5 if n_points < 500 else 1)
            plt.figure(figsize=(6, 6))
            plt.scatter(Y[:, 0], Y[:, 1], c=colors_array[index_all[ii]], s=point_size)
            plt.title(Title_all[ii])
            plt.axis('off')
            plt.savefig(os.path.join(save_path, Title_all[ii]+'.pdf'), dpi=300, bbox_inches='tight', format='pdf')
            plt.close()

    total_plots = 1 + len(Y_1) + len(Y_all)
    n_cols = 3
    n_rows = (total_plots + n_cols - 1) // n_cols
    figsize = (15, n_rows * 5)
    
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    axes[0, 0].scatter(Y_initial[:, 0], Y_initial[:, 1], c=colors_array, s=1)
    axes[0, 0].set_title('Initial')
    axes[0, 0].axis('off')
    
    for ii in range(len(Y_1)):
        Y = Y_1[ii]
        ax = axes[(ii + 1) // n_cols, (ii + 1) % n_cols]
        n_points = len(Y)
        point_size = 10 if n_points < 200 else (5 if n_points < 500 else 1)
        ax.scatter(Y[:, 0], Y[:, 1], c=colors_array[index_1[ii]], s=point_size)
        ax.set_title(Title_1[ii])
        ax.axis('off')
    
    for ii in range(len(Y_all)):
        Y = Y_all[ii]
        ax = axes[(len(Y_1) + ii + 1) // n_cols, (len(Y_1) + ii + 1) % n_cols]
        n_points = len(Y)
        point_size = 10 if n_points < 200 else (5 if n_points < 500 else 1)
        ax.scatter(Y[:, 0], Y[:, 1], c=colors_array[index_all[ii]], s=point_size)
        ax.set_title(Title_all[ii])
        ax.axis('off')
    
    for ax in axes.flatten()[total_plots:]:
        fig.delaxes(ax)
    
    plt.tight_layout()
    plt.show()
    
    return 