import numpy as np
from sklearn.metrics import (accuracy_score, normalized_mutual_info_score, 
                             adjusted_rand_score, f1_score, jaccard_score)
from scipy.optimize import linear_sum_assignment


import numpy as np
from sklearn.metrics import accuracy_score, normalized_mutual_info_score, adjusted_rand_score, f1_score, jaccard_score
from scipy.optimize import linear_sum_assignment

def process_labels(labels: np.ndarray) -> (np.ndarray, dict):
    """
    Convert labels to numpy array format, ensuring they are numeric.

    :param labels: Input labels, can be a list or numpy array.
    :return: Processed labels as a numpy array of integers, and the label mapping dictionary.
    :raises ValueError: If labels are not of valid type or contain invalid values.
    """
    # Ensure labels is a numpy array
    if not isinstance(labels, (np.ndarray, list)):
        raise ValueError("Input 'labels' must be a list or a numpy array.")
    
    # Convert list to numpy array if needed
    if not isinstance(labels, np.ndarray):
        labels = np.array(labels)
    
    # Handle non-numeric labels
    if labels.dtype.kind not in {'i', 'f'}:
        unique_labels = np.unique(labels)
        label_mapping = {label: idx for idx, label in enumerate(unique_labels)}
        labels = np.array([label_mapping[label] for label in labels])
    else:
        label_mapping = {}

    # Ensure labels are non-negative integers
    labels = labels - np.min(labels)
    labels = labels.astype(int)
    
    return labels, label_mapping

def best_map(L1: np.ndarray, L2: np.ndarray) -> np.ndarray:
    """
    Permute labels of L2 to match L1 as closely as possible.

    :param L1: True labels.
    :param L2: Predicted labels to be permuted.
    :return: Permuted labels of L2.
    :raises ValueError: If L1 and L2 are not the same length.
    """
    # Ensure L1 and L2 are numpy arrays
    if not isinstance(L1, np.ndarray) or not isinstance(L2, np.ndarray):
        raise ValueError("Inputs 'L1' and 'L2' must be numpy arrays.")
    
    # Ensure L1 and L2 have the same length
    if L1.size != L2.size:
        raise ValueError("Inputs 'L1' and 'L2' must have the same length.")

    # Ensure labels are non-negative integers
    L1 = L1 - np.min(L1)
    L2 = L2 - np.min(L2)

    D = max(L1.max(), L2.max()) + 1
    w = np.zeros((D, D), dtype=np.int64)
    for i in range(L1.size):
        w[L1[i], L2[i]] += 1

    ind = linear_sum_assignment(w.max() - w)
    ind = np.asarray(ind).T
    new_L2 = np.zeros(L2.size, dtype=L2.dtype)
    for i in range(ind.shape[0]):
        new_L2[L2 == ind[i, 1]] = ind[i, 0]

    return new_L2

def calculate_metrics(true_labels: np.ndarray, pred_labels: np.ndarray) -> dict:
    """
    Calculate and return clustering evaluation metrics.

    :param true_labels: True labels of the data.
    :param pred_labels: Predicted labels from the clustering algorithm.
    :return: A dictionary containing various clustering metrics.
    :raises ValueError: If inputs are not of the same length or invalid type.
    """
    # Ensure true_labels and pred_labels are numpy arrays
    if not isinstance(true_labels, np.ndarray) or not isinstance(pred_labels, np.ndarray):
        raise ValueError("Inputs 'true_labels' and 'pred_labels' must be numpy arrays.")

    # Ensure true_labels and pred_labels have the same length
    if true_labels.size != pred_labels.size:
        raise ValueError("Inputs 'true_labels' and 'pred_labels' must have the same length.")

    # Process and map labels
    true_labels, _ = process_labels(true_labels)
    pred_labels, _ = process_labels(pred_labels)

    # Permute pred_labels to best match true_labels
    pred_labels = best_map(true_labels, pred_labels)

    # Calculate clustering metrics
    acc = accuracy_score(true_labels, pred_labels)
    nmi = normalized_mutual_info_score(true_labels, pred_labels)
    ari = adjusted_rand_score(true_labels, pred_labels)
    f_measure = f1_score(true_labels, pred_labels, average='weighted')
    jaccard = jaccard_score(true_labels, pred_labels, average='micro')

    return {
        "Accuracy": acc,
        "NMI": nmi,
        "ARI": ari,
        "F1 Score": f_measure,
        "Jaccard": jaccard
    }
