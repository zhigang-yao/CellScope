import numpy as np
from scipy.stats import mannwhitneyu, ttest_ind
from scipy.sparse import csr_matrix, issparse
from concurrent.futures import ThreadPoolExecutor, as_completed


def FindMarker(fea: np.ndarray or csr_matrix, cluster1: np.ndarray, cluster2: np.ndarray, 
                            selected_number: int = 5, selected_method: str = 'diff pct') -> tuple:
    """
    Perform differential expression analysis between two clusters.

    :param cluster1: Indices of the first cluster.
    :param cluster2: Indices of the second cluster.
    :param selected_number: Number of top markers to select.
    :param selected_method: Method for selecting markers ('diff pct', 'diff mean', 'FC', 'Wilcoxon', 't-test').
    :return: Tuple of indices of selected markers and their corresponding values.
    """

    if not isinstance(fea, (np.ndarray, csr_matrix)):
        raise ValueError("Input 'fea' must be a NumPy array or CSR matrix.")

    # Check if fea is 2-dimensional
    if fea.ndim != 2:
        raise ValueError("Input 'fea' must be a 2-dimensional matrix.")

    # Check if cluster1 and cluster2 are NumPy arrays
    if not isinstance(cluster1, np.ndarray) or not isinstance(cluster2, np.ndarray):
        raise ValueError("Inputs 'cluster1' and 'cluster2' must be NumPy arrays.")

    # Check if selected_number is a positive integer
    if not (isinstance(selected_number, int) and selected_number > 0):
        raise ValueError("'selected_number' must be a positive integer.")

    # Check if selected_method is a valid string
    valid_methods = {'diff pct', 'diff mean', 'FC', 'Wilcoxon', 't-test'}
    if selected_method not in valid_methods:
        raise ValueError(f"'selected_method' must be one of {valid_methods}.")
    fea_1 = fea[cluster1, :]
    fea_2 = fea[cluster2, :]

    # Convert sparse matrices to dense arrays if needed
    if issparse(fea_1):
        fea_1 = fea_1.toarray()
        fea_2 = fea_2.toarray()

    if selected_method == 'diff pct':
        return _diff_pct(fea_1, fea_2, selected_number)

    elif selected_method == 'diff mean':
        return _diff_mean(fea_1, fea_2, selected_number)

    elif selected_method == 'FC':
        return _fold_change(fea_1, fea_2, selected_number)

    elif selected_method == 'Wilcoxon':
        return _wilcoxon_parallel(fea_1, fea_2, selected_number)

    elif selected_method == 't-test':
        return _t_test_parallel(fea_1, fea_2, selected_number)

def _diff_pct(fea_1: np.ndarray, fea_2: np.ndarray, selected_number: int) -> tuple:
    """
    Calculate the differential percentage of expression between two clusters.

    :param fea_1: Feature matrix for the first cluster.
    :param fea_2: Feature matrix for the second cluster.
    :param selected_number: Number of top markers to select.
    :return: Tuple of indices of selected markers and their corresponding values.
    """
    pct_1 = np.sum(fea_1 != 0, axis=0) / fea_1.shape[0]
    pct_2 = np.sum(fea_2 != 0, axis=0) / fea_2.shape[0]
    Res = pct_1 - pct_2
    return _get_sorted_indices(Res, selected_number)

def _diff_mean(fea_1: np.ndarray, fea_2: np.ndarray, selected_number: int) -> tuple:
    """
    Calculate the differential mean expression between two clusters.

    :param fea_1: Feature matrix for the first cluster.
    :param fea_2: Feature matrix for the second cluster.
    :param selected_number: Number of top markers to select.
    :return: Tuple of indices of selected markers and their corresponding values.
    """
    mean_1 = np.mean(fea_1, axis=0)
    mean_2 = np.mean(fea_2, axis=0)
    Res = mean_1 - mean_2
    return _get_sorted_indices(Res, selected_number)

def _fold_change(fea_1: np.ndarray, fea_2: np.ndarray, selected_number: int) -> tuple:
    """
    Calculate fold change (FC) between two clusters.

    :param fea_1: Feature matrix for the first cluster.
    :param fea_2: Feature matrix for the second cluster.
    :param selected_number: Number of top markers to select.
    :return: Tuple of indices of selected markers and their corresponding FC values.
    """
    mean_1 = np.mean(fea_1, axis=0)
    mean_2 = np.mean(fea_2, axis=0)
    FC = mean_1 / mean_2
    case1 = np.where(mean_1 == 0)[0]
    case2 = np.where(mean_2 == 0)[0]
    FC = np.log2(FC + 1)
    FC[case1] = -100
    FC[case2] = 100
    FC[np.intersect1d(case1, case2)] = 0
    return np.argsort(np.abs(FC))[-selected_number:], FC

def _wilcoxon_parallel(fea_1: np.ndarray, fea_2: np.ndarray, selected_number: int) -> tuple:
    """
    Perform Wilcoxon rank-sum test in parallel.

    :param fea_1: Feature matrix for the first cluster.
    :param fea_2: Feature matrix for the second cluster.
    :param selected_number: Number of top markers to select.
    :return: Tuple of indices of selected markers and their corresponding p-values.
    """
    Res = np.zeros(fea_1.shape[1])
    with ThreadPoolExecutor() as executor:
        futures = {executor.submit(mannwhitneyu, fea_1[:, ii], fea_2[:, ii], alternative='two-sided', use_continuity=True): ii for ii in range(fea_1.shape[1])}
        for future in as_completed(futures):
            ii = futures[future]
            _, Res[ii] = future.result()
    return np.argsort(Res)[:selected_number], Res

def _t_test_parallel(fea_1: np.ndarray, fea_2: np.ndarray, selected_number: int) -> tuple:
    """
    Perform t-test in parallel.

    :param fea_1: Feature matrix for the first cluster.
    :param fea_2: Feature matrix for the second cluster.
    :param selected_number: Number of top markers to select.
    :return: Tuple of indices of selected markers and their corresponding p-values.
    """
    Res = np.zeros(fea_1.shape[1])
    with ThreadPoolExecutor() as executor:
        futures = {executor.submit(ttest_ind, fea_1[:, ii], fea_2[:, ii]): ii for ii in range(fea_1.shape[1])}
        for future in as_completed(futures):
            ii = futures[future]
            _, Res[ii] = future.result()
    return np.argsort(np.abs(Res))[:selected_number], Res

def _get_sorted_indices(Res: np.ndarray, selected_number: int) -> tuple:
    """
    Get the indices of the top selected markers based on results.

    :param Res: Results array from differential analysis.
    :param selected_number: Number of top markers to select.
    :return: Indices of selected markers and the corresponding results.
    """
    if Res.shape[0] == 1:
        return np.argsort(np.abs(Res))[0, -selected_number:], Res
    else:
        return np.argsort(np.abs(Res))[-selected_number:], Res
