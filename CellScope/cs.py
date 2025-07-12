import numpy as np
import scipy.stats as stats
from scipy.sparse import csr_matrix, issparse, isspmatrix_csr, coo_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import fcluster
from scipy.stats import f_oneway
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors
from sklearn.utils.extmath import randomized_svd
from numba import njit
import fastcluster
import random
from scipy.sparse import diags
import gc
from scipy.stats import mode
from joblib import Parallel, delayed
from sklearn.manifold import TSNE
from umap.umap_ import UMAP

seed = 83
random.seed(seed)
np.random.seed(seed)

def deterministic_svd(X, n_components, random_state):
    X = X.astype(np.float64)
    U, S, Vt = randomized_svd(
        X,
        n_components=n_components,
        n_iter=5,
        random_state=random_state
    )
    sign_flip = np.sign(Vt[:, 0]).astype(np.int8)
    np.multiply(U, sign_flip[np.newaxis, :], out=U)
    np.multiply(Vt, sign_flip[:, np.newaxis], out=Vt)
    
    return U, S

def compute_umap_similarity_graph(knn_index, knn_dist):
    n_samples = knn_index.shape[0]
    n_neighbors = knn_index.shape[1]
    sigmas, rhos = smooth_knn_dist(knn_dist, n_neighbors)
    rows, cols, vals = compute_membership_strengths(
        knn_index, knn_dist, sigmas, rhos, n_neighbors
    )
    W = coo_matrix((vals, (rows, cols)), shape=(n_samples, n_samples))
    W.eliminate_zeros()
    transpose = W.transpose()
    prod_matrix = W.multiply(transpose)
    W = W + transpose - prod_matrix
    W.eliminate_zeros()
    return W

@njit
def smooth_knn_dist(distances, k, local_connectivity=1.0):
    n_samples = distances.shape[0]
    rhos = np.zeros(n_samples)
    sigmas = np.zeros(n_samples)
    
    for i in range(n_samples):
        rhos[i] = distances[i, 0] if k > 0 else 0.0
        lo = 0.0
        hi = np.inf
        mid = 1.0
        for _ in range(50): 
            psum = 0.0
            for j in range(1, k):
                d = distances[i, j] - rhos[i]
                if d > 0:
                    psum += np.exp(-(d / mid))
                else:
                    psum += 1.0
            if np.abs(psum - local_connectivity) < 1e-5:
                break
                
            if psum > local_connectivity:
                hi = mid
                mid = (lo + hi) / 2
            else:
                lo = mid
                if hi == np.inf:
                    mid *= 2
                else:
                    mid = (lo + hi) / 2
                    
        sigmas[i] = mid
    
    return sigmas, rhos

@njit
def compute_membership_strengths(knn_index, knn_dist, sigmas, rhos, k):
    n_samples = knn_index.shape[0]
    rows = np.zeros(n_samples * k, dtype=np.int32)
    cols = np.zeros(n_samples * k, dtype=np.int32)
    vals = np.zeros(n_samples * k, dtype=np.float32)
    
    position = 0
    for i in range(n_samples):
        for j in range(k):
            if knn_index[i, j] == -1:  # 无效邻居
                continue
            dist = max(0.0, knn_dist[i, j] - rhos[i])
            val = np.exp(-dist / sigmas[i])
            rows[position] = i
            cols[position] = knn_index[i, j]
            vals[position] = val
            position += 1
    return rows[:position], cols[:position], vals[:position]


def large_scale_knn(X, k=5,metric='euclidean'):
    nn = NearestNeighbors(n_neighbors=k, algorithm='auto', 
                         metric=metric, n_jobs=-1)
    nn.fit(X)
    return nn.kneighbors(X)

def knn_search(data,k=5,metric = 'euclidean'):
    knn_dist,knn_index = large_scale_knn(data, k,metric = metric)
    return knn_index,knn_dist
def jaccard_distance(A, B):
    intersection = np.logical_and(A, B).sum()
    union = np.logical_or(A, B).sum()
    if union == 0: 
        return 1.0
    return 1 - (intersection / union)

def construct_graph(fea,metric = 'euclidean',jaccard = False):
    if jaccard:
        fea[fea>0] = 1
        fea = fea.astype(bool)
        knn = 10*int(np.log(fea.shape[0]))
        knn_index,knn_dist = knn_search(fea,k=knn,metric=metric)
        num_cell = knn_index.shape[0]
        knn_dist = np.zeros((num_cell, knn))
        
        for i in range(num_cell):
            vec_i = fea[i,:]
            for k in range(knn):
                j = knn_index[i, k]
                vec_j = fea[j]
                knn_dist[i, k] =jaccard_distance(vec_i, vec_j)
        W = compute_umap_similarity_graph(knn_index, knn_dist)
    else:
        knn = int(np.log(fea.shape[0]))
        knn_index,knn_dist = knn_search(fea,k=knn,metric=metric)
        num_cell = knn_index.shape[0]
        W = compute_umap_similarity_graph(knn_index, knn_dist)
        # neighbors = umap.fuzzy_simplicial_set(
        #     fea,
        #     n_neighbors=knn, 
        #     random_state=83,
        #     metric=metric
        # )
        # W = neighbors[0]
    return W

def Normalization(fea_raw: np.ndarray or csr_matrix) -> (np.ndarray or csr_matrix, np.ndarray or csr_matrix, np.ndarray or csr_matrix):
    """
    Normalize the raw feature data.

    :param fea_raw: Raw feature data as a NumPy array or CSR matrix.
    :return: Tuple of normalized features and log-transformed features.
    """
    if not isinstance(fea_raw, (np.ndarray, csr_matrix)):
        raise ValueError("Input must be a NumPy array or CSR matrix.")

    if issparse(fea_raw):
        if fea_raw.min() < 0:
            raise ValueError("Sparse matrix contains negative values, which is not allowed.")
    else:
        if np.any(fea_raw < 0):
            raise ValueError("Array contains negative values, which is not allowed.")
        
    if issparse(fea_raw):
        fea_log = csr_matrix(fea_raw)
        if np.max(fea_raw.data) <= 1000:
            fea_log = fea_raw
        else:
            fea_log.data = np.log2(fea_raw.data + 1)
        
        row_squared_sum = np.array(fea_log.power(2).sum(axis=1)).flatten()
        row_scaling_factors = np.sqrt(row_squared_sum)
        row_scaling_factors[row_scaling_factors == 0] = 1
        row_inv_scaling_factors = 1.0 / row_scaling_factors
        
        n = fea_log.shape[0]
        row_indices = np.arange(n)
        row_inv_diag = csr_matrix((row_inv_scaling_factors, (row_indices, row_indices)), shape=(n, n))
        fea = row_inv_diag.dot(fea_log)

    else:
        if np.max(fea_raw) <= 1000:
            fea_log = fea_raw
            fea_raw = 2 ** (fea_raw) - 1
        else:
            fea_log = np.log2(fea_raw + 1)
        
        fea = fea_log / np.sqrt(np.sum(fea_log ** 2, axis=1)[:, np.newaxis])
        fea[np.isnan(fea)] = 0

    return fea_raw, fea_log, fea


def findCenters(rho: np.ndarray, delta: np.ndarray) -> (np.ndarray, np.ndarray):
    """
    Find centers based on the density and distance metrics.

    :param rho: Density values.
    :param delta: Distance values.
    :return: Centers and delta-rho values.
    :raises ValueError: If the input arrays are not 1-dimensional, 
                        or if they contain negative values, 
                        or if they have zero variance.
    """
    # Check if inputs are NumPy arrays
    if not isinstance(rho, np.ndarray) or not isinstance(delta, np.ndarray):
        raise ValueError("Inputs must be NumPy arrays.")

    # Check if inputs are 1-dimensional
    if rho.ndim != 1 or delta.ndim != 1:
        raise ValueError("Inputs must be 1-dimensional arrays.")

    # Check if rho and delta have the same length
    if rho.shape[0] != delta.shape[0]:
        raise ValueError("Inputs 'rho' and 'delta' must have the same length.")
    
    # Check for negative values in rho and delta
    if np.any(rho < 0) or np.any(delta < 0):
        raise ValueError("Inputs 'rho' and 'delta' cannot contain negative values.")

    # Check for zero variance in rho and delta
    if np.var(rho) == 0 or np.var(delta) == 0:
        raise ValueError("Inputs 'rho' and 'delta' must have non-zero variance.")

    # Normalize rho and delta
    rho = (rho - np.min(rho)) / (np.max(rho) - np.min(rho))
    delta = (delta - np.min(delta)) / (np.max(delta) - np.min(delta))
    
    # Calculate deltarho and sort it
    deltarho = delta * rho
    sorted_deltarho = np.sort(deltarho)[::-1]
    idx = np.argsort(deltarho)[::-1]
    
    # Calculate gradient and find potential centers
    g1 = np.gradient(sorted_deltarho)
    temp1 = idx[g1 < np.mean(g1)]
    temp2 = np.where(deltarho > np.mean(deltarho))[0]
    temp3 = np.where(delta > np.mean(delta))[0]
    temp4 = np.where(rho > np.mean(rho))[0]
    temp1 = np.intersect1d(temp1, temp4)
    centers = np.intersect1d(temp1, np.intersect1d(temp2, temp3))
    
    return centers, deltarho

def calculate_p_value(i: int, fea: np.ndarray or csr_matrix, temp: np.ndarray, label: np.ndarray) -> float:
    """
    Calculate p-values for the Genes using ANOVA.

    :param i: Index of the Gene.
    :param fea: Nomalization Feature matrix.
    :param temp: Temporary array for indexing.
    :param label: Labels for different classes.
    :return: p-value for the Gene.
    """
    if isspmatrix_csr(fea):
        data = [fea[temp[label==lbl], i].toarray().reshape(-1) for lbl in np.unique(label)]
    else:
        data = [fea[temp[label==lbl], i] for lbl in np.unique(label)]
    _, p_value = f_oneway(*data)
    return p_value


def Manifold_Fitting_1(fea: np.ndarray or csr_matrix, num_pca: int = 100, num_Selected_Gene: int = 500,
                        knn: int = 20, num_center: int = 0, random_seed: int = 83) -> (np.ndarray, np.ndarray, np.ndarray):
    """
    Perform manifold fitting to identify Signal Space.

    :param fea: Feature matrix.
    :param num_pca: Number of PCA components.
    :param num_Selected_Gene: Number of selected genes.
    :param knn: Number of nearest neighbors.
    :param num_center: Number of centers to find.
    :return: Tuple of selected features, p-values, and indices of significant features.
    """
    if not isinstance(fea, (np.ndarray, csr_matrix)):
        raise ValueError("Input 'fea' must be a NumPy array or CSR matrix.")
    
    # Check if fea has at least 2 dimensions
    if fea.ndim != 2:
        raise ValueError("Input 'fea' must be a 2-dimensional matrix.")
    
    # Check if num_pca, num_Selected_Gene, and knn are positive integers
    if not (isinstance(num_pca, int) and num_pca > 0):
        raise ValueError("'num_pca' must be a positive integer.")
    if not (isinstance(num_Selected_Gene, int) and num_Selected_Gene > 0):
        raise ValueError("'num_Selected_Gene' must be a positive integer.")
    if not (isinstance(knn, int) and knn > 0):
        raise ValueError("'knn' must be a positive integer.")
    
    # Check if num_center is non-negative
    if not (isinstance(num_center, int) and num_center >= 0):
        raise ValueError("'num_center' must be a non-negative integer.")
    num_cell,num_gene = fea.shape
    num_pca = min(num_pca,num_cell-1)
    if issparse(fea):
        U, S = deterministic_svd(fea, num_pca, random_seed)
        fea_pca = U * S
    else:
        pca = PCA(n_components=num_pca,random_state=random_seed)
        fea_pca = pca.fit_transform(fea)

    nbrs = NearestNeighbors(n_neighbors=knn, algorithm='auto',n_jobs = -1).fit(fea_pca)
    D_NB,ID_NB = nbrs.kneighbors(fea_pca)
    
    # Select the top k nearest neighbor distances and indices
    rho = 1. / np.sum(D_NB, axis=1)
    delta = np.zeros(num_cell)
    higher_density_mask = rho.reshape(-1, 1) < rho.reshape(1, -1)
    
    for ii in range(num_cell):
        temp = np.where(higher_density_mask[ii])[0]
        inter_temp = np.intersect1d(temp, ID_NB[ii,:])
        if len(inter_temp) == 0:
            delta[ii] = np.max(D_NB[ii,:])
        else:
            ib = np.where(np.isin(ID_NB[ii], inter_temp))[0]
            delta[ii] = np.min(D_NB[ii, ib])
    if num_center == 0:
        id_max_deltarho, deltarho = findCenters(rho, delta)
        id_max_deltarho = np.int32(id_max_deltarho)
        if len(id_max_deltarho) < 3:
            id_max_deltarho = np.argsort(-deltarho)[:3]
        num_center = len(id_max_deltarho)
    else:
        deltarho = delta*rho
        id_max_deltarho = np.argsort(-deltarho)[:num_center]
    Clusters = np.zeros((num_cell, num_center))
    for jj in range(num_center):
        idx = id_max_deltarho[jj]
        Clusters[ID_NB[idx, 0:5], jj] = 1
    clusters = np.max(Clusters, axis=1)
    Clusters = Clusters[clusters > 0]
    Q = np.dot(Clusters, Clusters.T)
    _, label = connected_components(Q, directed=False)
    kk = 5
    while len(np.unique(label))==1:
        Clusters = np.zeros((num_cell, num_center))
        kk -= 1
        for jj in range(num_center):
            idx = id_max_deltarho[jj]
            Clusters[ID_NB[idx, 0:kk], jj] = 1
        clusters = np.max(Clusters, axis=1)
        Clusters = Clusters[clusters > 0]
        Q = np.dot(Clusters, Clusters.T)
        _, label = connected_components(Q, directed=False)
    temp = np.where(clusters>0)[0]
    unique_labels = np.unique(label)
    p_values = np.zeros(num_gene)
    temp = np.where(clusters>0)[0]
    p_values = Parallel(n_jobs=-1)(delayed(calculate_p_value)(i, fea, temp, label) for i in range(num_gene))
    p_values = np.array(p_values)
    sorted_indices = np.argsort(p_values)
    significant_features_index = sorted_indices[:num_Selected_Gene]
    fea_selected = fea[:, significant_features_index]
    return fea_selected, significant_features_index, id_max_deltarho

def Manifold_Fitting_2( fea_selected: np.ndarray, num_neighbor: int = 5, fitting_prop: float = 0.05,
                        coeff: float = 0.1, op_Outlier: bool = False) -> (np.ndarray, np.ndarray, list):
    """
    Perform second manifold fitting and outlier removal.

    :param fea_selected: Selected feature matrix.
    :param num_neighbor: Number of nearest neighbors.
    :param fitting_prop: Proportion of fitting points.
    :param coeff: Coefficient for updating outlier points.
    :param op_Outlier: Option to remove outliers.
    :return: Tuple of new feature matrix, fitting indices, and remaining indices.
    """
    if not (isinstance(num_neighbor, int) and num_neighbor > 0):
        raise ValueError("'num_neighbor' must be a positive integer.")

    # Check if fitting_prop is between 0 and 1
    if not (0 < fitting_prop <= 1):
        raise ValueError("'fitting_prop' must be a float between 0 and 1.")

    # Check if coeff is between 0 and 1
    if not (0 <= coeff <= 1):
        raise ValueError("'coeff' must be a float between 0 and 1.")

    # Check if op_Outlier is a boolean
    if not isinstance(op_Outlier, bool):
        raise ValueError("'op_Outlier' must be a boolean.")
    
    if issparse(fea_selected):
        fea_selected = fea_selected.toarray()
    n = fea_selected.shape[0]
    nbrs = NearestNeighbors(n_neighbors=num_neighbor, algorithm='auto').fit(fea_selected)
    distances, indices = nbrs.kneighbors(fea_selected)
    D_NB = distances[:, :num_neighbor]
    rho = 1. / np.sum(D_NB, axis=1)
    Purity_Matrix = rho
    fitting_index = np.argsort(Purity_Matrix)
    fitting_index = fitting_index[:int(np.floor(fitting_prop*n))]
    fea_new = fea_selected.copy()
    outlier = []
    for ii in range(len(fitting_index)):
        for jj in range(num_neighbor):  
            if indices[fitting_index[ii], jj] not in fitting_index:
                fea_new[fitting_index[ii], :] = coeff*fea_new[fitting_index[ii], :] + (1-coeff)*fea_selected[indices[fitting_index[ii], jj], :]
                break
        if op_Outlier:                
            if jj == num_neighbor - 1:
                outlier.append(fitting_index[ii])
    index= set(range(n)) 
    index = index-set(outlier)
    index = [i for i in index]
    fea_new = fea_new[index,:]
    fea_selected = fea_new
    return fea_selected,fitting_index,index

def graph_importance_sampling(W, num_samples):
    centrality = np.sum(W, axis=1)       # 使用度中心性
    p = centrality / np.sum(centrality)
    return np.random.choice(range(len(W)), num_samples, p=p, replace=False)


def Trans_W_D(W_sparse):
    n = W_sparse.shape[0]
    D_vector = np.ones(n * (n - 1) // 2, dtype=np.float32)  # 全部初始化为1（即 D[i,j]=1）

    W_sparse = W_sparse - diags(W_sparse.diagonal())
    W_sparse = W_sparse.tocoo()

    for i, j, val in zip(W_sparse.row, W_sparse.col, W_sparse.data):
        if i > j:
            idx = j * n - (j * (j + 1)) // 2 + (i - j - 1)
        else:
            continue  # 跳过对角线
        D_vector[idx] = 1 - val  # 更新真实值
    return D_vector

def GraphCluster(fea_selected: np.ndarray, metric: str = None, 
                    num_cell_thre: int = 100000, index: list = [],random_seed = 83) -> np.ndarray:
    """
    Perform graph-based clustering on the selected features.

    :param fea_selected: Selected feature matrix.
    :param metric: Metric for distance calculation (None or others).
    :param num_cell_thre: Threshold for number of cells.
    :param index: Optional index for selected cells.
    :return: Array of cluster labels for all cells.
    """
    if fea_selected.ndim != 2:
        raise ValueError("Input 'fea_selected' must be a 2-dimensional array.")

    # Check if num_cell_thre is a positive integer
    if not (isinstance(num_cell_thre, int) and num_cell_thre > 0):
        raise ValueError("'num_cell_thre' must be a positive integer.")
    
    # Check if index is a list
    if not isinstance(index, list):
        raise ValueError("'index' must be a list.")
    
    # Check if metric is a string
    if metric!=None:
        if not isinstance(metric, str):
            raise ValueError("'metric' must be a string.")
        
    # Check if random_seed is an integer
    if not isinstance(random_seed, int):
        raise ValueError("'random_seed' must be an integer.")
    if fea_selected.shape[0]<10000:
        jaccard = False
    else:
        jaccard = True
    num_cell = fea_selected.shape[0]
    if num_cell>num_cell_thre:
        if index==[]:
            random.seed(random_seed) 
            index = random.sample(range(num_cell), num_cell_thre)
        
        index_unselected = set(range(num_cell))-set(index)
        index_unselected = list(index_unselected)
        num_remaining_samples = len(index_unselected)
        fea_selected_index = fea_selected[index,:]
        fea_selected_unindex = fea_selected[index_unselected,:]
        if 'fea_selected' in locals() and jaccard:
            del fea_selected
            gc.collect()
        nn = NearestNeighbors(n_neighbors=5, metric=metric)
        nn.fit(fea_selected_index)
        knn_dist, knn_index = nn.kneighbors(fea_selected_unindex)
        W1 = construct_graph(fea_selected_index ,jaccard = jaccard)
        
        W1 = Trans_W_D(W1)
        Z = fastcluster.linkage(W1, method='average')
        del W1; gc.collect()
        T_all = np.zeros((num_cell,49))
        
        for ii in range(2,51):
            T = 100*np.ones((num_cell))
            T1 = fcluster(Z, ii, criterion='maxclust')
            T1 = T1-1
            T[index] = T1
            knn_classes = T1[knn_index]  # 获取每个近邻的类别 (n_samples_unselected, k)
            remaining_samples_labels, _ = mode(knn_classes, axis=1)
            remaining_samples_labels = remaining_samples_labels.ravel()  # 展平为一维数组
            T[index_unselected] = remaining_samples_labels  # 未选取样本使用KNN多数投票结果
            T_all[:, ii-2] = T
    else:
        W = construct_graph(fea_selected,jaccard = jaccard)
        W = Trans_W_D(W)
        Z = fastcluster.linkage(W, method='average')
        T_all = np.zeros((num_cell,49))
    
        for ii in range(2,51):
            T = fcluster(Z, ii, criterion='maxclust')
            T_all[:,ii-2] = T
    T_all -= np.min(T_all, axis=0)
    T_all = T_all.astype(int)
    return T_all


def Visualization( fea: np.ndarray, Visualization_Method: str = "UMAP", random_seed: int = 83) -> np.ndarray:
    """
    Visualize the feature data using specified visualization method.

    :param fea: Feature matrix to visualize.
    :param Visualization_Method: Method for visualization ('PCA', 'tsne', or 'UMAP').
    :return: 2D representation of the feature data.
    """
    if issparse(fea):
        fea = fea.toarray()
    if Visualization_Method.lower() == "pca":
        pca = PCA(n_components=2)
        Y = pca.fit_transorm(fea)
    elif Visualization_Method.lower() == "tsne":
        tsne = TSNE(n_components=2, random_state=random_seed)
        Y = tsne.fit_transform(fea)
    elif Visualization_Method.lower() == "umap":
        umap_model = UMAP(n_components=2, random_state=random_seed)
        Y = umap_model.fit_transform(fea)
    else:
        raise ValueError("Visualization_Method should be 'PCA', 'tsne', or 'UMAP'")
    return Y