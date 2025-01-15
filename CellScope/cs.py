import numpy as np
from scipy.cluster.hierarchy import linkage, fcluster
import umap.umap_ as umap
from sklearn.decomposition import PCA, TruncatedSVD
from scipy.stats import f_oneway
from scipy.sparse.csgraph import connected_components
from joblib import Parallel, delayed
from sklearn.neighbors import NearestNeighbors
from scipy.sparse import issparse,csr_matrix,isspmatrix_csr
import random
from scipy.spatial.distance import squareform
from sklearn.manifold import TSNE
from umap.umap_ import UMAP
import warnings
from scipy.stats import ConstantInputWarning
warnings.filterwarnings("ignore", category=ConstantInputWarning)

np.random.seed(83)


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
        pca = PCA(n_components=num_pca)
        svd = TruncatedSVD(n_components=num_pca,random_state=random_seed)
        fea_pca = svd.fit_transform(fea)
    else:
        pca = PCA(n_components=num_pca,random_state=random_seed)
        fea_pca = pca.fit_transform(fea)
    nbrs = NearestNeighbors(n_neighbors=min(num_cell, 1000), algorithm='auto').fit(fea_pca)

    # Find k nearest neighbors
    distances, indices = nbrs.kneighbors(fea_pca)
    D_NB = distances[:, :knn]
    ID_NB = indices[:, :knn]
    # Select the top k nearest neighbor distances and indices
    rho = 1. / np.sum(D_NB, axis=1)
    delta = np.zeros(num_cell)
    for ii in range(num_cell):
        temp = np.where(rho > rho[ii])[0]
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
    knn = np.ceil(np.log2(fea_selected.shape[0]))

    if metric==None:
        if fea_selected.shape[0]<10000:
            metric = 'euclidean'
        else:
            metric = 'jaccard'
    else:
        metric = metric
    neighbors = umap.fuzzy_simplicial_set(
    fea_selected,
    n_neighbors=knn, 
    random_state=random_seed,
    metric=metric
)
    W = neighbors[0]
    num_cell = fea_selected.shape[0]
    if num_cell>num_cell_thre:
        if index==[]:
            index = random.sample(range(num_cell), num_cell_thre)
        index_unselected = set(range(num_cell))-set(index)
        index_unselected = list(index_unselected)
        knn_index = np.ceil(np.log2(len(index)))
        fea_selected_index = fea_selected[index,:]
        neighbors_index = umap.fuzzy_simplicial_set(
        fea_selected_index,
        n_neighbors=knn_index, 
        random_state=random_seed,
        metric=metric
    )   
        W1 = neighbors_index[0].toarray()
        D1 = 1-W1
        D1 = D1-np.diag(np.diag(D1))
        D1 = squareform(D1)
        Z = linkage(D1, method='average')
        T_all = np.zeros((num_cell,49))

        for ii in range(2,51):
            T = 100*np.ones((num_cell))
            T1 = fcluster(Z, ii, criterion='maxclust')
            T1 = T1-1
            T[index] = T1
            num_classes = len(np.unique(T1))
            num_remaining_samples = len(index_unselected)
            class_avg_similarities = np.zeros((num_remaining_samples, num_classes))
            for class_label in range(num_classes):
                class_indices = np.where(T1 == class_label)[0]
                class_indices = [index[i] for i in class_indices]
                class_similarity_matrix = W[np.ix_(index_unselected, class_indices)]
                class_avg_similarity = np.mean(class_similarity_matrix, axis=1)
                class_avg_similarities[:, class_label] = class_avg_similarity.reshape(-1)
            remaining_samples_labels = np.argmax(class_avg_similarities, axis=1)
            T[index_unselected] = remaining_samples_labels
            T_all[:,ii-2] = T
    else:
        W = W.toarray()
        D = 1-W
        D = D-np.diag(np.diag(D))
        D = squareform(D)
        Z = linkage(D, method='average')
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
        Y = pca.fit_transform(fea)
    elif Visualization_Method.lower() == "tsne":
        tsne = TSNE(n_components=2, random_state=random_seed)
        Y = tsne.fit_transform(fea)
    elif Visualization_Method.lower() == "umap":
        umap_model = UMAP(n_components=2, random_state=random_seed)
        Y = umap_model.fit_transform(fea)
    else:
        raise ValueError("Visualization_Method should be 'PCA', 'tsne', or 'UMAP'")
    return Y