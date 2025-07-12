import matplotlib.pyplot as plt
import numpy as np
import umap
from sklearn.metrics import adjusted_rand_score
import os
from scipy.sparse import issparse
import networkx as nx
from matplotlib.patches import Circle
from matplotlib.animation import FuncAnimation
from IPython.display import HTML

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

def add_hyphen_between_digits(name):
    result = []
    for char in name:
        if char.isdigit() and result and result[-1].isdigit():
            result.append('-')  
        result.append(char)
    return ''.join(result)

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
        new_Title_all = [add_hyphen_between_digits(name) for name in Title_all]
        Title_all = []
        for ii in range(len(new_Title_all)):
            Title_all.append(new_Title_all[ii])
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

def tree_structure_visualization_static(T,step0,step1,Title_1,Title_all,Y_initial,Y_1,Y_all,index_1,index_all, save_figure=False, filename="tree_visualization_static.png"):
    celltype = T[:,step1]
    unique_celltypes = np.unique(celltype)
    colors = np.random.rand(len(unique_celltypes), 3)
    colors = np.array([
        [0.902, 0.098, 0.110],   # Red
        [0.235, 0.706, 0.294],   # Green
        [0.000, 0.459, 0.757],   # Blue
        [0.961, 0.510, 0.118],   # Orange
        [0.569, 0.118, 0.706],   # Purple
        [0.220, 0.557, 0.235],   # Forest Green
        [0.808, 0.361, 0.000],   # Brown
        [0.580, 0.404, 0.741],   # Violet
        [0.549, 0.337, 0.294],   # Dark Brown
        [0.890, 0.467, 0.761],   # Pink
        [0.498, 0.498, 0.000],   # Olive
        [0.737, 0.741, 0.133],   # Lime Green
        [0.055, 0.647, 0.647],   # Teal
        [0.941, 0.502, 0.502],   # Light Pink
        [0.000, 0.502, 0.502],   # Dark Teal
        [0.686, 0.933, 0.933],   # Light Blue
        [0.502, 0.000, 0.502],   # Dark Purple
        [0.933, 0.867, 0.510],   # Light Yellow
        [0.400, 0.400, 0.400],   # Medium Gray
        [0.000, 0.000, 0.000]    # Black
    ])
    color_map = {label: colors[i] for i, label in enumerate(unique_celltypes)}
    colors_array = np.array([color_map[label] for label in celltype])

    short_name_indices = []

    for ii in range(len(Title_all)):
        name = Title_all[ii]
        name = name[7:]  
        if len(name) <= 5: 
            short_name_indices.append(ii) 
    Title_all = [Title_all[ii] for ii in short_name_indices]
    Y_all = [Y_all[ii] for ii in short_name_indices]
    index_all = [index_all[ii] for ii in short_name_indices]

    label0 = T[:,step0]
    n = len(np.unique(label0))+len(Title_all)+1
    n = int(n)
    Name_all = ['Initial']
    for jj in range(len(np.unique(label0))):
        Name_all.append(str(jj))
    for jj in range(len(Title_all)):
        Name_all.append(Title_all[jj][7:])
    A = np.zeros((n,n))
    A[0,1:len(np.unique(label0))+1] = 1
    for ii in range(len(Title_all)):
        name = Title_all[ii]
        name = name[7:]
        temp = Name_all.index(name[:-2])
        A[ii+1+len(np.unique(label0)),temp] = 1
    A = np.maximum(A,A.T)

    G = nx.from_numpy_array(A)
    root = 0

    children = {n: list(G.neighbors(n)) for n in G.nodes()}
    def count_leaves(u, p):
        ch = [v for v in children[u] if v != p]
        if not ch:
            return 1
        return sum(count_leaves(v, u) for v in ch)

    size = {u: count_leaves(u, None) for u in G.nodes()}

    pos = {}
    R = 2.5
    def assign(u, p, theta0, theta1, depth):
        theta_mid = (theta0 + theta1) / 2
        pos[u] = (depth * R * np.cos(theta_mid), depth * R * np.sin(theta_mid))
        ch = [v for v in children[u] if v != p]
        total = sum(size[v] for v in ch)
        theta = theta0
        for v in ch:
            span = (theta1 - theta0) * size[v] / total
            assign(v, u, theta, theta + span, depth + 1)
            theta += span

    assign(root, None, 0, 2*np.pi, 0)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(30, 15))

    node_radius = 0.8
    for u, v in G.edges():
        if u==0:
            node_radius = 1.1
        else:
            node_radius = 1.1
        x0, y0 = pos[u]
        x1, y1 = pos[v]
        dx, dy = x1 - x0, y1 - y0
        dist = np.sqrt(dx*dx + dy*dy)
        x0e = x0 + dx/dist * node_radius
        y0e = y0 + dy/dist * node_radius
        x1e = x1 - dx/dist * node_radius
        y1e = y1 - dy/dist * node_radius
        ax1.plot([x0e, x1e], [y0e, y1e], '-', lw=1, color='gray')
        ax2.plot([x0e, x1e], [y0e, y1e], '-', lw=1, color='gray')

    for ii in range(len(pos)):
        if ii==0:
            temp = 0
            Y = Y_initial
            index = range(T.shape[0])
            node_radius = 1.1
        elif ii<len(Title_1)+1:
            temp = 1
            Y = Y_1[ii-1]
            index = index_1[ii-1]
            node_radius = 1.1
        else:
            Y = Y_all[ii-len(Title_1)-1]
            temp = 2
            index = index_all[ii-len(Title_1)-1]
            node_radius =1.1

        [x,y] = pos[ii]
        circ = Circle((x, y), radius=node_radius, facecolor='none', edgecolor='black', lw=1.2)
        ax1.add_patch(circ)
        circ = Circle((x, y), radius=node_radius, facecolor='none', edgecolor='black', lw=1.2)
        ax2.add_patch(circ)

        
        main_center = np.mean(Y, axis=0)
        Y_initial_centered = Y - main_center
        main_scale = 0.85 *node_radius / np.max(np.linalg.norm(Y_initial_centered, axis=1))
        n_points = Y.shape[0]
        point_size = 10 if n_points < 200 else (5 if n_points < 500 else 1)
        if len(np.unique(colors_array[index], axis=0)) == 1 and temp==1:
            ax2.scatter(*(Y_initial_centered * main_scale + np.array(pos[ii])).T,s = 8,c=colors_array[index],label=f'Cluster {ii-1}')
        elif len(np.unique(colors_array[index], axis=0)) == 1 and temp==2:
            ax2.scatter(*(Y_initial_centered * main_scale + np.array(pos[ii])).T,s = 8,c=colors_array[index],label=Title_all[ii-len(Title_1)-1])
        else:
            ax2.scatter(*(Y_initial_centered * main_scale + np.array(pos[ii])).T,s = 8,c=colors_array[index])
        if ii == 0:
            title = "Initial" 
        elif ii < len(Title_1) + 1:
            title = f"Cluster {ii-1}" 
        else:
            title = Title_all[ii - len(Title_1) - 1]  

        if ii<len(Title_1)+1:
            ax1.text(
                x, y, title,
                ha='center', va='center',
                fontsize=12,  
                color='black'
            )
        else:
            ax1.text(
                x, y, title,
                ha='center', va='center',
                fontsize=12,  
                color='black'
            )
        
    ax1.set_aspect('equal')
    ax1.axis('off')
    ax2.set_aspect('equal')
    ax2.axis('off')
    

    if save_figure:
        # Validate file format
        supported_formats = {'.png', '.jpg', '.jpeg', '.pdf'}
        file_ext = filename[filename.rfind('.'):].lower()
        
        if file_ext not in supported_formats:
            raise ValueError(f"Unsupported file format: {file_ext}. Use one of {supported_formats}")
        
        # Create directory if needed
        os.makedirs(os.path.dirname(filename), exist_ok=True)

        if file_ext == '.pdf':
            plt.savefig(filename, 
                    format='pdf', 
                    bbox_inches='tight', 
                    dpi=1200,  # High-res vector format
                    metadata={'CreationDate': None})  # Better reproducibility

        elif file_ext in ('.jpg', '.jpeg'):
            plt.savefig(filename, 
                    format='jpeg', 
                    bbox_inches='tight', 
                    quality=100,  # Maximum quality
                    dpi=600,     # High DPI for raster
                    optimize=True, 
                    progressive=True)  # Enhanced web display

        elif file_ext == '.tif':
            plt.savefig(filename,
                    format='tiff',
                    bbox_inches='tight',
                    dpi=1200,
                    compression='lzw')  # Lossless compression

        else:  # .png
            plt.savefig(filename, 
                    format='png', 
                    bbox_inches='tight', 
                    dpi=600,            # High DPI
                    transparent=False,
                    pil_kwargs={'compress_level': 1})  # Minimal PNG compression

    plt.show()
    plt.close() 

def tree_structure_visualization_dynamic(T,step0,step1,Title_1,Title_all,Y_initial,Y_1,Y_all,index_1,index_all,interval=800,save_gif=False, 
                                       gif_filename="tree_animation.gif", dpi=100):
    """分层动画的树形结构可视化"""
    
    # ===================== 数据准备阶段 =====================
    celltype = T[:, step1]
    unique_celltypes = np.unique(celltype)
    colors = np.array([
        [0.902, 0.098, 0.110], [0.235, 0.706, 0.294], [0.000, 0.459, 0.757],
        [0.961, 0.510, 0.118], [0.569, 0.118, 0.706], [0.220, 0.557, 0.235],
        [0.808, 0.361, 0.000], [0.580, 0.404, 0.741], [0.549, 0.337, 0.294],
        [0.890, 0.467, 0.761], [0.498, 0.498, 0.000], [0.737, 0.741, 0.133],
        [0.055, 0.647, 0.647], [0.941, 0.502, 0.502], [0.000, 0.502, 0.502],
        [0.686, 0.933, 0.933], [0.502, 0.000, 0.502], [0.933, 0.867, 0.510],
        [0.400, 0.400, 0.400], [0.000, 0.000, 0.000]
    ])
    color_map = {label: colors[i % len(colors)] for i, label in enumerate(unique_celltypes)}
    colors_array = np.array([color_map[label] for label in celltype])

    # 标题处理
    short_name_indices = [ii for ii, name in enumerate(Title_all) if len(name[7:]) <= 5]
    Title_all = [Title_all[ii] for ii in short_name_indices]
    Y_all = [Y_all[ii] for ii in short_name_indices]
    index_all = [index_all[ii] for ii in short_name_indices]

    # 构建邻接矩阵
    label0 = T[:, step0]
    n = len(np.unique(label0)) + len(Title_all) + 1
    Name_all = ['Initial'] + [str(jj) for jj in range(len(np.unique(label0)))] + [t[7:] for t in Title_all]
    
    A = np.zeros((n, n))
    A[0, 1:len(np.unique(label0)) + 1] = 1
    for ii in range(len(Title_all)):
        name = Title_all[ii][7:-2]
        temp = Name_all.index(name)
        A[ii + 1 + len(np.unique(label0)), temp] = 1
    A = np.maximum(A, A.T)

    # ===================== 图形结构计算 =====================
    G = nx.from_numpy_array(A)
    root = 0

    # 计算节点位置
    children = {n: list(G.neighbors(n)) for n in G.nodes()}
    
    def count_leaves(u, p=None):
        ch = [v for v in children[u] if v != p]
        return 1 if not ch else sum(count_leaves(v, u) for v in ch)
    
    size = {u: count_leaves(u) for u in G.nodes()}
    
    pos = {}
    R = 2.5
    def assign(u, p, theta0, theta1, depth):
        theta_mid = (theta0 + theta1) / 2
        pos[u] = (depth * R * np.cos(theta_mid), depth * R * np.sin(theta_mid))
        ch = [v for v in children[u] if v != p]
        total = sum(size[v] for v in ch)
        theta = theta0
        for v in ch:
            span = (theta1 - theta0) * size[v] / total
            assign(v, u, theta, theta + span, depth + 1)
            theta += span
            
    assign(root, None, 0, 2 * np.pi, 0)

    # ===================== 动画初始化 =====================
    plt.ioff()
    
    # 分层节点顺序
    nodes_order = [root]  # 根节点
    middle_nodes = []      # 中间层节点
    leaf_nodes = []       # 叶节点
    
    # 广度优先遍历分类节点
    current_level = [root]
    while current_level:
        next_level = []
        for node in current_level:
            neighbors = sorted([n for n in G.neighbors(node) if n not in nodes_order])
            for n in neighbors:
                if n < len(Title_1) + 1 and n != root:
                    middle_nodes.append(n)
                else:
                    leaf_nodes.append(n)
            next_level.extend(neighbors)
            nodes_order.extend(neighbors)
        current_level = next_level
    
    paired_leaf_nodes = [leaf_nodes[i:i+2] for i in range(0, len(leaf_nodes), 2)]
    middle_nodes = [n for n in G.nodes() if 0 < n <= len(Title_1)]
    nodes_order = [root] + [middle_nodes] + paired_leaf_nodes

    all_pos = np.array(list(pos.values()))
    padding = 1.2  # 边距系数
    xmin, xmax = all_pos[:,0].min()-padding, all_pos[:,0].max()+padding
    ymin, ymax = all_pos[:,1].min()-padding, all_pos[:,1].max()+padding

    # ===================== 动画初始化 =====================
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(24, 12))
    for ax in (ax1, ax2):
        ax.set_xlim(xmin, xmax)
        ax.set_ylim(ymin, ymax)

        ax.set_aspect('equal')
        ax.axis('off')

    # 状态管理系统
    state = {
        'plotted_edges': set(),
        'plotted_nodes': set(),
        'all_scatters': [],
        'all_texts': [],
        'all_circles': []
    }

    # ===================== 动画更新函数 =====================
    def update(frame):
        current_element = nodes_order[frame]
        
        # 根节点处理
        if frame == 0:
            current_node = current_element
            # 绘制根节点
            node_radius = 1.1
            circ = Circle(pos[current_node], radius=node_radius, 
                        facecolor='none', edgecolor='black', lw=1.2)
            ax1.add_patch(circ)
            ax2.add_patch(Circle(pos[current_node], radius=node_radius, 
                               facecolor='none', edgecolor='black', lw=1.2))
            state['plotted_nodes'].add(current_node)
            
            # 根节点文字
            text = ax1.text(*pos[current_node], "Initial", 
                          ha='center', va='center', fontsize=14, color='black')
            state['all_texts'].append(text)
            
            # 根节点散点
            Y = Y_initial
            index = range(T.shape[0])
            main_center = np.mean(Y, axis=0)
            Y_centered = Y - main_center
            scale = 0.85 * 1.1 / np.max(np.linalg.norm(Y_centered, axis=1))
            scatter = ax2.scatter(*(Y_centered * scale + pos[current_node]).T,
                                s=10, c=colors_array[index], alpha=0.8)
            state['all_scatters'].append(scatter)
        
        # 中间层批量处理
        elif frame == 1:
            for current_node in current_element:
                # 绘制边
                for u, v in G.edges(current_node):
                    edge = frozenset({u, v})
                    if edge not in state['plotted_edges']:
                        x0, y0 = pos[u]
                        x1, y1 = pos[v]
                        dx, dy = x1 - x0, y1 - y0
                        dist = np.hypot(dx, dy)
                        x0e = x0 + dx/dist * 1.1
                        y0e = y0 + dy/dist * 1.1
                        x1e = x1 - dx/dist * 1.1
                        y1e = y1 - dy/dist * 1.1
                        
                        line1, = ax1.plot([x0e, x1e], [y0e, y1e], '-', lw=1, color='gray')
                        line2, = ax2.plot([x0e, x1e], [y0e, y1e], '-', lw=1, color='gray')
                        state['plotted_edges'].add(edge)
                
                # 绘制节点
                node_radius = 1.1
                circ = Circle(pos[current_node], radius=node_radius, 
                            facecolor='none', edgecolor='black', lw=1.2)
                ax1.add_patch(circ)
                ax2.add_patch(Circle(pos[current_node], radius=node_radius, 
                                   facecolor='none', edgecolor='black', lw=1.2))
                state['plotted_nodes'].add(current_node)
                
                # 中间层文字
                text = ax1.text(*pos[current_node], f"Cluster {current_node-1}", 
                              ha='center', va='center', fontsize=12, color='black')
                state['all_texts'].append(text)
                
                # 中间层散点
                Y = Y_1[current_node-1]
                index = index_1[current_node-1]
                main_center = np.mean(Y, axis=0)
                Y_centered = Y - main_center
                scale = 0.85 * 1.1 / np.max(np.linalg.norm(Y_centered, axis=1))
                scatter = ax2.scatter(*(Y_centered * scale + pos[current_node]).T,
                                    s=8, c=colors_array[index], alpha=0.7)
                state['all_scatters'].append(scatter)
        
        else:
            for current_node in current_element:
                # 绘制边（从父节点到当前叶节点）
                for u, v in G.edges(current_node):
                    edge = frozenset({u, v})
                    if edge not in state['plotted_edges']:
                            # 计算边坐标
                            x0, y0 = pos[u]
                            x1, y1 = pos[v]
                            dx, dy = x1 - x0, y1 - y0
                            dist = np.hypot(dx, dy)
                            x0e = x0 + dx/dist * 1.1
                            y0e = y0 + dy/dist * 1.1
                            x1e = x1 - dx/dist * 1.1
                            y1e = y1 - dy/dist * 1.1
                            
                            # 绘制边
                            line1, = ax1.plot([x0e, x1e], [y0e, y1e], '-', lw=1, color='gray')
                            line2, = ax2.plot([x0e, x1e], [y0e, y1e], '-', lw=1, color='gray')
                            state['plotted_edges'].add(edge)
                
                # 绘制叶节点（如果未绘制过）
                if current_node not in state['plotted_nodes']:
                    # 绘制节点圆圈
                    node_radius = 1.1
                    circ = Circle(pos[current_node], radius=node_radius, 
                                facecolor='none', edgecolor='black', lw=1.2)
                    ax1.add_patch(circ)
                    ax2.add_patch(Circle(pos[current_node], radius=node_radius, 
                                       facecolor='none', edgecolor='black', lw=1.2))
                    state['plotted_nodes'].add(current_node)
                    
                    # 添加文字标签
                    title = Title_all[current_node - len(Title_1) - 1][7:]
                    text = ax1.text(*pos[current_node], title, 
                                  ha='center', va='center', fontsize=12, color='black')
                    state['all_texts'].append(text)
                    
                    # 添加散点图
                    Y = Y_all[current_node - len(Title_1) - 1]
                    index = index_all[current_node - len(Title_1) - 1]
                    main_center = np.mean(Y, axis=0)
                    Y_centered = Y - main_center
                    scale = 0.85 * 1.1 / np.max(np.linalg.norm(Y_centered, axis=1))
                    scatter = ax2.scatter(*(Y_centered * scale + pos[current_node]).T,
                                        s=8, c=colors_array[index], alpha=0.7)
                    state['all_scatters'].append(scatter)

        return state['all_scatters'] + state['all_texts'] + state['all_circles']

    ani = FuncAnimation(
        fig,
        update,
        frames=len(nodes_order),
        interval=interval,
        blit=False,
        repeat=False,
        cache_frame_data=False
    )
    if save_gif:
        ani.save(gif_filename, writer='pillow', dpi=dpi, 
                fps=1000/interval,  # Calculate FPS based on interval
                progress_callback=lambda i, n: print(f'Saving progress: {i+1}/{n} frames', end='\r'))

    plt.close(fig)
    return HTML(ani.to_jshtml())

def tree_structure_visualization_step(T,step0,step1,Title1,Title_all,Y_initial,Y_1,Y_all,index_1,index_all,
    main_radius=2.0,
    orbit_radius=1.5,
    orbit_distance=4.0,
    branch_params=None,
    save_fig = False,
    save_path = ''
):
    if branch_params is None:
        branch_params = {
            'length': lambda d: 3.0/(1 + d), 
            'radius': lambda d: 1.5/(1 + d),
            'angle_offset': lambda d: np.pi/(4 + d)
        }
    Title_all_new = []
    for ii in range(len(Title_all)):
        Title_all_new.append(Title_all[ii][7:])
    Title_all = Title_all_new
    celltype = T[:,step1]
    unique_celltypes = np.unique(celltype)
    colors = np.random.rand(len(unique_celltypes), 3)
    color_map = {label: colors[i] for i, label in enumerate(unique_celltypes)}
    colors_array = np.array([color_map[label] for label in celltype])

    fig, ax = plt.subplots(figsize=(5, 5))
    ax.set_aspect('equal')
    ax.axis('off')
    
    main_center = np.mean(Y_initial, axis=0)
    Y_initial_centered = Y_initial - main_center
    main_scale = 0.95 * main_radius / np.max(np.linalg.norm(Y_initial_centered, axis=1))
    
    ax.scatter(*(Y_initial_centered * main_scale + main_center).T, c=colors_array, s = 1)
    ax.add_patch(plt.Circle(main_center, main_radius, edgecolor='black', fill=False, lw=3, zorder=11))
    
    angles = np.linspace(0, 2*np.pi, len(Y_1), endpoint=False)
    R_dict = {}
    C_dict = {}
    A_dict = {}
    Data_dict = {}
    for i, angle in enumerate(angles):
        orbit_center = (
            main_center[0] + orbit_distance * np.cos(angle),
            main_center[1] + orbit_distance * np.sin(angle)
        )
        if Y_1[i] is not None:

            orbit_data = Y_1[i] - np.mean(Y_1[i], axis=0)
            orbit_scale = orbit_radius / np.max(np.linalg.norm(orbit_data, axis=1))
            n_points = len(Y_1[i])
            point_size = 10 if n_points < 200 else (5 if n_points < 500 else 1)
            ax.scatter(*(orbit_data * orbit_scale + orbit_center).T, c=colors_array[index_1[i]], s=point_size)
            ax.add_patch(plt.Circle(orbit_center, orbit_radius, 
                                  edgecolor='black', fill=False, linewidth=2, zorder=9))
            ax.plot([main_center[0] + main_radius * np.cos(angle), 
                    orbit_center[0] - orbit_radius * np.cos(angle)],
                   [main_center[1] + main_radius * np.sin(angle), 
                    orbit_center[1] - orbit_radius * np.sin(angle)],
                   'black', linewidth=2, zorder=8)
            R_dict[str(i)] = orbit_radius
            C_dict[str(i)] = orbit_center
            A_dict[str(i)] = angle
            Data_dict[str(i)] = orbit_data * orbit_scale + orbit_center
    if save_fig:
        plt.savefig(save_path+'Initial.pdf', format='pdf', bbox_inches='tight', dpi=300)
    plt.show()

    for i, angle in enumerate(angles):
        
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.set_aspect('equal')
        ax.axis('off')
        
        orbit_center = (
            main_center[0] + orbit_distance * np.cos(angle),
            main_center[1] + orbit_distance * np.sin(angle)
        )
        orbit_data = Y_1[i] - np.mean(Y_1[i], axis=0)
        orbit_scale = orbit_radius / np.max(np.linalg.norm(orbit_data, axis=1))
        n_points = len(Y_1[i])
        point_size = 10 if n_points < 200 else (5 if n_points < 500 else 1)
        ax.scatter(*(orbit_data * orbit_scale + orbit_center).T, c=colors_array[index_1[i]], s=point_size)
        ax.add_patch(plt.Circle(orbit_center, orbit_radius, 
                                edgecolor='black', fill=False, linewidth=2, zorder=9))
        if save_fig:
            plt.savefig(save_path+'Cluster'+str(i)+'.pdf', format='pdf', bbox_inches='tight', dpi=300)
        print('Cluster'+str(i))
        plt.show()
    def tail_circle(
        orbit_center,      # 轨道圆圆心 (x, y)
        orbit_radius,       # 轨道圆半径
        orbit_data,
        tail_data,          # 尾部数据 [tail1_data, tail2_data])
        tail_length=1.5,    # 尾部连接线长度
        tail_radius=0.8,    # 尾部圆半径
        angle_offset=np.pi/4, # 分叉角度偏移
        base_angle=0.0,
        colors_orbit= None,
        colors_tail = None,
        orbit_title = None
    ):  
        C = []
        R = []
        Data = []
        angles = [base_angle + angle_offset, base_angle - angle_offset]
        for i, angle in enumerate(angles[:2]):
            colors_i = colors_tail[i]
            if i >= len(tail_data) or tail_data[i] is None:
                continue
            tail_center = (
                orbit_center[0] + (orbit_radius + tail_length) * np.cos(angle),
                orbit_center[1] + (orbit_radius + tail_length) * np.sin(angle)
            )

            ax.plot(
                [orbit_center[0],
                    tail_center[0] - tail_radius * np.cos(angle)],
                [orbit_center[1] ,
                    tail_center[1] - tail_radius * np.sin(angle)],
                color='black', linewidth=1.5
            )
            ax.add_patch(plt.Circle(tail_center, tail_radius,
                                    edgecolor='black', fill=False, linewidth=1.5))
            R.append(tail_radius)
            C.append([orbit_center[0] + (orbit_radius + tail_length) * np.cos(angle),
                orbit_center[1] + (orbit_radius + tail_length) * np.sin(angle)])
            
            if isinstance(tail_data[i], np.ndarray) and len(tail_data[i]) > 0:
                data = tail_data[i] - np.mean(tail_data[i], axis=0)
                scale = 0.95 * tail_radius / np.max(np.linalg.norm(data, axis=1))
                scaled_data = data * scale + tail_center
                n_points = len(scaled_data)
                point_size = 10 if n_points < 200 else (5 if n_points < 500 else 1)
                ax.scatter(*scaled_data.T, c=colors_i, s=point_size)
            Data.append(scaled_data)
        return R,C,Data
    for ii in range(int(len(Title_all)/2)):

        title = Title_all[2*ii]
        
        tail_data = [Y_all[2*ii],Y_all[2*ii+1]]
        orbit_index = title[:-2]
        orbit_center = C_dict[orbit_index]
        orbit_radius = R_dict[orbit_index]
        orbit_data = Data_dict[orbit_index]

        angle_offset = np.pi/6
        temp = np.union1d(index_all[2*ii],index_all[2*ii+1])
        fig, ax = plt.subplots(figsize=(5, 5))
        ax.set_aspect('equal')
        ax.axis('off')
        R,C,Data = tail_circle(orbit_center,orbit_radius,orbit_data,tail_data,base_angle=A_dict[title[0]],colors_orbit = colors_array[temp],colors_tail=[colors_array[index_all[2*ii],:],colors_array[index_all[2*ii+1],:]],angle_offset=angle_offset,orbit_title = orbit_index)
        R_dict[Title_all[2*ii]] = R[0]
        R_dict[Title_all[2*ii+1]] = R[1]
        C_dict[Title_all[2*ii]] = C[0]
        C_dict[Title_all[2*ii+1]] = C[1]
        Data_dict[Title_all[2*ii]] = Data[0]
        Data_dict[Title_all[2*ii+1]] = Data[1]
        if save_fig:
            plt.savefig(save_path+'Cluster'+title+'&2.pdf', format='pdf', bbox_inches='tight', dpi=300)
        print('Cluster'+title+'&2')
        plt.show()
    return R_dict,C_dict