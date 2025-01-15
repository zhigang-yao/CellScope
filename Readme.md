# CellScope

CellScope is an innovative framework designed for constructing high-resolution, multi-scale cell atlases. It provides detailed information about gene expression, molecular characteristics, and functional states of individual cells, enabling researchers to understand cellular heterogeneity, development, and function in unprecedented detail. Traditional methods for developing comprehensive cell atlases often face challenges, such as reliance on conventional unsupervised learning algorithms and difficulties in identifying rare cell populations. CellScope addresses these limitations by employing a refined gene selection strategy, utilizing manifold fitting for cell clustering, and introducing a tree-structured visualization to represent cellular hierarchy intuitively.

## Workflow Overview (**Figure c**)

1. **Manifold Assumption**
    - CellScope assumes that true biological structures in single-cell data lie on a **low-dimensional manifold** that captures cellular states and subtypes.
    - However, observed gene expression is affected by two types of noise (**Figure a**):
        - **Housekeeping gene noise** (ubiquitous genes that do not distinguish cell types).
        - **Technical noise** (mRNA loss, inefficient capture, sequencing errors).

2. **Two-Step Manifold Fitting (**Figure b**)**
    - **Step 1: Noise Mitigation**
        - Separates **signal space** (informative genes) from **noise space** (housekeeping genes).
        - Identifies **manifold seeds** and **highly reliable cliques** to filter out irrelevant variation.
    - **Step 2: Cell-Type Stratification**
        - Assigns low-density cells (transitional states, noisy profiles) to the nearest **cell-type submanifold**.
        - Enhances biological signal while reducing artifacts.

3. **Neighborhood Graph & Clustering**
    - Constructs a **cell-to-cell similarity graph** based on gene expression.
    - Performs **agglomerative clustering** to define biologically meaningful cell types.

4. **Tree-Structured Visualization**
    - Integrates **UMAP** for intuitive low-dimensional visualization.
    - Generates a **hierarchical tree** to depict relationships between cell types.
    - Annotates key regulatory genes driving cell fate decisions (**Figure d**).

<p align="center">
  <img src="/Workflow.jpg" width="800">
</p>

---

## 📂 Test Data

To facilitate testing and benchmarking, we provide several datasets that can be directly accessed from this repository or downloaded from external sources.

### ✅ Directly Available in This Repository:
| Dataset Name       | Tissue            | #Cells  | #Genes  | #Cell Types |
|-------------------|------------------|---------|---------|------------|
| **GSE36552.mat**  | Mouse Lumbar      | 90      | 20,214  | 6          |
| **SRP041736.mat** | Human Pancreatic  | 249     | 14,805  | 11         |
| **GSE67835.mat**  | Human Brain       | 457     | 19,950  | 7          |
| **GSE58739.mat**  | Mouse Pancreas    | 466     | 22,088  | 9          |

### 🌍 Download from External Sources:
| Dataset Name           | Tissue                | #Cells  | #Genes  | #Cell Types | Download Link |
|------------------------|----------------------|---------|---------|------------|---------------|
| **SCR015820-MidBrain** | Human Brain          | 4,714   | 59,357  | 11         | [Download Here](https://datasets.cellxgene.cziscience.com/5488ff72-58ed-4f0d-913c-1b6d4d8412b1.h5ad) |
| **SCR015820-LNC**      | Human Brain          | 6,877   | 59,357  | 10         | [Download Here](https://datasets.cellxgene.cziscience.com/160cef00-39e7-49a3-a882-da7eb0e215fa.h5ad) |
| **PBMC**              | PBMC                  | 49,139  | 22,827  | 25         | [Download Here](https://datasets.cellxgene.cziscience.com/fbe23743-b3b5-4e2c-9bb2-95ee14d36783.h5ad) |
| **GSE178101**         | Human Fallopian Tubes | 77,536  | 36,960  | 12         | [Download Here](https://datasets.cellxgene.cziscience.com/26f36ff7-17b6-4285-8b35-9512dcae307b.h5ad) |

### 📖 Data Usage
To use these datasets in **CellScope**, refer to **`Test Code.ipynb`**, which provides a step-by-step tutorial on:
- **Loading `.mat` files** using `scipy.io.loadmat()`
- **Loading `.h5ad` files** using `scanpy.read_h5ad()`
- **Running CellScope analysis on these datasets**

## Documentation

For detailed tutorials, installation instructions, and usage examples, please visit the [CellScope Documentation](https://cellscope.readthedocs.io/en/latest/).


## Features
- **Refined Gene Selection**: Implements an advanced strategy for selecting relevant genes to improve the accuracy of cell clustering.
- **Manifold Fitting for Clustering**: Utilizes manifold fitting techniques to enhance the clustering performance of cells, allowing for better identification of cell populations.
- **Tree-Structured Visualization**: Introduces a novel visualization method that intuitively represents the hierarchical relationships between cells, facilitating easier interpretation of complex data.
- **Superior Performance**: Demonstrates superior performance across 36 single-cell atlas datasets, outshining state-of-the-art analysis pipelines such as Seurat and Scanpy.
- **Enhanced Biological Insights**: Capable of uncovering novel biological insights, including the identification of new cell subtypes with distinct functions.
- **Exceptional Computational Efficiency**: Optimized for scalability to large datasets with minimal dependence on hyper-parameter tuning.
- **Disease Analysis Pipeline**: Includes a specialized pipeline for analyzing disease-control single-cell atlases, revealing distinct cellular expression patterns in conditions such as COVID-19.

## Dependencies
- numpy
- pandas
- matplotlib
- seaborn
- scikit-learn
- scipy
- umap-learn
- joblib
- anndata
- plotly

## License
This project is licensed under the MIT License. See the LICENSE file for details.

## Citation

If you use CellScope in your research, please cite:

```bibtex
@article{Li2025CellScope
  title={CellScope: High-Performance Cell Atlas Workflow with Tree-Structured Representation},
  author={Bingjie Li, Runyu Lin, Tianhao Ni, Guanao Yan, Mannix Burns, Jingyi Jessica Li and Zhigang Yao},
  journal={},
  year={2025}
}
```

## Contact

Contact Email: zhigang.yao@nus.edu.sg, bjlistat@nus.edu.sg, and thni@zju.edu.cn  

Copyright (c) 2025 Yao group

