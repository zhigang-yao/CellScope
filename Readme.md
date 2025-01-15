# CellScope

## Overview
CellScope is an innovative framework designed for constructing high-resolution, multi-scale cell atlases. It provides detailed information about gene expression, molecular characteristics, and functional states of individual cells, enabling researchers to understand cellular heterogeneity, development, and function in unprecedented detail. Traditional methods for developing comprehensive cell atlases often face challenges, such as reliance on conventional unsupervised learning algorithms and difficulties in identifying rare cell populations. CellScope addresses these limitations by employing a refined gene selection strategy, utilizing manifold fitting for cell clustering, and introducing a tree-structured visualization to represent cellular hierarchy intuitively.

<p align="center">
  <img src="/Workflow.jpg" width="800">
</p>

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

## License
This project is licensed under the MIT License. See the LICENSE file for details.

## Contact

Contact Email: zhigang.yao@nus.edu.sg, bjlistat@nus.edu.sg, and thni@zju.edu.cn  

Copyright (c) 2025 Yao group

