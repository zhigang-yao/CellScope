# Test Data

This folder contains datasets for testing **CellScope**, including directly available `.mat` files and `.h5ad` files, as well as external download links for additional datasets.

## Available Datasets

| File Name         | Tissue            | #Cells | #Genes  | #Cell Types |
|------------------|-----------------|--------|--------|------------|
| **GSE36552.mat**  | Mouse Lumbar     | 90     | 20,214 | 6          |
| **SRP041736.mat** | Human Pancreatic | 249    | 14,805 | 11         |
| **GSE67835.mat**  | Human Brain      | 457    | 19,950 | 7          |
| **GSE58739.mat**  | Mouse Pancreas   | 466    | 22,088 | 9          |


### Download from External Sources:
For larger datasets, please download them from the following links:

| Dataset Name | File Format | Download Link |
|-------------|------------|---------------|
| **Human Cell Atlas** | `.h5ad` | [Download Here](https://example.com/human_cell_atlas) |
| **Mouse Brain Atlas** | `.mat` | [Download Here](https://example.com/mouse_brain_atlas) |
| **COVID-19 Single-Cell Atlas** | `.h5ad` | [Download Here](https://example.com/covid19_atlas) |

## Usage Guide

To use these datasets in CellScope, refer to **`Test Code.ipynb`**, which provides a step-by-step tutorial on:

- Loading `.mat` files using `scipy.io.loadmat()`
- Loading `.h5ad` files using `scanpy.read_h5ad()`
- Running CellScope analysis on these datasets
