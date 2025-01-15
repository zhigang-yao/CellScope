# Test Data

This folder contains datasets for testing **CellScope**, including directly available `.mat` files and `.h5ad` files, as well as external download links for additional datasets.

## Available Datasets

### Directly Available in This Repository:
| File Name     | Tissue |
|--------------|-------------|
| **GSE59739.mat** | **Mouse Lumbar** |
| **GSE83139.mat** | **Human Pancreatic** |
| **PHS000424V9P2_Brain.mat** | **Human Brain** |
| **GSE132042_Pancreas.mat** | **Mouse Pancreas* |
| **GSE67835.mat** | ** Human Brain** |
| **GSE103322.mat** | **Human oral cavity** |

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
