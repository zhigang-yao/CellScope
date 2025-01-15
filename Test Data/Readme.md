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

| Dataset Name           | Tissue                | #Cells | #Genes | #Cell Types | Download Link |
|------------------------|----------------------|--------|--------|-------------|---------------|
| **SCR015820-MidBrain** | Human Brain         |  4,714  | 59,357 | 11          | [Download Here](https://datasets.cellxgene.cziscience.com/5488ff72-58ed-4f0d-913c-1b6d4d8412b1.h5ad) |
| **SCR015820-LNC**      | Human Brain         |  6,877  | 59,357 | 10          | [Download Here](https://datasets.cellxgene.cziscience.com/160cef00-39e7-49a3-a882-da7eb0e215fa.h5ad) |
| **PBMC**              | PBMC                 | 49,139  | 22,827 | 25          | [Download Here](https://datasets.cellxgene.cziscience.com/fbe23743-b3b5-4e2c-9bb2-95ee14d36783.h5ad) |
| **GSE178101**         | Human Fallopian Tubes | 77,536  | 36,960 | 12          | [Download Here](https://datasets.cellxgene.cziscience.com/26f36ff7-17b6-4285-8b35-9512dcae307b.h5ad) |


## Usage Guide

To use these datasets in CellScope, refer to **`Test Code.ipynb`**, which provides a step-by-step tutorial on:

- Loading `.mat` files using `scipy.io.loadmat()`
- Loading `.h5ad` files using `scanpy.read_h5ad()`
- Running CellScope analysis on these datasets
