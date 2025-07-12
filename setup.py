from setuptools import setup, find_packages

setup(
    name='CellScope-RNA',
    version='0.1.4',
    author='Tianhao Ni',
    author_email='thni@zju.edu.cn',
    description='A package for analyzing and visualizing gene expression data',
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.7, <3.11',
    install_requires = [
    'anndata>=0.11.4',
    'fastcluster>=1.2.6',
    'joblib>=1.4.2',
    'kaleido>=0.2.1',
    'matplotlib>=3.10.1',
    'nbformat>=5.10.4',
    'networkx>=3.4.2',
    'numba>=0.61.2',
    'numpy>=1.26.4,<2.2',  
    'pandas>=2.2.3',
    'plotly>=6.0.1',
    'requests>=2.32.3',
    'scikit-learn>=1.6.1',
    'scipy>=1.13.1',
    'seaborn>=0.13.2',
    'statsmodels>=0.14.4',
    'threadpoolctl>=3.6.0',
    'umap-learn>=0.5.7'
    ],
)
