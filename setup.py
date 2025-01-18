from setuptools import setup, find_packages

setup(
    name='CellScope',
    version='0.1.0',
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
    install_requires=[
        'numpy>=1.24,<2.2',
        'pandas',
        'matplotlib',
        'seaborn',
        'scikit-learn',
        'scipy',
        'umap-learn',
        'joblib',
        'requests',
        'anndata',
        'plotly',
        'kaleido'
    ],
)
