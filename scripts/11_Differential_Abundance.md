# Getting started with differential abundance analysis

This guide will walk you through setting up the environment and installing the necessary packages to perform differential abundance (DA) analysis.

## 1. Create the environment

First, create a new Conda environment named `milo` with the required R and Python versions:

```
conda create -n milo -c conda-forge r-base=4.4.0 python=3.9 
```

If you encounter the following error:

```
SafetyError: The package for r-base located at /opt/anaconda3/pkgs/r-base-4.4.0-hba25a81_1 ......
```

then remove the problematic package and redo the last step. 

```
rm -r /opt/anaconda3/pkgs/r-base-4.4.0-hba25a81_1
```

Otherwise, activate the environment:

```
conda activate milo
```

## 2. Install R packages

Before proceeding to install milopy, start R and install the necessary R packages:

```
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("edgeR") 
install.packages('statmod')
```

## 3. Install and set up milopy

Exit R and continue with the following steps:

```
git clone https://github.com/emdann/milopy.git
cd milopy
pip install .
```

## 4. Enable milo in JupyterLab

To use milo within JupyterLab, install the IPython kernel:

```
pip install ipykernel
python -m ipykernel install --user --name=milo
```