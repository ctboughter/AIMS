#!/bin/bash

# Note to advanced users, only exact Biopython and Scikit-Learn versions should be necessary
# All other packages should be compatible with multiple versions, but have not been tested
# Note, these new package versions are compatible with python3.9

conda install -c conda-forge umap-learn=0.5.3
conda install -c conda-forge biopython=1.79
conda install -c conda-forge scipy=1.4.1
conda install pandas=1.5.3
conda install numpy=1.24.1
conda install matplotlib=3.7.1
conda install scikit-learn=1.3.0
conda install seaborn=0.12.2

# conda install -c conda-forge kivy=1.11.1
# Note, it seems that kivy 1.11.1 is buggy on newer macOS
# If you install kivy last, you can use kivy 2.1.0 instead
conda install -c conda-forge kivy=2.1.0
