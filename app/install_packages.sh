#!/bin/bash

# Note to advanced users, only exact Biopython and Scikit-Learn versions should be necessary
# All other packages should be compatible with multiple versions, but have not been tested

conda install -c conda-forge umap-learn=0.5.1
conda install -c conda-forge biopython=1.76
conda install -c conda-forge scipy=1.4.1
conda install pandas=1.0.3
conda install numpy=1.18.1
conda install matplotlib=3.1.3
conda install scikit-learn=0.22.1
conda install seaborn=0.10.1
# conda install -c conda-forge kivy=1.11.1
# Note, it seems that kivy 1.11.1 is buggy on newer macOS
# If you install kivy last, you can use kivy 2.0.0 instead
conda install -c conda-forge kivy=2.0.0
