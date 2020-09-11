# AIMS - An Automated Immune Molecule Separator

# Description
The primary goal of AIMS is to identify discriminating factors between two distinct sets of immune molecules. Currently the software can be applied to
immunoglobulin (Ig) molecules such as T cell receptors and antibodies, and major histocompatibility complex (MHC) and MHC-like molecules. 
AIMS is a python package distributed in both a notebook and GUI format. Those wishing to use the notebooks
just need to download this repository and the necessary packages identified in AIMS/app/install_packages.sh. Those
wishing to use the GUI, particularly those relatively new to programming, can follow the installation instructions below.
Example data is provided in AIMS/app/ab_testData and AIMS/app/mhc_testData, and an example of an application of AIMS can be seen in
this preprint: https://www.biorxiv.org/content/10.1101/2020.07.30.229013v1

When publishing analysis from this software, please cite the preprint:

Boughter CT, Borowska MT, Guthmiller JJ, Bendelac A, Wilson PC, Roux B, Adams EJ. Biochemical Patterns of Antibody Polyreactivity Revealed Through a Bioinformatics-
 Based Analysis of CDR Loops. BioRxiv. 2020.

And be on the lookout for the published version.

# Installation
For more advanced users familiar with python modules, skip to step 5 to find a complete list of python dependencies. For beginners, these instructions will help install all necessary packages and programs needed to run the GUI. Mac/Linux OS preferred, these steps outline installation via Anaconda. Other installations should be supported but have not been tested. 

1) Install Anaconda (https://www.anaconda.com/products/individual) to manage the python packages we're going to be using. This can be a fairly large package, so if space is at a premium for your computer, you can instead install miniconda (https://docs.conda.io/en/latest/miniconda.html). NOTE: Windows users should likely install the full Anaconda package, for a contained environment to run python programs from.

2) Test that your conda install is working properly by creating a conda environment. Windows OS users, you will likely do this within the Anaconda application. Mac/Linux users, open the terminal application. Once terminal is open, type:

conda create -n aims-env python=3.7

If anaconda/miniconda is installed properly, a Y/N prompt should appear. Type "y" then hit the "enter key" and you will create a conda environment.

3) If you haven't already, you should download the code from this repository.

4) Next, navigate to the new folder created from this repository. First, open up terminal and enter the environment created earlier by typing:

conda activate aims-env

Again, windows users will have to activate the environment from within Anaconda.
You should now see a little extra bit of text on your terminal command line that looks something like "(aims-env)". If this didn't work for some reason, an error message should pop up, otherwise assume you're fine.

Use terminal to navigate into the newly downloaded folder. If you've never used terminal before, you can type in "cd" and then drag and drop the folder into the terminal. Doing so should automatically populate the "path" to the folder. Then hit enter.

When I do this, my terminal line reads "cd /Users/boughter/Desktop/AIMS" hopefully you see something similar, and when you hit enter you can move to the next step.

5) In terminal, type:

./install_packages.sh

and hit enter. It should just work, and you'll be prompted with a bunch of [y]/n questions, for which you should consistently hit "y" then the "enter key" for.

If this doesn't work, and you get some kind of an error instead of the prompts, type each of these lines (or copy/paste) one by one, hitting enter after each one:

conda install -c conda-forge kivy
conda install -c conda-forge biopython
conda install -c conda-forge scipy
conda install pandas
conda install numpy
conda install matplotlib
conda install scikit-learn

6) Everything should now be installed, you should now be able to open up the software! Navigate to the app in terminal by typing: 

cd app

and then finally enter the command

python AIMS.py

From there, the GUI should open. A step by step instruction guide for GUI usage will be shortly uploaded.
