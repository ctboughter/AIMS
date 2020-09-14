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

```
 conda create -n aims-env python=3.7
```

 If anaconda/miniconda is installed properly, a Y/N prompt should appear. Type "y" then hit the "enter key" and you will create a conda environment.

3) If you haven't already, you should download the code from this repository.

4) Next, navigate to the new folder created from this repository. First, open up terminal and enter the environment created earlier by typing:

```
 conda activate aims-env
```

Again, windows users will have to activate the environment from within Anaconda. You should now see a little extra bit of text on your terminal command line that looks something like "(aims-env)". If this didn't work for some reason, an error message should pop up, otherwise assume you're fine.

Use terminal to navigate into the newly downloaded folder. If you've never used terminal before, you can type in "cd" and then drag and drop the folder into the terminal. Doing so should automatically populate the "path" to the folder. Then hit enter.

 When I do this, my terminal line reads "cd /Users/boughter/Desktop/AIMS" hopefully you see something similar, and when you hit enter you can move to the next step.

5) In terminal, type:

```
./install_packages.sh 
```

It should just work, and you'll be prompted with a bunch of [y]/n questions, for which you should consistently hit "y" then the "enter key" for. If this doesn't work, and you get some kind of an error instead of the prompts, type each of these lines (or copy/paste) one by one, hitting enter after each one:

```
 conda install -c conda-forge kivy
 conda install -c conda-forge biopython=1.76
 conda install -c conda-forge scipy
 conda install pandas
 conda install numpy
 conda install matplotlib
 conda install scikit-learn
```

6) Everything should now be installed, you should now be able to open up the software! Navigate to the app in terminal by typing: 

```
 cd app
```

7) Finally, enter the command

```
 python AIMS.py
```
From there, the GUI should open. A step by step instruction guide for GUI usage will be shortly uploaded.

# Acknowledgements
This initial stable build has only been possible through a highly collaborative effort focused on creating a functional pipeline with the incorporated features necessary for thorough repertoire analysis. The following people have made major contributions to this project:

#### Anna Borowska - UX/UI Designer
As the sole graphic designing consultant on the project, Anna helped ensure the initial application was functional, efficient, and as easy to use as possible. Multiple quality control features within the application were added based on suggestions from Anna, and her comprehensive user testing reports ensured critical functionalities were improved. To contact Anna or view her portfolio, see her webiste.

#### Marta Borowska - Scientific Consultant/Beta Tester
The initial scientific inquiry behind this project was initiated by Marta, looking at the biochemical mechanisms behind antibody polyreactivity. In addition to her scientific insights, Marta was an early tester of the antibody/TCR software module.

#### Caitlin Castro - Scientific Consultant/Beta Tester
Caitlin, the resource for all things immunology, helped with the initial design and testing of the MHC-like software module. Caitlin also helped with interpretation of results and had critical insights for the improvement of the user experience.

#### Erin Adams & Lab - Mentorship/Funding
Erin's expertise and helpful mentoring allowed this project to reach its full potential. Erin was fully supportive of the addition of a new project to the lab, and had great comments for making the analysis pipeline more accessible to experimentalists. Erin's advice has been invaluable throughout. For more information on Erin's lab, see http://ejadamslab.bsd.uchicago.edu/. This project was supported with funding from Erin's National Institutes of Health NIAID grants R01 AI115471 and R01 AI147954.

#### Benoit Roux & Lab - Mentorship/Funding
Benoit provided consistent feedback throughout the project, providing his expertise and insights for all matters biophyscial. Specifically, his insights into the general biophysics of protein-protein interactions helped guide this project. More information on Benoit's lab can be found at http://thallium.bsd.uchicago.edu/RouxLab/. This project was supported with funding from Benoit's National Science Foundation grant MCB-1517221.

#### Martin Flajnik & Lab - MHC Module Collaborators
The MHC module of this project was started due to a collaboration with the Flajnik lab, and sequences used in the mhc_testData directory were isolated by this research group. More information on the Flajnik lab can be found at: https://www.medschool.umaryland.edu/profiles/Flajnik-Martin/ . If utilizing any of the MHC data provided in the mhc_testData directory, please cite:

Almeida T., Esteves P.J., Flajnik M.F., Ohta Y., and Veríssimo A. An Ancient, MHC-Linked, Nonclassical Class I Lineage in Cartilaginous Fish. Journal of Immunology. 2020
