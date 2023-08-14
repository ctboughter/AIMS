Installation Instructions
=====

.. _beg.install:

Initial Installation for Beginners
------------

These instructions will help novice programmers install all necessary packages and programs needed to run a python program. Mac/Linux OS preferred, these steps outline installation via Anaconda. Other installations should be supported but have had limited testing.

.. note::
    If you are a more advanced programmer, skip to the :ref:`aims.install` section

1. Install Anaconda (https://www.anaconda.com/products/individual) to manage the python packages we're going to be using. This can be a fairly large package, so if space is at a premium for your computer, you can instead install miniconda (https://docs.conda.io/en/latest/miniconda.html). NOTE: Windows users should likely install the full Anaconda package, for a contained environment to run python programs from.

2. Test that your conda install is working properly by creating a conda environment. Windows OS users, you will likely do this within the Anaconda application. Mac/Linux users, open the terminal application. Once terminal is open, type:

.. code-block:: python

    conda create -n aims-env python=3.9

If anaconda/miniconda is installed properly, a Y/N prompt should appear. Type "y" then hit the "enter key" and you will create a conda environment.

3. If you haven't already, you should download the code from https://github.com/ctboughter/AIMS. To download from GitHub, click the green "Code" button, and then "Download Zip". Then, unzip the folder and move the AIMS directory to whichever location you would like to run the analysis from. Alternatively you can also download via terminal using this command:

.. code-block:: python

    git clone https://github.com/ctboughter/AIMS.git


4. Next, navigate to the new folder created from this repository. First, open up terminal and enter the environment created earlier by typing:

.. code-block:: python

    conda activate aims-env

.. note::
    Windows users will have to activate the environment from within Anaconda

You should now see a little extra bit of text on your terminal command line that looks something like "(aims-env)". If this didn't work for some reason, an error message should pop up, otherwise assume you're fine.

Use terminal to navigate into the newly downloaded folder. If you've never used terminal before, you can type in "cd" and then drag and drop the folder into the terminal. Doing so should automatically populate the "path" to the folder. Then hit enter.

When I do this, my terminal line reads: 

.. code-block:: python
    cd /Users/boughter/Desktop/AIMS-master

Hopefully you see something similar (replacing my user name with your own, and noting that "AIMS-master" may simply be "AIMS" if you cloned the repository rather than downloading the zip), and when you hit enter you can move to the next step.

.. _aims.install:

Installation of AIMS & Required Packages
------------

The above section focused on the installation of the necessary software to run python programming, as well as basic navigational instructions. In this section, instructions for installing requisite packages and running the GUI Are discussed.

1. From within the AIMS master directory, type in the terminal:

.. code-block:: python

    ./app/install_packages.sh

This bash script should run after typing in this command, and you'll be prompted with a bunch of [y]/n prompts, for which you should consistently enter "y" then the "enter key". 

.. note::
    For more advanced users, you shouldn't need to use these EXACT package versions, save for Biophython and SciKit-Learn. However, using these versions should guarantee proper functionality of AIMS

If the install_packages.sh script doesn't work, and you get some kind of an error instead of the prompts, type each of these lines (or copy/paste) one by one, hitting enter after each one:

.. code-block:: python

    conda install -c conda-forge umap-learn=0.5.3
    conda install -c conda-forge biopython=1.79
    conda install -c conda-forge scipy=1.4.1
    conda install pandas=1.5.3
    conda install numpy=1.24.1
    conda install matplotlib=3.7.1
    conda install scikit-learn=1.3.0
    conda install seaborn=0.12.2
    conda install -c conda-forge kivy=2.1.0

.. warning::
    The versions of these apps have been updated as of AIMS v0.8 to ensure AIMS is not using outdated packages. However, not every function has been tested, so please do not hesitate to raise issues on GitHub if something is non-functional with these new packages.

    Further, if you do not plan on using the GUI, do not install Kivy. It seems to be the source of trouble for most installs, and is only used to run the GUI.

2. Everything should now be installed, you should now be able to open up the software! If you'd like to run the notebook or the command-line interface, please see the :doc:`AIMS_notebooks` or the :doc:`AIMS_CLI` sections. If you would like to run the GUI, navigate to the app in terminal by typing:

.. code-block:: python

    cd app

3. Launch the GUI with:

.. code-block:: python 

    python aims.py

From there, the GUI should open. A step by step instruction guide for GUI usage can be found in the :doc:`AIMS_GUI` section. If you don't want to be bothered reading instructions, the app should prevent most major errors. If a "next" button is grayed out, make sure you've pressed all of the analysis buttons on the bottom of the current AIMS app screen.

If you're a more advanced user and would prefer a more customizable experience, check out the :doc:`AIMS_notebooks` section.

Lastly, if you're generally interested in an overview of what AIMS does and how it works, refer to the :doc:`AIMS_basics`.