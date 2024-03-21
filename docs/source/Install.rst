Installation & Startup Instructions
=====

.. _new.install:

Installing AIMS via PyPI
------------

.. note::
    If you are new to programming and/or python, it might be worth reading :ref:`beg.install` section first.

As of AIMSv0.9, users can now simply install AIMS using "pip" to download the software from PyPI. This will install all the necessary packages and the software itself, and allow users to launch the AIMS GUI, CLI, or notebook from the command line without additional complicated steps. The quick install step is simply entering:

.. code-block:: python

    pip install aims-immune

into your terminal. If you would prefer a specific version of AIMS, say this first version available on PyPI (v0.9), you can instead specify the version using:

.. code-block:: python

    pip install aims-immune==0.9

.. warning::
    Currently, AIMS requires python v3.9. Future AIMS versions will hopefully allow for more flexible installs.

Note that there are versions available on PyPI before v0.9, but none of these contain the completed package. In other words, do not use pip to install versions of AIMS older than v0.9. Instead, go to the AIMS github (https://github.com/ctboughter/AIMS) and use the "tags" feature to select older versions.

The different AIMS wrappers can conveniently be launched directly from the command line using one of the below commands:

.. code-block:: python

    aims-gui

.. code-block:: python

    aims-notebook

.. code-block:: python

    aims-cli

For more details on how to use each of these wrappers, see :doc:`AIMS_GUI`, :doc:`AIMS_notebooks`, or :doc:`AIMS_CLI`, for each of these commands, respectively.

.. note::
    Unfortunately, Windows users do not have a nice Linux-based terminal to enter these commands into. Instead, read below to learn how to install Anaconda, where you can likely use the Qt Console to effectively emulate the terminal. Sadly I don't have a windows machine, so can't test this.

.. _beg.install:

Installation Notes for Beginners
------------

While users can simply use pip to install AIMS directly, it is best practice to install AIMS in a self-contained environment. The python package Kivy, which is used to run the GUI, tends to cause issues when installing other python packages. The self-contained AIMS environment will help alleviate this issue. Read this section *before* installing AIMS using pip. These steps will show you how to do this using Anaconda. Mac/Linux OS preferred. Other installations should be supported but have had limited testing.

1. Install Anaconda (https://www.anaconda.com/products/individual) to manage the python packages we're going to be using. This can be a fairly large package, so if space is at a premium for your computer, you can instead install miniconda (https://docs.conda.io/en/latest/miniconda.html). Windows users should likely install the full Anaconda package, for a contained environment to run python programs from.

2. Test that your conda install is working properly by creating a conda environment. Windows OS users, you will likely do this within the Anaconda application (probably using Qt Console). Mac/Linux users, open the terminal application. Once terminal is open, type:

.. code-block:: python

    conda create -n aims-env python=3.7

If anaconda/miniconda is installed properly, a Y/N prompt should appear. Type "y" then hit the "enter key" and you will create a conda environment.

3. Next, "enter" the environment you just created by typing in the terminal:

.. code-block:: python

    conda activate aims-env

You should now see a little extra bit of text on your terminal command line that looks something like "(aims-env)". If this didn't work for some reason, an error message should pop up, otherwise assume you're fine.

4. Use terminal to navigate into the directory with the data you'd like to analyze. If you've never used terminal before, you can type in "cd" and then drag and drop the folder into the terminal. Doing so should automatically populate the "path" to the folder. Then hit enter.

When I do this, my terminal line reads: 

.. code-block:: python
    
    cd /Users/boughter/Desktop/myData

Hopefully you see something similar (replacing my user name with your own, and noting that "myData" is of course replaced with your data folder name).

5. Install AIMS and run the analysis! As highlighted in the above :ref:`new.install` section.

Best of luck with your programming journey! Hope this was a useful introduction to using Anaconda to create environments.

.. _aims.install:

Installing AIMS the Old Fashioned Way
------------

While the pip install is very useful and convenient, some users may want more control over their installation or would prefer to install a version of AIMS that predates v0.9. The "old" steps for installing via GitHub are included here.

1. Start by downloading the code from https://github.com/ctboughter/AIMS. To download from GitHub, click the green "Code" button, and then "Download Zip". Then, unzip the folder and move the AIMS directory to whichever location you would like to run the analysis from. Alternatively you can also download via terminal using this command:

.. code-block:: python

    git clone https://github.com/ctboughter/AIMS.git


Accessing previous versions of AIMS using the "git clone" option is a little tricky, so it is recommended you download the zip from the website by navigating to your version of interest using the "tags".

The dependencies are as follows (for python3.9). See previous versions of this ReadTheDocs page for the versions that are compatible with python3.7.

.. warning::
    The versions of these apps have been updated as of AIMS v0.8 to ensure AIMS is not using outdated packages. However, not every function has been tested, so please do not hesitate to raise issues on GitHub if something is non-functional with these new packages.

    Further, if you do not plan on using the GUI, do not install Kivy. It seems to be the source of trouble for most installs, and is only used to run the GUI.

.. code-block:: python

    conda install -c conda-forge kivy=2.1.0
    conda install -c conda-forge umap-learn=0.5.3
    conda install -c conda-forge biopython=1.79
    conda install -c conda-forge scipy=1.4.1
    conda install pandas=1.5.3
    conda install numpy=1.24.1
    conda install matplotlib=3.7.1
    conda install scikit-learn=1.3.0
    conda install seaborn=0.12.2

2. If you are installing via GitHub, then the functions for calling the CLI, GUI, or notebook directly from the terminal will not work. In AIMS v0.9 and above, you can launch these wrappers in the following ways (assuming the downloaded GitHub directory is called "AIMS", and you have navigated into the directory which holds AIMS):

.. code-block:: python

    python AIMS/aims_immune/aims_cli.py

.. code-block:: python

    python AIMS/aims_immune/aims.py

.. code-block:: python

    jupyter lab AIMS/aims_immune/AIMS_notebook.ipynb

3. If you are instead downloading AIMS v0.8, look at the individual versions for the directory structure, as this has changed a bit over the different versions.

A step by step instruction guide for GUI usage can be found in the :doc:`AIMS_GUI` section. If you don't want to be bothered reading instructions, the app should prevent most major errors. If a "next" button is grayed out, make sure you've pressed all of the analysis buttons on the bottom of the current AIMS app screen.

If you're a more advanced user and would prefer a more customizable experience, check out the :doc:`AIMS_notebooks` section.

If you're really comfortable with AIMS, check out the :doc:`AIMS_CLI`.

Lastly, if you're generally interested in an overview of what AIMS does and how it works, refer to the :doc:`AIMS_basics`.