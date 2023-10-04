AIMS Jupyter Notebooks
=====
While the GUI is useful for those with little to no experience using python, unfortunately given the time commitment required to create a functional GUI the available features incorporated can lag behind significantly. For the most up-to-date analysis pipelines, fastest code, and maximal flexibility, it is recommended that the Jupyter Notebooks are used. I have tried to explain step by step within the notebook how to run it. It may take a while, but reading through the markdowns and comments should give the user a good feel for what each code block is doing.

It should be noted that previous versions of AIMS had a separate notebook for every possible application. To reduce clutter of the repository, these separate notebooks have been combined into a single notebook (as of AIMS v0.7.5), which is in turn front loaded with instructions for each specific application of AIMS. There are example files included for every possible analysis mode of AIMS. Here we will include some general instructions for using notebooks, and provide some details for the "best practices" of using the AIMS Notebooks.

The AIMS notebook comes pre-loaded with scripts to test AIMS functionalities with example data, which can be found by following the instructions available in the :doc:`Testing` section.

.. _notes:

Notes on Notebooks
------------
Now that AIMS can be installed using pip, there is no longer a need to separately install jupyterlab and ipython. If you'd still like to see the precise package versions we use, go to the bottom of this page. You can now directly launch the jupyter notebook simply entering into the terminal:

.. code-block:: python

    aims-notebook

Since notebooks are a little tricky to launch from a package distributed via pip, AIMS uses something of a cheat by simply copying the notebook from the AIMS package directory to your current directory. So, after you use the above command once, it is recommended that you launch the notebook using:

.. code-block:: python

    jupyter lab

This will open up whatever directory you are currently in, where you can now hopefully open the file AIMS_notebook.ipynb, which is in your current directory because you ran the "aims-notebook" command at some point in the past.

.. warning::
    If you make any custom changes to your copied AIMS_notebook.ipynb file, it is *very important* that you either change the name of the notebook or *not* run the aims-notebook command. There is a chance this will overwrite your changes.

There may be an added "trick" which you need in order to access your AIMS environment in Jupyter lab. Again with your AIMS environment activated (see :doc:`Install` for a reminder on how to do this), and assuming you named your environment "aims-env" as is suggested in the install instructions, enter the following command:

.. code-block:: python
    
    python -m ipykernel install --user --name aims-env --display-name "aims-env"

When you launch the notebook (either using aims-notebook or jupyter lab) a new window should open automatically, and you should see the "AIMS_notebook" file on the left-hand side. Double click to open and start running the notebook! You may also need to change the kernel, which can be found near the upper-right corner of the window. The default is likely "Python 3" which you can click to open up a dropdown window where you will (hopefully) find your aims-env kernel listed (assuming you did the above "trick". Select this and you should be fully ready. Remember that "ctrl+Enter" is how you can run each individual cell, or if this is your first time running the notebook you could also select the option to run all the code at once so you can test if anything in your particular build is broken. For more instructions using Jupyter lab and notebooks, see some random help pages such as (https://www.datacamp.com/tutorial/installing-jupyter-notebook). I unfortunately did not find the Jupyter documentation to be terribly helpful, so it seems like third-party instruction may be better here.

.. _bookOptions:

Notebook Usage
------------

The notebooks are the most frequently updated parts of this software because they are probably the most useful form of the entire analysis package. AIMS is largely meant to be an exploratory tool used for interrogating large repertoire datasets. The flexibility and ability to dig deeper into certain aspects of the data is unique to the Notebooks, and is something that users should take advantage of. Once users are comfortable with their given dataset, the hope is that they can then automate much of the analysis using the CLI (see :doc:`AIMS_CLI`), which is currently under active development.

**Precise Package Versions**
Not exactly where else to put this, so we'll leave it here for now (sorry for the brutal honesty):

.. code-block:: python

    conda install -c conda-forge jupyterlab=4.0.4
    conda install -c anaconda ipython=8.14.0
    conda install -c anaconda ipykernel=6.25.1