AIMS Jupyter Notebooks
=====
While the GUI is useful for those with little to no experience using python, unfortunately given the time commitment required to create a functional GUI the available features incorporated can lag behind significantly. For the most up-to-date analysis pipelines, fastest code, and maximal flexibility, it is recommended that the Jupyter Notebooks are used. I have tried to explain step by step within the notebook how to run it. It may take a while, but reading through the markdowns and comments should give the user a good feel for what each code block is doing.

It should be noted that previous versions of AIMS had a separate notebook for every possible application. To reduce clutter of the repository, these separate notebooks have been combined into a single notebook (as of AIMS v0.7.5), which is in turn front loaded with instructions for each specific application of AIMS. There are example files included for every possible analysis mode of AIMS. Here we will include some general instructions for using notebooks, and provide some details for the "best practices" of using the AIMS Notebooks.

.. _notes:

Notes on Notebooks
------------
If using Anaconda-Navigator, you should just be able to launch Jupyter Lab from the main application. If using terminal, you can install Jupyter lab (and associated packages) within your AIMS environment (discussed in :doc:`Install`) using:

.. code-block:: python

    conda install -c conda-forge jupyterlab=4.0.4
    conda install -c anaconda ipython=8.14.0
    conda install -c anaconda ipykernel=6.25.1

Copying each of these commands in line-by-line. There is an added "trick" which you need in order to access your AIMS environment in Jupyter lab. Again with your AIMS environment activated, and assuming you named your environment "aims-env" as is suggested in the install instructions, enter the following command:

.. code-block:: python
    
    python -m ipykernel install --user --name aims-env --display-name "aims-env"

You should now be ready to use the notebook! Start using the notebook by typing in:

.. code-block:: python
    
    jupyter lab


A new window should open automatically, and you should see the "AIMS_notebook" file on the left-hand side. Double click to open and start running the notebook! You may also need to change the kernel, which can be found near the upper-right corner of the window. The default is likely "Python 3" which you can click to open up a dropdown window where you will (hopefully) find your aims-env kernel listed. Select this and you should be fully ready. Remember that "ctrl+Enter" is how you can run each individual cell, or if this is your first time running the notebook you could also select the option to run all the code at once so you can test if anything in your particular build is broken. For more instructions using Jupyter lab and notebooks, see some random help pages such as (https://www.datacamp.com/tutorial/installing-jupyter-notebook). I unfortunately did not find the Jupyter documentation to be terribly helpful, so it seems like third-party instruction may be better here.

.. _bookOptions:

Notebook Usage
------------

The notebooks are the most frequently updated parts of this software because they are probably the most useful form of the entire analysis package. AIMS is largely meant to be an exploratory tool used for interrogating large repertoire datasets. The flexibility and ability to dig deeper into certain aspects of the data is unique to the Notebooks, and is something that users should take advantage of. Once users are comfortable with their given dataset, the hope is that they can then automate much of the analysis using the CLI (see :doc:`AIMS_CLI`), which is currently under active development.