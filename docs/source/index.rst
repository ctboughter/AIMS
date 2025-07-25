Using AIMS to Analyze Immune Repertoires
===================================

.. figure:: images/ania1.png
   :alt: AIMS Graphic By Anna Borowska

   AIMS Graphic created by Anna Borowska (see https://annazofiaborowska.com/)

.. note::
   This readthedocs site is still under development. Hopefully what is included thus far helps users, but more information is coming soon!

AIMS, an Automated Immune Molecule Separator, was originally developed to identify differences between two distinct antibody repertoires, but has since expanded to become a multi-purpose repertoire analysis tool. Currently the AIMS analytical tools can be applied to immunoglobulin (Ig) molecules such as T cell receptors and antibodies, major histocompatibility complex (MHC) and MHC-like molecules, immunopeptidomic data, and broadly to any data presented in a multi-sequence alignment format. This documentation will teach users how to use AIMS to:

- **Get started** following the :doc:`Install`.


- **Analyze repertoires with no programming experience required** in a user-friendly format with the :doc:`AIMS_GUI`.


- **Characterize the key features of immune repertoires** using the AIMS biophysical characterization tools through the :doc:`AIMS_basics` section. This biophysical characterization is the central feature of AIMS, taking properly formatted files (see :ref:`formatting`) and applying analytical tools from the :ref:`core` of AIMS for downstream manipulation.


- **Identify biophysically distinct repertoire clusters** within single repertoires or comparing across multiple repertoires using the :doc:`AIMS_cluster` features.


- **Build off the AIMS analysis with your own custom features** by taking advantage of the :doc:`AIMS_notebooks`.

AIMS is a python package distributed in a notebook, CLI, and GUI format. Those wishing to use the GUI, particularly those relatively new to programming, can follow the installation instructions. Example data is provided in the test_data directories, and an example of an application of AIMS can be seen in this peer-reviewed article: https://elifesciences.org/articles/61393
   
.. note::
   When publishing analysis from this software, please cite:

   Boughter CT, Borowska MT, Guthmiller JJ, Bendelac A, Wilson PC, Roux B, Adams EJ. Biochemical Patterns of Antibody Polyreactivity Revealed Through a Bioinformatics-Based Analysis of CDR Loops. eLife. 2020. DOI: 10.7554/eLife.61393

   &

   Boughter CT, Meier-Schellersheim M. An Integrated Approach to the Characterization of Immune Repertoires Using AIMS\: An Automated Immune Molecule Separator. PLoS Computational Biology. 2023. DOI: 10.1371/journal.pcbi.1011577

Contents
--------

.. toctree::
   :maxdepth: 2

   Install
   AIMS_GUI
   AIMS_basics
   AIMS_cluster
   AIMS_notebooks
   AIMS_CLI
   Testing
   Acknowledgements
   
A Note on Exploration with AIMS
--------

The AIMS software should be considered as both a means for exploratory searches through data to generate hypotheses and as a tool for rigorous quantification of differences between molecular subsets. In the former application, users can freely explore their data, tuning different AIMS parameters and seeing how these changes alter identified clusters or comparison groups. However, in the latter application, users should carefully record the setting of each tuned parameter. Analysis using AIMS should be considered akin to modern RNAseq analysis, where the rigor of a given analytical tool depends on proper implementation by the user. Reproducibility is key!

Further Reading
--------

Now that AIMS has been out and in the wild for around five years, there have been additional published peer-reviewed manuscripts or posted preprints that highlight the capabilities of AIMS! I'll try to keep this list relatively up to date. Manuscripts thus far include:

- The manuscript that started it all, using AIMS to assess antibody polyreactivity: https://elifesciences.org/articles/61393 
- A targeted application of AIMS to identify a clonally expanded, phenotypically homogenous T cell population: https://www.nature.com/articles/s41590-025-02198-4
- A [preprint] textbook chapter giving step-by-step instructions for running the analysis on peptides: https://www.biorxiv.org/content/10.1101/2024.09.05.611486v1 
- An application of AIMS to non-immune molecules using multi-sequence alignment (MSA) encoding: https://pubs.acs.org/doi/abs/10.1021/acs.jpcb.2c02173
- A unique paper highlighting the flexiblity of AIMS, with an application to an analysis of SARS-CoV-2 binding epitopes: https://www.nature.com/articles/s42003-023-05332-w 
- Technically a paper that doesn't use AIMS, but does provide experimental validation of AIMS-based predictions: https://www.cell.com/cell-reports/fulltext/S2211-1247(23)01202-0
- The AIMS bible, with a thorough explanation of the rationale behind the AIMS analysis: https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1011577
- An investigation of the nature of the germline interactions between TCR CDR loops and MHC: https://elifesciences.org/articles/90681