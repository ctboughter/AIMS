AIMS Command Line Interface
=====
After a few years of using AIMS, I have found that there are a *lot* of use cases where I would like to repeat an analysis with a few slight tweaks in either cluster size, alignment scheme, or data subset analyzed. This becomes a bit of a chore in the AIMS Notebooks, because you need to repeatedly find the precise lines in the code that you want to change and repeatedly tweak those, waiting for your analysis to finish. Which is where the command line interface (CLI) comes in: for when you're past the data analysis step and want to move on to the "production" phase and start making publication-quality figures by scanning through some of the metadata space. This is still a work in progress, so if anyone out there has suggestions please drop them on the GitHub or contact me directly.

.. _cliIntro:

An Introduction to the AIMS CLI
------------

As of AIMS v0.8, there *is* a functioning CLI, and some example use cases. AIMSv0.8 had "aims_run.sh" gives a bash script example of some of the options that are available to use the "aims_cli.py" for every type of analysis that AIMS offers. 

With AIMS v0.9 and the pip-based installation, the CLI now runs directly from the terminal rather than needing to call the python script explicitly. The examples that were in the aims_run.sh script can now be found below using the data discussed in the :doc:`Testing` section:

**Antibody Analysis**

.. code-block:: python
    
    aims-cli \
    --datDir test_data/abs \
    --outputDir AIMS_ab \
    --fileNames flu_mono.csv \
    --datNames flu \
    --numLoop 6 \
    --molecule ig > aims_ab.out

**TCR Analysis**

.. code-block:: python
    
    aims-cli \
    --datDir test_data/tcrs \
    --outputDir AIMS_tcr \
    --fileNames siv_tl8.csv siv_cm9.csv \
    --datNames TL8 CM9 \
    --DOstats True \
    --DOboot True \
    --Plotprops True \
    --numLoop 1 \
    --REnorm False \
    --analysisSel metadata \
    --molecule ig > aims_tcr.out

**Peptide Analysis**

.. code-block:: python

    aims-cli \
    --molecule peptide \
    --align bulge \
    --bulgePad 6 \
    --DOstats True \
    --DOboot True \
    --AAorder 'WFMLIVPYHAGSTDECNQRK' \
    --datDir test_data/peptides \
    --outputDir AIMS_pep \
    --analysisSel metadata \
    --fileNames pancreas_hla_atlas.csv kidney_hla_atlas.csv \
    --datNames Pancreas Kidney > aims_pep.out

**MSA Analysis**

.. code-block:: python
    
    aims-cli \
    --datDir test_data/mhcs \
    --outputDir AIMS_mhc \
    --molecule MSA \
    --fileNames cd1.fasta classIa.fasta fish.fasta \
    --datNames CD1 ClassIa Fish \
    --align center \
    --subset True \
    --subStart 164 214 275 327 \
    --subEnd 214 275 327 376 \
    --dropDup True \
    --normProp True \
    --clustData avg \
    --projAlg pca \
    --umapSeed 42 \
    --clustAlg kmean \
    --clustSize 3 \
    --metaForm category \
    --metaName Dset \
    --showProj both \
    --showClust both \
    --normBar False \
    --analysisSel metadata \
    --saveSeqs True \
    --selDat 0 1 2 \
    --prop1 3 \
    --prop2 4 \
    --colors purple blue red \
    --matSize 5 > aims_msa.out

There are a LOT of options and flags to be used with this CLI (33 by my last count), and I don't even think that there are enough programmed in yet. For now, I would advise opening up "aims_cli.py" in your favorite editor to look at all of the possible options, and using the above scripts as templates to make your own scripts. A little bash programming would go a long way (on the user end) to help run through a bunch of possible CLI options.

And, of course, you can see all of the options that are available in the AIMS CLI using:

.. code-block:: python
    
    aims-cli --help

.. _cliOptions:

CLI Options
------------

Further info on the AIMS command line interface coming soon!

**Automating analysis**

**Generating reproducible analysis**

**Restarting from Saved Analysis**

**Parallelization**