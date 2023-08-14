AIMS Command Line Interface
=====
After a few years of using AIMS, I have found that there are a *lot* of use cases where I would like to repeat an analysis with a few slight tweaks in either cluster size, alignment scheme, or data subset analyzed. This becomes a bit of a chore in the AIMS Notebooks, because you need to repeatedly find the precise lines in the code that you want to change and repeatedly tweak those, waiting for your analysis to finish. Which is where the command line interface (CLI) comes in, for when you're past the data analysis step and want to move on to the "production" phase and start making publication-quality figures by scanning through some of the metadata space. This is still a work in progress, so if anyone out there has suggestions please drop them on the GitHub or contact me directly.

.. _cliIntro:

An Introduction to the AIMS CLI
------------

As of AIMS v0.8, there *is* a functioning CLI, and some example use cases. "aims_run.sh" gives a bash script example of some of the options that are available to use the "aims_cli.py" for every type of analysis that AIMS offers. The AIMS CLI is closer to up-to-date with the AIMS Notebooks compared to the GUI, but as of this update is currently missing statistical analysis (Number 1 priority on the to-do list right now). 

There are a LOT of options and flags to be used with this CLI (33 by my last count), and I don't even think that there are enough programmed in yet. For now, I would advise opening up "aims_cli.py" in your favorite editor to look at all of the possible options, and using "aims_run.sh" as a template to make your own scripts. A little bash programming would go a long way (on the user end) to help run through a bunch of possible CLI options.

.. _cliOptions:

CLI Options
------------

Further info on the AIMS command line interface coming soon!

**Automating analysis**

**Generating reproducible analysis**

**Restarting from Saved Analysis**

**Parallelization**