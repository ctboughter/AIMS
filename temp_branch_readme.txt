Creating a special README just for this branch to keep track of what I'm working on:

NEW NOTE: I don't need to use this "dev" version of the files
And, remember, before that it is "python -m build" to create the pip package to then install from.
Instead, I can use "python -m pip install -e ." within my AIMS directory.
For doing this, I've created an AIMS_dev environment

Extra notes from a pip explainer (see here https://setuptools.pypa.io/en/latest/userguide/development_mode.html):
"When you install a package in editable mode, you’re creating a link in the site-packages to the local project path:"

As a key note though, when we are dealing with these editable installs, they are not guaranteed to be the same as a normal install
So, you should still test in a fresh test_env without the editable flag.

COMPLETED:
Fix the mat_size bug (get_sequence_dimension command)
Added numba to "get bigass matrix" but there are still multiple other places to add it 

TODO:
Get numba parallelization to work
Continuous AIMS embedding (i.e. "stretch out" amino acids to make individual sequences "blurrier")
Get on to conda
Create a bash script for automating sequence clustering