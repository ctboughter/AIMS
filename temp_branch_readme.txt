Creating a special README just for this branch to keep track of what I'm working on:

NOTE:
For now, we need to keep the "dev" files on the Github, because the non-dev versions rely on the pip install.
Obviously we don't want to be pushing changes to the pip install from this branch.
When we're done developing, we can delete the "dev" files after merging with the originals.

COMPLETED:
Fix the mat_size bug (get_sequence_dimension command)

TODO:
Make the code cleaner: Shouldn't need to deal with these "dev" versions of *every* python file. Should be independent of each other.
Get numba parallelization to work
Continuous AIMS embedding (i.e. "stretch out" amino acids to make individual sequences "blurrier")
Get on to conda
Create a bash script for automating sequence clustering