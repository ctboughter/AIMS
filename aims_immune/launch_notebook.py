import jupyterlab
# Kind of a wild thing but this will be a nice way 
# to be able to launch the notebook automatically without
# our users having to search for it... maybe.
from jupyterlab import labapp
import aims_immune
import os
import shutil

def launchit():
    aims_dir = aims_immune.__file__[:-11]
    startDir = os.getcwd()
    shutil.copyfile(aims_dir+'AIMS_notebook.ipynb',startDir+'/AIMS_notebook.ipynb')
    labapp.main()