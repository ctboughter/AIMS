# Yet another wild way to do this, but way easier than
# trying to configure things within each app/cli to find
# the test data... Plus it lets users play with it a bit more
import aims_immune
import os
import shutil

def getit():
    aims_dir = aims_immune.__file__[:-11]
    startDir = os.getcwd()
    shutil.copytree(aims_dir+'app_data/test_data',startDir+'/test_data')