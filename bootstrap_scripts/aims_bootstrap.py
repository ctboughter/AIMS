#!/bin/python
import numpy as np
import pandas
import math
from sklearn.utils import resample
import time
import sys
import aims_analysis as aims
import random
import aims_loader as aimsLoad
import multiprocessing as mp

# I might want to get rid of this stuff to make it a more tunable script...
import seq_loader

# references to ntasks and "which_num" are for doing this over multiple
# CPUs, which we don't really need to worry about on local machines
#ntasks = int(sys.argv[1])
#which_num = int(sys.argv[2])

num_RE = 10

# How many times will this replicant run
#num_iter = int(np.ceil(num_RE/ntasks))
#print('We will run '+str(num_iter)+' iterations in process '+str(which_num))

ALL_mono = aimsLoad.Ig_loader('/Users/boughterct/Desktop/AIMS-immunopep/app/ab_testData/flu_mono.csv',
label='treg', loops = 6, drop_degens = True)
ALL_poly = aimsLoad.Ig_loader('/Users/boughterct/Desktop/AIMS-immunopep/app/ab_testData/flu_poly.csv',
                              label='conv',loops = 6, drop_degens = True)

newnew=pandas.read_csv('app_data/new_props')
oldold=pandas.read_csv('app_data/old_props')

# Again, ugly to hard code in the number of properties (62) but 
# For now no harm no foul
AA_key=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
properties=np.zeros((62,20))
for i in np.arange(len(AA_key)):
    properties[0:16,i]=oldold[AA_key[i]]
    properties[16:,i]=newnew[AA_key[i]]

AA_num_key_new=properties[1]
AA_num_key=np.arange(20)+1

# Might want to change this to like half of the min,
# in order to actually subsample the data.
#NUMsamples = min(len(ALL_mono[0]),len(ALL_poly[0]))

mono_size = len(np.transpose(ALL_mono))
poly_size = len(np.transpose(ALL_poly))

test_MATRIX1 = ALL_mono
test_MATRIX2 = ALL_poly

# Really need to pre-define the matrix size based on the largest matrix                                         
poly_matrix_MI,mono_matrix_MI,matSize=aims.gen_tcr_matrix(np.array(ALL_poly),pre_mono=np.array(ALL_mono),key = AA_num_key,binary=True,return_Size=True)

def boot_it(RE):
    # REPLACE NUMsamples with the actual size of each dataset... (They're close anyway)
    remono=np.array(np.transpose(resample(np.transpose(test_MATRIX1),random_state=int(1000*random.random()),n_samples = mono_size)))
    repoly=np.array(np.transpose(resample(np.transpose(test_MATRIX2),random_state=int(1000*random.random()),n_samples = poly_size)))
    mono_matrix_MI,poly_matrix_MI=aims.gen_tcr_matrix(np.array(remono),pre_mono=np.array(repoly),key = AA_num_key,binary = True,giveSize = matSize)
    len_distM=np.zeros((6,len(remono[0])))
    len_distP=np.zeros((6,len(repoly[0])))
    for j in [0,1,2,3,4,5]: # This one for light
        for k in np.arange(len(remono[0])):
            len_distM[j,k]=len(remono[j][k])
        for k in np.arange(len(repoly[0])):
            len_distP[j,k]=len(repoly[j][k])
    len_meanM=np.average(len_distM,axis=1)
    len_meanP=np.average(len_distP,axis=1)
    len_stdM=np.std(len_distM,axis=1)
    len_stdP=np.std(len_distP,axis=1)

    clone_meansM=aims.gen_clone_props(mono_matrix_MI)
    clone_meansP=aims.gen_clone_props(poly_matrix_MI)
    mean_varsM=np.average(clone_meansM,axis=1)
    mean_varsP=np.average(clone_meansP,axis=1)
    std_varsM=np.std(clone_meansM,axis=1)
    std_varsP=np.std(clone_meansP,axis=1)

    mono_dset = aims.gen_dset_props(mono_matrix_MI,stdev=True); poly_dset = aims.gen_dset_props(poly_matrix_MI,stdev=True)
    mono_MI,mono_count_cond,counted_mono = aims.calculate_MI(mono_matrix_MI)
    poly_MI,poly_count_cond,counted_poly = aims.calculate_MI(poly_matrix_MI)
    # No need to do these anymore, but leave it as is in the hopes that
    # this will keep it consistent with supercomputer bootstraps.
    dsetM = [mono_dset]; dsetP = [poly_dset]
    MI_m = [mono_MI]; MI_p = [poly_MI]
    count_m = [counted_mono]; count_p = [counted_poly]
    sin_propM = [mean_varsM]; sin_propP = [mean_varsP]
    sin_propMstd = [std_varsM]; sin_propPstd = [std_varsP]
    lenM = [len_meanM]; lenP = [len_meanP]
    lenMstd = [len_stdM]; lenPstd = [len_stdP]
    
    return(dsetM, MI_m, count_m, sin_propM, sin_propMstd, lenM, lenMstd,
        dsetP, MI_p, count_p, sin_propP, sin_propPstd, lenP, lenPstd)

def do_boot(numbers):
    with mp.Pool() as pool:
        results=pool.map(boot_it, numbers)
        # results is broken down as results[iteration][calculated property]
        return(results)

if __name__ == "__main__":
    start_time = time.time()
    x=np.arange(0,num_RE)
    aa = do_boot(x)
    #    mono_props,mono_MI,poly_props,poly_MI,len_mono,len_poly,mean_mono,mean_poly = do_boot(x)
    duration = time.time() - start_time
    print(duration)

np.save('bootstrap_ALL',aa)
np.save('bootstrap_size_ALL',matSize)