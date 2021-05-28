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

# I might want to get rid of this shit to make it a more tunable script...
import seq_loader

ntasks = int(sys.argv[1])
which_num = int(sys.argv[2])

num_RE = 1000

# How many times will this replicant run
num_iter = int(np.ceil(num_RE/ntasks))
print('We will run '+str(num_iter)+' iterations in process '+str(which_num))

#[mono_all,poly_all,mono,poly]=seq_loader.getBunker()
#[mono_jennaAll,poly_jennaAll,mono_jenna,poly_jenna]=seq_loader.getJenna()
#[mono_hugoAll,poly_hugoAll,mono_Hugo,poly_Hugo]=seq_loader.getHugo()
#[mono_HugNat_all,poly_HugNat_all,mono_HugNat,poly_HugNat]=seq_loader.getHugo_Nature()
#[mono_PLOS,poly_PLOS]=seq_loader.getHugo_PLOS()
#[mono_HugNat_allCNT,poly_HugNat_allCNT,mono_HugNatCNT,poly_HugNatCNT] = seq_loader.getHugo_NatCNTRL()

# flu
#ALL_mono = mono_jenna
#ALL_poly = poly_jenna

#mouse
#ALL_mono = mono
#ALL_poly = poly

# HIV ALL
#ALL_mono=np.hstack((mono_hugoAll,mono_HugNat_all,mono_PLOS,mono_HugNat_allCNT))
#ALL_poly=np.hstack((poly_hugoAll,poly_HugNat_all,poly_PLOS,poly_HugNat_allCNT))

# HIV
#ALL_mono=np.hstack((mono_Hugo,mono_HugNat,mono_HugNatCNT,mono_PLOS)) 
#ALL_poly=np.hstack((poly_Hugo,poly_HugNat,poly_HugNatCNT,poly_PLOS))

# All Data 
#ALL_mono=np.hstack((mono_all,mono_jennaAll,mono_hugoAll,mono_HugNat_all,mono_HugNat_allCNT,mono_PLOS))
#ALL_poly=np.hstack((poly_all,poly_jennaAll,poly_hugoAll,poly_HugNat_all,poly_HugNat_allCNT,poly_PLOS))

# Human Data                                                                
#ALL_mono=np.hstack((mono_jennaAll,mono_hugoAll,mono_HugNat_all,mono_HugNat_allCNT,mono_PLOS))
#ALL_poly=np.hstack((poly_jennaAll,poly_hugoAll,poly_HugNat_all,poly_HugNat_allCNT,poly_PLOS))

# Parsed Data
#ALL_mono=np.hstack((mono,mono_jenna,mono_Hugo,mono_HugNat,mono_HugNatCNT,mono_PLOS))
#ALL_poly=np.hstack((poly,poly_jenna,poly_Hugo,poly_HugNat,poly_HugNatCNT,poly_PLOS))

ALL_mono = aimsLoad.Ig_loader('/project2/ejadams/boughter/cdrds/WT_Treg_AIMS_allCDRs.csv',
label='treg', loops = 6, drop_degens = False)
ALL_poly = aimsLoad.Ig_loader('/project2/ejadams/boughter/cdrds/WT_Tconv_AIMS_allCDRs.csv',
                              label='conv',loops = 6, drop_degens = False)

newnew=pandas.read_csv('new_props')
oldold=pandas.read_csv('old_props')

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

tot_size = int(poly_size + mono_size)

stack = np.transpose(np.hstack((ALL_mono,ALL_poly)))
np.random.shuffle(stack)

test_MATRIX1 = stack
# Really need to pre-define the matrix size based on the largest matrix                                         
poly_matrix_MI,mono_matrix_MI,matSize=aims.gen_tcr_matrix(np.array(ALL_poly),pre_mono=np.array(ALL_mono),key = AA_num_key,binary=True,return_Size=True)

for i in np.arange(num_iter):
    # REPLACE NUMsamples with the actual size of each dataset... (They're close anyway)
    remono=np.array(np.transpose(resample(np.transpose(stack),random_state=int(1000*random.random()),n_samples = mono_size)))

    shuff_mono = remono[:,0:mono_size]
    shuff_poly = remono[:,mono_size:poly_size+mono_size]

    mono_matrix_MI,poly_matrix_MI=aims.gen_tcr_matrix(np.array(shuff_mono),pre_mono=np.array(shuff_poly),key = AA_num_key,binary = True,giveSize = matSize)

    len_distM=np.zeros((6,len(shuff_mono[0])))
    len_distP=np.zeros((6,len(shuff_poly[0])))
    for j in [0,1,2,3,4,5]: # This one for light
        for k in np.arange(len(shuff_mono[0])):
            len_distM[j,k]=len(shuff_mono[j][k])
        for k in np.arange(len(shuff_poly[0])):
            len_distP[j,k]=len(shuff_poly[j][k])
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
    # THIS IS NOW ALL CHANGED TO INCLUDE STDEV
    print(i)
    mono_dset = aims.gen_dset_props(mono_matrix_MI,stdev=True); poly_dset = aims.gen_dset_props(poly_matrix_MI,stdev=True)
    mono_MI,mono_count_cond,counted_mono = aims.calculate_MI(mono_matrix_MI)
    poly_MI,poly_count_cond,counted_poly = aims.calculate_MI(poly_matrix_MI)
    if i == 0:
        dsetM = [mono_dset]; dsetP = [poly_dset]
        MI_m = [mono_MI]; MI_p = [poly_MI]
        count_m = [counted_mono]; count_p = [counted_poly]
        sin_propM = [mean_varsM]; sin_propP = [mean_varsP]
        sin_propMstd = [std_varsM]; sin_propPstd = [std_varsP]
        lenM = [len_meanM]; lenP = [len_meanP]
        lenMstd = [len_stdM]; lenPstd = [len_stdP]
    else:
        dsetM.append(mono_dset); dsetP.append(poly_dset)
        MI_m.append(mono_MI); MI_p.append(poly_MI)
        count_m.append(counted_mono); count_p.append(counted_poly)
        sin_propM.append(mean_varsM); sin_propP.append(mean_varsP)
        sin_propMstd.append(std_varsM); sin_propPstd.append(std_varsP)
        lenM.append(len_meanM); lenP.append(len_meanP)
        lenMstd.append(len_stdM); lenPstd.append(len_stdP)

np.save('dset_mono_boot.'+str(which_num),dsetM); np.save('dset_poly_boot.'+str(which_num),dsetP); 
np.save('MI_mono_boot.'+str(which_num),MI_m); np.save('MI_poly_boot.'+str(which_num),MI_p);
np.save('count_mono_boot.'+str(which_num),count_m); np.save('count_poly_boot.'+str(which_num),count_p);
np.save('sinAvg_mono_boot.'+str(which_num),sin_propM); np.save('sinAvg_poly_boot.'+str(which_num),sin_propP);
np.save('sinStd_mono_boot.'+str(which_num),sin_propMstd); np.save('sinStd_poly_boot.'+str(which_num),sin_propPstd);
np.save('len_mono_boot.'+str(which_num),lenM); np.save('len_poly_boot.'+str(which_num),lenP);
np.save('lenStd_mono_boot.'+str(which_num),lenMstd); np.save('lenStd_poly_boot.'+str(which_num),lenPstd);
