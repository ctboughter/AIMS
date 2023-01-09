### This script is for loading VERY specific antibody data from the eLife manuscript ###
# To see how more general data is loaded in, look through "aims_loader.py"
# Let's start off by loading in Jeff's CDR3's
import numpy as np
import pandas

def getBunker():
    total_Abs=pandas.read_csv('app_data/mouse_IgA.dat',sep='\s+',header=None,names=['cdrL1_aa','cdrL2_aa','cdrL3_aa','cdrH1_aa','cdrH2_aa','cdrH3_aa','react'])
    total_abs1 = total_Abs.where((pandas.notnull(total_Abs)), '')
    # Remove X's in sequences... Should actually get a count of these at some point...
    total_abs2=total_abs1[~total_abs1['cdrL1_aa'].str.contains("X")]
    total_abs3=total_abs2[~total_abs2['cdrL2_aa'].str.contains("X")]
    total_abs4=total_abs3[~total_abs3['cdrL3_aa'].str.contains("X")]
    total_abs5=total_abs4[~total_abs4['cdrH1_aa'].str.contains("X")]
    total_abs6=total_abs5[~total_abs5['cdrH2_aa'].str.contains("X")]
    total_abs7=total_abs6[~total_abs6['cdrH3_aa'].str.contains("X")]

    mono_all=total_abs7[total_abs7['react'].isin([0.0,1.0])].values
    poly_all=total_abs7[total_abs7['react'].isin([2.0,3.0,4.0,5.0,6.0,7.0])].values

    mono=total_abs7[total_abs7['react'].isin([0.0])].values
    poly=total_abs7[total_abs7['react'].isin([5.0,6.0,7.0])].values

    a=0
    del_these=[]
    for i in np.arange(len(mono_all[:,5])):
        if mono_all[i,5] == '' or mono_all[i,4] == '' or mono_all[i,3] == '' or mono_all[i,2] == '' or mono_all[i,1] == '' or mono_all[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    mono_all2=np.delete(mono_all,del_these,axis=0)

    a=0
    del_these=[]
    for i in np.arange(len(poly_all[:,5])):
        if poly_all[i,5] == '' or poly_all[i,4] == '' or poly_all[i,3] == '' or poly_all[i,2] == '' or poly_all[i,1] == '' or poly_all[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    poly_all2=np.delete(poly_all,del_these,axis=0)
    
    a=0
    del_these=[]
    for i in np.arange(len(mono[:,5])):
        if mono[i,5] == '' or mono[i,4] == '' or mono[i,3] == '' or mono[i,2] == '' or mono[i,1] == '' or mono[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    mono2=np.delete(mono,del_these,axis=0)

    a=0
    del_these=[]
    for i in np.arange(len(poly[:,5])):
        if poly[i,5] == '' or poly[i,4] == '' or poly[i,3] == '' or poly[i,2] == '' or poly[i,1] == '' or poly[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    poly2=np.delete(poly,del_these,axis=0)

    return(np.transpose(mono_all2[:,0:6]),np.transpose(poly_all2[:,0:6]),np.transpose(mono2[:,0:6]),np.transpose(poly2[:,0:6]))
#####################################################################################

def getJenna():
    total_Abs=pandas.read_csv('app_data/flu_IgG.dat',sep='\s+',header=None,
    names=['cdrL1_aa','cdrL2_aa','cdrL3_aa','cdrH1_aa','cdrH2_aa','cdrH3_aa','react'])

    total_abs1 = total_Abs.where((pandas.notnull(total_Abs)), '')
    # Remove X's in sequences... Should actually get a count of these at some point...
    total_abs2=total_abs1[~total_abs1['cdrL1_aa'].str.contains("X")]
    total_abs3=total_abs2[~total_abs2['cdrL2_aa'].str.contains("X")]
    total_abs4=total_abs3[~total_abs3['cdrL3_aa'].str.contains("X")]
    total_abs5=total_abs4[~total_abs4['cdrH1_aa'].str.contains("X")]
    total_abs6=total_abs5[~total_abs5['cdrH2_aa'].str.contains("X")]
    total_abs7=total_abs6[~total_abs6['cdrH3_aa'].str.contains("X")]

    # Having this and the above lines as "if" options could make this loader more generalizable...
    mono_all=total_abs7[total_abs7['react'].isin([0,1])].values
    poly_all=total_abs7[total_abs7['react'].isin([2,3,4,5,6,7])].values

    mono=total_abs7[total_abs7['react'].isin([0])].values
    poly=total_abs7[total_abs7['react'].isin([5,6,7])].values

    a=0
    del_these=[]
    for i in np.arange(len(mono_all[:,5])):
        if mono_all[i,5] == '' or mono_all[i,4] == '' or mono_all[i,3] == '' or mono_all[i,2] == '' or mono_all[i,1] == '' or mono_all[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    mono_all2=np.delete(mono_all,del_these,axis=0)

    a=0
    del_these=[]
    for i in np.arange(len(poly_all[:,5])):
        if poly_all[i,5] == '' or poly_all[i,4] == '' or poly_all[i,3] == '' or poly_all[i,2] == '' or poly_all[i,1] == '' or poly_all[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    poly_all2=np.delete(poly_all,del_these,axis=0)

    a=0
    del_these=[]
    for i in np.arange(len(mono[:,5])):
        if mono[i,5] == '' or mono[i,4] == '' or mono[i,3] == '' or mono[i,2] == '' or mono[i,1] == '' or mono[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    mono2=np.delete(mono,del_these,axis=0)

    a=0
    del_these=[]
    for i in np.arange(len(poly[:,5])):
        if poly[i,5] == '' or poly[i,4] == '' or poly[i,3] == '' or poly[i,2] == '' or poly[i,1] == '' or poly[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    poly2=np.delete(poly,del_these,axis=0)

    return(np.transpose(mono_all2[:,0:6]),np.transpose(poly_all2[:,0:6]),np.transpose(mono2[:,0:6]),np.transpose(poly2[:,0:6]))

def getHugo():
    my_heavy=pandas.read_csv('app_data/hiv_igg_data/gut_heavy_aa.dat',sep='\s+')
    my_light=pandas.read_csv('app_data/hiv_igg_data/gut_light_aa.dat',sep='\s+')
    poly_YN=pandas.read_csv('app_data/hiv_igg_data/gut_num_react.dat',sep='\s+',header=None,names=['react'])
    total_abs=pandas.concat([my_light,my_heavy,poly_YN],axis=1)
    total_abs7 = total_abs.where((pandas.notnull(total_abs)), '')
    mono_all=total_abs7[total_abs7['react'].isin([0,1])].values
    poly_all=total_abs7[total_abs7['react'].isin([2,3,4])].values
    mono=total_abs7[total_abs7['react'].isin([0])].values
    poly=total_abs7[total_abs7['react'].isin([3,4])].values

    a=0
    del_these=[]
    for i in np.arange(len(mono_all[:,5])):
        if mono_all[i,5] == '' or mono_all[i,4] == '' or mono_all[i,3] == '' or mono_all[i,2] == '' or mono_all[i,1] == '' or mono_all[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    mono_all2=np.delete(mono_all,del_these,axis=0)

    a=0
    del_these=[]
    for i in np.arange(len(poly_all[:,5])):
        if poly_all[i,5] == '' or poly_all[i,4] == '' or poly_all[i,3] == '' or poly_all[i,2] == '' or poly_all[i,1] == '' or poly_all[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    poly_all2=np.delete(poly_all,del_these,axis=0)

    a=0
    del_these=[]
    for i in np.arange(len(mono[:,5])):
        if mono[i,5] == '' or mono[i,4] == '' or mono[i,3] == '' or mono[i,2] == '' or mono[i,1] == '' or mono[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    mono2=np.delete(mono,del_these,axis=0)

    a=0
    del_these=[]
    for i in np.arange(len(poly[:,5])):
        if poly[i,5] == '' or poly[i,4] == '' or poly[i,3] == '' or poly[i,2] == '' or poly[i,1] == '' or poly[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    poly2=np.delete(poly,del_these,axis=0)
    return(np.transpose(mono_all2[:,0:6]),np.transpose(poly_all2[:,0:6]),np.transpose(mono2[:,0:6]),np.transpose(poly2[:,0:6]))

def getHugo_Nature():
    my_heavy=pandas.read_csv('app_data/hiv_igg_data/nat_heavy_aa.dat',sep='\s+')
    my_light=pandas.read_csv('app_data/hiv_igg_data/nat_light_aa.dat',sep='\s+')
    poly_YN=pandas.read_csv('app_data/hiv_igg_data/nat_num_react.dat',sep='\s+',header=None,names=['react'])
    total_Abs=pandas.concat([my_light,my_heavy,poly_YN],axis=1)

    total_abs1 = total_Abs.where((pandas.notnull(total_Abs)), '')
    # Remove X's in sequences... Should actually get a count of these at some point...
    total_abs2=total_abs1[~total_abs1['cdrL1_aa'].str.contains("X")]
    total_abs3=total_abs2[~total_abs2['cdrL2_aa'].str.contains("X")]
    total_abs4=total_abs3[~total_abs3['cdrL3_aa'].str.contains("X")]
    total_abs5=total_abs4[~total_abs4['cdrH1_aa'].str.contains("X")]
    total_abs6=total_abs5[~total_abs5['cdrH2_aa'].str.contains("X")]
    total_abs7=total_abs6[~total_abs6['cdrH3_aa'].str.contains("X")]

    # And finish it up...
    mono_all=total_abs7[total_abs7['react'].isin([0.0,1.0])].values
    poly_all=total_abs7[total_abs7['react'].isin([2.0,3.0,4.0,5.0,6.0])].values

    mono=total_abs7[total_abs7['react'].isin([0.0])].values
    poly=total_abs7[total_abs7['react'].isin([5.0,6.0])].values

    a=0
    del_these=[]
    for i in np.arange(len(mono_all[:,5])):
        if mono_all[i,5] == '' or mono_all[i,4] == '' or mono_all[i,3] == '' or mono_all[i,2] == '' or mono_all[i,1] == '' or mono_all[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    mono_all2=np.delete(mono_all,del_these,axis=0)

    a=0
    del_these=[]
    for i in np.arange(len(poly_all[:,5])):
        if poly_all[i,5] == '' or poly_all[i,4] == '' or poly_all[i,3] == '' or poly_all[i,2] == '' or poly_all[i,1] == '' or poly_all[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    poly_all2=np.delete(poly_all,del_these,axis=0)
    
    a=0
    del_these=[]
    for i in np.arange(len(mono[:,5])):
        if mono[i,5] == '' or mono[i,4] == '' or mono[i,3] == '' or mono[i,2] == '' or mono[i,1] == '' or mono[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    mono2=np.delete(mono,del_these,axis=0)

    a=0
    del_these=[]
    for i in np.arange(len(poly[:,5])):
        if poly[i,5] == '' or poly[i,4] == '' or poly[i,3] == '' or poly[i,2] == '' or poly[i,1] == '' or poly[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    poly2=np.delete(poly,del_these,axis=0)

    return(np.transpose(mono_all2[:,0:6]),np.transpose(poly_all2[:,0:6]),np.transpose(mono2[:,0:6]),np.transpose(poly2[:,0:6]))

def getHugo_NatCNTRL():
    my_heavy=pandas.read_csv('app_data/hiv_igg_data/nat_cntrl_heavy_aa.dat',sep='\s+')
    my_light=pandas.read_csv('app_data/hiv_igg_data/nat_cntrl_light_aa.dat',sep='\s+')
    poly_YN=pandas.read_csv('app_data/hiv_igg_data/nat_cntrl_num_react.dat',sep='\s+',header=None,names=['react'])
    total_Abs=pandas.concat([my_light,my_heavy,poly_YN],axis=1)

    total_abs1 = total_Abs.where((pandas.notnull(total_Abs)), '')
    # Remove X's in sequences... Should actually get a count of these at some point...
    total_abs2=total_abs1[~total_abs1['cdrL1_aa'].str.contains("X")]
    total_abs3=total_abs2[~total_abs2['cdrL2_aa'].str.contains("X")]
    total_abs4=total_abs3[~total_abs3['cdrL3_aa'].str.contains("X")]
    total_abs5=total_abs4[~total_abs4['cdrH1_aa'].str.contains("X")]
    total_abs6=total_abs5[~total_abs5['cdrH2_aa'].str.contains("X")]
    total_abs7=total_abs6[~total_abs6['cdrH3_aa'].str.contains("X")]

    # And finish it up...
    mono_all=total_abs7[total_abs7['react'].isin([0.0,1.0])].values
    poly_all=total_abs7[total_abs7['react'].isin([2.0,3.0,4.0,5.0,6.0])].values

    mono=total_abs7[total_abs7['react'].isin([0.0])].values
    poly=total_abs7[total_abs7['react'].isin([5.0,6.0])].values

    a=0
    del_these=[]
    for i in np.arange(len(mono_all[:,5])):
        if mono_all[i,5] == '' or mono_all[i,4] == '' or mono_all[i,3] == '' or mono_all[i,2] == '' or mono_all[i,1] == '' or mono_all[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    mono_all2=np.delete(mono_all,del_these,axis=0)

    a=0
    del_these=[]
    for i in np.arange(len(poly_all[:,5])):
        if poly_all[i,5] == '' or poly_all[i,4] == '' or poly_all[i,3] == '' or poly_all[i,2] == '' or poly_all[i,1] == '' or poly_all[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    poly_all2=np.delete(poly_all,del_these,axis=0)
    
    a=0
    del_these=[]
    for i in np.arange(len(mono[:,5])):
        if mono[i,5] == '' or mono[i,4] == '' or mono[i,3] == '' or mono[i,2] == '' or mono[i,1] == '' or mono[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    mono2=np.delete(mono,del_these,axis=0)

    a=0
    del_these=[]
    for i in np.arange(len(poly[:,5])):
        if poly[i,5] == '' or poly[i,4] == '' or poly[i,3] == '' or poly[i,2] == '' or poly[i,1] == '' or poly[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    poly2=np.delete(poly,del_these,axis=0)

    return(np.transpose(mono_all2[:,0:6]),np.transpose(poly_all2[:,0:6]),np.transpose(mono2[:,0:6]),np.transpose(poly2[:,0:6]))


def getHugo_PLOS():
    my_heavy=pandas.read_csv('app_data/hiv_igg_data/plos_heavy_aa.dat',sep='\s+')
    my_light=pandas.read_csv('app_data/hiv_igg_data/plos_light_aa.dat',sep='\s+')
    poly_YN=pandas.read_csv('app_data/hiv_igg_data/plos_yn.dat',sep='\s+',header=None,names=['YN'])
    total_Abs=pandas.concat([my_light,my_heavy,poly_YN],axis=1)

    total_abs1 = total_Abs.where((pandas.notnull(total_Abs)), '')
    # Remove X's in sequences... Should actually get a count of these at some point...
    total_abs2=total_abs1[~total_abs1['cdrL1_aa'].str.contains("X")]
    total_abs3=total_abs2[~total_abs2['cdrL2_aa'].str.contains("X")]
    total_abs4=total_abs3[~total_abs3['cdrL3_aa'].str.contains("X")]
    total_abs5=total_abs4[~total_abs4['cdrH1_aa'].str.contains("X")]
    total_abs6=total_abs5[~total_abs5['cdrH2_aa'].str.contains("X")]
    total_abs7=total_abs6[~total_abs6['cdrH3_aa'].str.contains("X")]

    # And finish it up...
    mono_all=total_abs7[total_abs7['YN']=='N'].values
    poly_all=total_abs7[total_abs7['YN']=='Y'].values

    a=0
    del_these=[]
    for i in np.arange(len(mono_all[:,5])):
        if mono_all[i,5] == '' or mono_all[i,4] == '' or mono_all[i,3] == '' or mono_all[i,2] == '' or mono_all[i,1] == '' or mono_all[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    mono_all2=np.delete(mono_all,del_these,axis=0)

    a=0
    del_these=[]
    for i in np.arange(len(poly_all[:,5])):
        if poly_all[i,5] == '' or poly_all[i,4] == '' or poly_all[i,3] == '' or poly_all[i,2] == '' or poly_all[i,1] == '' or poly_all[i,0] == '':
            if a == 0:
                del_these=i
            else:
                del_these=np.vstack((del_these,i))
            a=a+1
    poly_all2=np.delete(poly_all,del_these,axis=0)
    return(np.transpose(mono_all2[:,0:6]),np.transpose(poly_all2[:,0:6]))

def getAdimab():
    heavy_Abs=pandas.read_csv('app_data/adimab_data/cdrs_H_final.txt',sep='\s+',header=None,names=['cdrH1_aa','cdrH2_aa','cdrH3_aa'])
    light_Abs=pandas.read_csv('app_data/adimab_data/cdrs_L_final.txt',sep='\s+',header=None,names=['cdrL1_aa','cdrL2_aa','cdrL3_aa'])
    outcomes=pandas.read_csv('app_data/adimab_data/drug_outcomes.csv',sep=',',header=0)
    assays=pandas.read_csv('app_data/adimab_data/drug_properties.csv',sep=',',header=0)
    
    names=outcomes['Name']
    clinical=outcomes['Clinical Status']
    phage=outcomes['Phagec']
    elisa_polyScores=assays['ELISA']
    psr_assayScore=assays['Poly-Specificity Reagent (PSR) SMP Score (0-1)']

    total_Abs=pandas.concat([names,heavy_Abs,clinical,light_Abs,phage,assays.loc[:, assays.columns != 'Unnamed: 13']],axis=1).dropna()

    # Let's not process this data, just return the matrix
    return(total_Abs)
#####################################################################################

def getSabDab():
    heavy_Abs=pandas.read_csv('app_data/SabDab_data/nonAdimab_igblast/cdrs_H_final.txt',sep='\s+',header=None,names=['cdrH1_aa','cdrH2_aa','cdrH3_aa'])
    light_Abs=pandas.read_csv('app_data/SabDab_data/nonAdimab_igblast/cdrs_L_final.txt',sep='\s+',header=None,names=['cdrL1_aa','cdrL2_aa','cdrL3_aa'])
    notAdi=pandas.read_csv('app_data/SabDab_data/non_adimab_dataHu.csv',sep=',',header=0)
    Adi=pandas.read_csv('app_data/SabDab_data/adimab_SabDabdata.csv',sep=',',header=0)
    
    adiName=Adi['Therapeutic']; NotadiName=notAdi['Therapeutic']
    adiOutcome=Adi["Highest_Clin_Trial (Jan '20)"]
    NotadiOutcome=notAdi["Highest_Clin_Trial (Jan '20)"]
    adiDeact=Adi['Est. Status']; NotadiDeact=notAdi['Est. Status']

    adimab_info=pandas.concat([adiName,adiOutcome,adiDeact],axis=1).dropna()
    notadimab_all=pandas.concat([NotadiName,heavy_Abs,light_Abs,NotadiOutcome,NotadiDeact],axis=1).dropna()
    # Let's not process this data, just return the matrix
    return(adimab_info,notadimab_all)
#####################################################################################
