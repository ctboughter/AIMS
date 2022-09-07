from Bio import AlignIO
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
from Bio.pairwise2 import format_alignment
import numpy as np
import pandas

thing = True

# Using the nomenclature of the GUI explanation, here are some example GUI start/end values
# As a reminder, it goes S1s, S1e/H1s, H1e/S2s, S2e/H2s, H2e
# For the ji_cartFish we have: 2,49,93,152,193
# For the cd1d.fasta we have: 124,167,209,262,303
# For the hlaA.fasta we have: 170,218,260,306,348
# For cd1_ufa_genes.fasta: 22,66,105,158,199

def mhc_loader(fastapath,mat_coords,label,drop_dups = False):

    thing = True
    xxx1 = fastapath.rfind('/')
    xxx2 = fastapath.rfind('.fasta')
    yyy = fastapath[xxx1+1:xxx2]
    a = 0
    for seq_record in SeqIO.parse(fastapath,'fasta'):
        seqV=str(seq_record.seq)
        fasta_id = str(seq_record.id)
        ori_titleV = yyy + ' - ' + fasta_id
        # Replace 
        titleV = label + '_' + str(a)
        
        seg1 = seqV[int(mat_coords[0]):int(mat_coords[1])].replace('-','')
        seg2 = seqV[int(mat_coords[1]):int(mat_coords[2])].replace('-','')
        seg3 = seqV[int(mat_coords[2]):int(mat_coords[3])].replace('-','')
        seg4 = seqV[int(mat_coords[3]):int(mat_coords[4])].replace('-','')
        
        segs = [seg1,seg2,seg3,seg4]
        if thing:
            final_Seg1 = segs
            final_title = [titleV]
            final_ori_title = [ori_titleV]
            thing = False
        else:
            final_Seg1 = np.vstack((final_Seg1,segs))
            final_title = final_title + [titleV]
            final_ori_title = final_ori_title + [ori_titleV]
        
        a = a+1
    
    ff_seg1 = np.transpose(final_Seg1)
    # Obviously don't need to worry about extra sequences if there is only one...
    # What a dumb f***ing way to do this, but it works...
    if np.shape(np.shape(ff_seg1))[0] == 1:
        drop_dups = False

    if drop_dups:
        # Extra bit here to delete duplicate sequences
        # REMOVE POINT MUTANTS AND SEQs TWO MUTATIONS OFF
        num_muts = 2
        aa,bb = np.shape(ff_seg1)
        indices = np.array([0,0])
        for i in np.arange(bb):
            for j in np.arange(bb):
                if i == j:
                    continue
                count = 0
                for k in np.arange(aa):
                    # SO THIS IS A NICE CODE I STOLE FROM ONLINE TO FIND NUMBER OF MATCHES
                    # ABSOLUTE VALUE COUNTS DIFF LENGTH AS A MISMATCH
                    count += sum(1 for a, b in zip(ff_seg1[k,i], ff_seg1[k,j]) if a != b) + abs(len(ff_seg1[k,i]) - len(ff_seg1[k,j]))
                if count < num_muts:        
                    indices = np.vstack((indices,[i,j]))

        thing = True
        for i in np.arange(len(indices)):
            if len(indices) < 3:
                break
            if indices[i,0] < indices[i,1]:
                if thing:
                    index_new = [indices[i,0]]
                    thing = False
                elif len(index_new) == 1:
                    if index_new == indices[i,0]:
                        continue
                    else:
                        index_new = np.vstack((index_new,indices[i,0]))
                elif len(index_new) > 1:
                    if index_new[len(index_new)-1] == indices[i,0]:
                        continue
                    else:
                        index_new = np.vstack((index_new,indices[i,0]))

        if len(indices) < 3:
            finalDF = pandas.DataFrame(ff_seg1,columns = final_title)
            title_key = np.vstack((final_title,final_ori_title))
            return(finalDF,title_key)
        else:
            seq_new = np.delete(ff_seg1,index_new,axis = 1)
            title_new = np.delete(final_title,index_new,axis = 0)
            title_ori_new = np.delete(final_ori_title,index_new,axis = 0)
            finalDF = pandas.DataFrame(seq_new,columns = title_new)
            title_key = np.vstack((title_new,title_ori_new))
            return(finalDF,title_key)
    else:
        seq_new = ff_seg1
        title_new = final_title
        title_ori_new = final_ori_title
        finalDF = pandas.DataFrame(seq_new,columns = title_new)
        title_key = np.vstack((title_new,title_ori_new))
        return(finalDF,title_key)
# So in the main version of the script, we have a special loader for each data subset
# Can we make just a generalizable one? Let's give it a try...
def Ig_loader(fastapath,label,loops=6,drop_degens = False):
    if loops == 6:
        total_Abs=pandas.read_csv(fastapath,sep=',',header=0,names=['cdrL1_aa','cdrL2_aa','cdrL3_aa','cdrH1_aa','cdrH2_aa','cdrH3_aa'])
    elif loops == 3:
        total_Abs=pandas.read_csv(fastapath,sep=',',header=0,names=['cdr1_aa','cdr2_aa','cdr3_aa'])
    elif loops == 2:
        total_Abs=pandas.read_csv(fastapath,sep=',',header=0,names=['cdrH3_aa','cdrL3_aa'])
    elif loops == 1:
        total_Abs=pandas.read_csv(fastapath,sep=',',header=0,names=['cdr_aa'])
    # Remove empty entries
    total_abs1 = total_Abs.where((pandas.notnull(total_Abs)), '')
    # Remove X's in sequences... Should actually get a count of these at some point...
    if loops == 6:
        total_abs2=total_abs1[~total_abs1['cdrL1_aa'].str.contains("X")]
        total_abs3=total_abs2[~total_abs2['cdrL2_aa'].str.contains("X")]
        total_abs4=total_abs3[~total_abs3['cdrL3_aa'].str.contains("X")]
        total_abs5=total_abs4[~total_abs4['cdrH1_aa'].str.contains("X")]
        total_abs6=total_abs5[~total_abs5['cdrH2_aa'].str.contains("X")]
        totalF=total_abs6[~total_abs6['cdrH3_aa'].str.contains("X")].values
    elif loops == 3:
        total_abs5=total_abs1[~total_abs1['cdr1_aa'].str.contains("X")]
        total_abs6=total_abs5[~total_abs5['cdr2_aa'].str.contains("X")]
        totalF=total_abs6[~total_abs6['cdr3_aa'].str.contains("X")].values
    elif loops == 2:
        total_abs5=total_abs1[~total_abs1['cdrH3_aa'].str.contains("X")]
        totalF=total_abs5[~total_abs5['cdrL3_aa'].str.contains("X")].values
    elif loops == 1:
        totalF=total_abs1[~total_abs1['cdr_aa'].str.contains("X")].values
    # Remove incomplete entries
    a=0
    del_these=[]
    if loops == 6:
        for i in np.arange(len(totalF[:,5])):
            if totalF[i,5] == '' or totalF[i,4] == '' or totalF[i,3] == '' or totalF[i,2] == '' or totalF[i,1] == '' or totalF[i,0] == '':
                if a == 0:
                    del_these=i
                else:
                    del_these=np.vstack((del_these,i))
                a=a+1
    elif loops == 3:
        for i in np.arange(len(totalF[:,2])):
            if totalF[i,2] == '' or totalF[i,1] == '' or totalF[i,0] == '':
                if a == 0:
                    del_these=i
                else:
                    del_these=np.vstack((del_these,i))
                a=a+1
    elif loops == 2:
        for i in np.arange(np.shape(totalF)[0]):
            if totalF[i,1] == '' or totalF[i,0] == '':
                if a == 0:
                    del_these=i
                else:
                    del_these=np.vstack((del_these,i))
                a=a+1
    elif loops == 1:
        for i in np.arange(len(totalF[:])):
            if totalF[i] == '':
                if a == 0:
                    del_these=i
                else:
                    del_these=np.vstack((del_these,i))
                a=a+1

    final_Ig=np.delete(totalF,del_these,axis=0)

    # Remove degeneracies in the dataset (optional)
    if drop_degens:
        aa = np.shape(final_Ig)[0]
        for i in np.arange(aa):
            degen = False
            for j in np.arange(i):
                # ignore diagonal
                if i == j:
                    continue
                # to get around changing number of loops,
                if loops == 1:
                    test1 = final_Ig[i,0]
                    test2 = final_Ig[j,0]
                elif loops == 2:
                    test1 = final_Ig[i,0] + final_Ig[i,1]
                    test2 = final_Ig[j,0] + final_Ig[j,1]
                elif loops == 3:
                    test1 = final_Ig[i,0] + final_Ig[i,1] + final_Ig[i,2]
                    test2 = final_Ig[j,0] + final_Ig[j,1] + final_Ig[j,2]
                elif loops == 6:
                    test1 = final_Ig[i,0] + final_Ig[i,1] + final_Ig[i,2] + final_Ig[i,3] + final_Ig[i,4] + final_Ig[i,5]
                    test2 = final_Ig[j,0] + final_Ig[j,1] + final_Ig[j,2] + final_Ig[j,3] + final_Ig[j,4] + final_Ig[j,5]
                # if the sequences are of a different length, clearly they aren't identical
                if len(test1) - len(test2) != 0:
                    continue
                # Sum zip here counts the number of matched residues (position senstive)
                # So by subtracting the length, identical sequences should have count = 0
                count = sum(1 for a, b in zip(test1, test2) if a == b) - len(test1)
                # as soon as you find an identical sequence, break out
                if count == 0:
                    degen = True
                    break
            if i == 0 and not degen:
                indices = np.array([0])
            elif not degen:        
                indices = np.vstack((indices,i))
        if loops == 1:
            f_Ig = final_Ig[indices,:].reshape(len(indices),1)
        elif loops == 2:
            f_Ig = final_Ig[indices,:].reshape(len(indices),2)
        elif loops == 3:
            f_Ig = final_Ig[indices,:].reshape(len(indices),3)
        elif loops == 6:
            f_Ig = final_Ig[indices,:].reshape(len(indices),6)
    else:
        f_Ig = final_Ig

    final_title = [label + '_' + str(a) for a in np.arange(len(f_Ig))]
    final_Df = pandas.DataFrame(np.transpose(f_Ig),columns = final_title)

    return(final_Df)
#####################################################################################
def pep_loader(fastapath,label, scrape=False, start_label=0):

    thing = True
    xxx1 = fastapath.rfind('/')
    xxx2 = fastapath.rfind('.csv')
    yyy = fastapath[xxx1+1:xxx2]
    a = 0
    # Alright so for now there are abolsutely no standards here
    if scrape:
        csv_file = pandas.read_csv(fastapath,sep=',',header=0)
        # Need to do this because not every file calls their "MHC class"
        # the same thing. Some must predict, some must control for it...
        headers = csv_file.columns
        data = csv_file['search_hit']
        # I believe the MHC allele is ALWAYS the last column,
        # but I should probably make sure of that at some point
        allele = csv_file[headers[-1]]
    else:
        data = pandas.read_csv(fastapath,sep=',',header=1)['sequence']
    for i in np.arange(len(data)):
        # Replace 
        titleV = label + '_' + str(a+start_label)
        
        if thing:
            final_title = [titleV]
            thing = False
        else:
            final_title = final_title + [titleV]
        a = a+1

    finalDF = np.transpose(pandas.DataFrame(np.array(data)))
    finalDF.columns=final_title

    if scrape:
        finalAllele = np.transpose(pandas.DataFrame(np.array(allele)))
        finalAllele.columns=final_title
        return(finalDF,finalAllele)
    else:
        return(finalDF)
