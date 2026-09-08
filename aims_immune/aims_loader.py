# ADD A CATCH TO SUPPRESS WARNINGS JUST FOR PIP VERSION OF THE APP
# catch here was for Biopython pairwise2. Remove that function and rely on AIMS_manuscripts for that instead.
#import warnings
#warnings.simplefilter("ignore")

from Bio import SeqIO
import numpy as np
import pandas

# More special stuff for the pip version of the script....
import aims_immune
# The -11 is because the filepath includes '__init__.py'
# So we need to remove that to get our data path.
datPath = aims_immune.__file__[:-11]

# Believe that "thing" here is trying to keep track of things for the GUI
# How the GUI keeps track of variables is... questionable at best...
thing = True

# For the msa_sub, this helps you get the subsets
def get_msa_sub(seqF,loc_start,loc_end):
    if len(loc_start) != len(loc_end):
        print("ERROR: Don't have same number of start and end entries")
        return()
    all_feat = []
    for i in np.transpose(seqF).values:
        pre_feat = []
        for j in np.arange(len(loc_start)):    
            s1 = loc_start[j]
            s2 = loc_end[j]
            struct = i[0][s1:s2]
            pre_feat.append(struct)
        all_feat.append([pre_feat])

    seqNEW = np.transpose(pandas.DataFrame(np.array(all_feat).reshape(len(all_feat),len(loc_start))))
    seqNEW.columns = seqF.columns
    return(seqNEW)

################################################
# Brand new loading function! No longer need to worry about formatting or chosing which "molecule type"
# The function really should just be smart enough to figure it out for you!
################################################
def seq_loader(seqPath,label,drop_dups=False,return_index=False,dataLoc=[],datType='',
               subset=False,subset_starts=[],subset_ends=[]):
    if datType.lower()!='fasta' or datType.lower()!='csv':
        with open(seqPath, 'r', encoding='utf-8') as file:
            first_line = file.readline()
        if first_line.find('>')!=-1:
            datType = 'fasta'
        else:
            datType = 'csv'
    ################################################
    # Code for loading in msa (mhc or other)
    ################################################
    if datType.lower()=='fasta':
        a = 0
        fin_seq = []; pre_id = []; fin_title=[]
        for seq_record in SeqIO.parse(seqPath,'fasta'):
            seqV=str(seq_record.seq)
            fasta_id = str(seq_record.description)
            # Keep this as a sort of index to match up
            titleV = label + '_' + str(a)
            fin_seq.append(seqV)
            pre_id.append(fasta_id)
            fin_title.append(titleV)
            a+=1

        seq_pre = pandas.DataFrame([fin_seq])
        seq_pre.columns = fin_title

        if subset:
            # basically need to allow for either version of subsetting
            # the newer version is CLEARLY better
            if len(subset_ends) == 0 and len(subset_starts)== 0:
                print('ERROR: did not define subset start/end!')
                print('Set subset=False if you want full sequence')
                return()
            elif len(subset_ends)==0:
                startF = []; subset_ends = []
                # Lotta ways to mess this up so probably need to discuss in methods
                for x in np.arange(len(subset_starts)-1):
                    startF = startF + [subset_starts[x]]
                    subset_ends = subset_ends + [subset_starts[x+1]]
            else:
                startF = subset_starts

            pre_fin_out = get_msa_sub(seq_pre,startF,subset_ends)
            if len(pre_fin_out)==0:
                print('Poorly defined subset bounds. Double check for errors')
                return()
        else:
            pre_fin_out = seq_pre

        # Think that the orientation of this matrix is not correct
        if drop_dups:
            fin_out = np.transpose(np.transpose(pre_fin_out).drop_duplicates())
            pre_id_df = np.transpose(pandas.DataFrame(pre_id))
            pre_id_df.columns = fin_title
            fin_id = pre_id_df[fin_out.columns]
        else:
            fin_out = pre_fin_out
            fin_id = pre_id

        if return_index:
            return(fin_out,fin_id)
        else:
            return(fin_out)

    ################################################
    # Code for loading in csvs (peptides or Igs)
    ################################################
    if datType.lower()=='csv':
        # We do, here, assume that we are 1. working with a csv and 2. that there is a header
        # think we have to ALWAYS assume that.
        tempAbs = pandas.read_csv(seqPath,sep=',',header=0)
        # Try to be a little more flexible about the data inputs we are accepting.
        if len(dataLoc) != 0:
            total_Abs = tempAbs[dataLoc]
        else:
            total_Abs = tempAbs

        # Remove empty entries
        total_abs1 = total_Abs.where((pandas.notnull(total_Abs)), '')

        # Remove X's in sequences... Should actually get a count of these at some point...
        totalF = total_abs1[~total_abs1.stack().str.contains('X').unstack().any(axis=1)]
        # Remove incomplete entries (i.e. missing cdr loops)
        final_Ig = totalF[~(totalF=='').any(axis=1)]

        # Remove degeneracies in the dataset (optional)
        if drop_dups:
            f_Ig = final_Ig.drop_duplicates()
        else:
            f_Ig = final_Ig

        final_title = [label + '_' + str(a) for a in np.arange(len(f_Ig))]
        final_Df = np.transpose(f_Ig)
        final_Df.columns = final_title

        if return_index:
            return(final_Df,f_Ig.index)
        else:
            return(final_Df)

def convert_3Let(inp):
    first = True
    three_let = ['ALA','GLY','ARG','LYS','ASP','GLU','ASN','GLN','MET','CYS','PHE','TYR','THR','TRP','PRO','SER','LEU','VAL','HIS','ILE']
    sin_let = [  'A',  'G',  'R',  'K',  'D',  'E',  'N',  'Q',  'M',  'C',  'F',  'Y',  'T',  'W',  'P',  'S',  'L',  'V',  'H',  'I']
    for i in inp:
        for scan in np.arange(len(three_let)):
            if i.lower() == three_let[scan].lower():
                hold = sin_let[scan]
                break
        if first:
            sin_final = hold
            first = False
        else:
            sin_final = np.hstack((sin_final,hold))
    return(sin_final)