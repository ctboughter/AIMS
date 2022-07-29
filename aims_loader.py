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

def mhc_loader(fastapath,mat_coords,label):

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
    if np.shape(np.shape(ff_seg1))[0] != 1:
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
        data = pandas.read_csv(fastapath,sep=',',header=0)['sequence']
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

# Special Loaders for MHC Germline Analysis
####################################################################################
# A Temporary Change Here: Load in HLA A, B, C platform domains for analysis
def get_HLA():
    from Bio import SeqIO
    mhc_dir = 'germline_data/'
    records = list(SeqIO.parse(mhc_dir+'ABC_prot.fasta','fasta'))

    full_id = []; concat_seq = []
    for i in np.arange(len(records)):
        rec_desc = records[i].description
        seqs_pre = str(records[i].seq)
        hla_loc = rec_desc.find('*')
        partial_id = rec_desc[hla_loc-1:]
        id_end = partial_id.find(' ')
        # So only the first two entries actually mean anything.
        # We always want the "01", so remove any entries that don't have "01"
        # at the very end. Pretty sure that's how these are organized
        if partial_id[:id_end][-2:] != '01':
            continue
        # So it turns out that there are some truncated sequences that ONLY
        # express the platform domain or a bit of the alpha3 domain.
        # We can continue to include these in the database if we truncate ALL sequences
        # to be this size (which would be fine by me). Also auto-include full-seqs
        start_plat = seqs_pre.find('SHS')
        if start_plat == -1:
            if len(seqs_pre) == 365:
                start_plat = 25
            elif len(seqs_pre) == 362:
                # oddly, HLA-B*15:16:01 has a slightly truncated alpha3
                # But seems to have the same length platform
                start_plat = 25
            elif len(seqs_pre) == 366:
                # oddly, HLA-C*05:171:01 has a slightly extended alpha3
                # But seems to have the same length platform
                # Both these oddities start plat with "SHF" end "KTL"
                start_plat = 25
            elif len(seqs_pre) == 181:
                start_plat = 0
            elif len(seqs_pre) == 273:
                start_plat = 0
        end_plat = seqs_pre.find('QRT')
        if end_plat == -1:
            if start_plat == 25:
                end_plat = 203
            else:
                end_plat = start_plat + 178
        seqs = seqs_pre[start_plat:end_plat]
        full_id = full_id +  [partial_id[:id_end]]
        concat_seq = concat_seq + [seqs]
    sequence_frame = np.transpose(pandas.DataFrame(np.vstack((full_id,concat_seq))))
    # remove duplicate entries
    seqFramePre1 = sequence_frame.drop_duplicates(subset=[1])
    # remove the "N", which is "null" according to IMGT
    seqFrameF = seqFramePre1.loc[~seqFramePre1[0].str.contains("N")]
    return(seqFrameF)

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


# A Temporary Change Here: Load in Class II HLA-DQ, DR, and DP
def get_classII():
    mhc_dir = 'germline_data/'
    classIIa = list(SeqIO.parse(mhc_dir+'classII_alpha.fasta','fasta'))
    classIIb = list(SeqIO.parse(mhc_dir+'classII_beta.fasta','fasta'))

    concat_seq = []
    fin_vect = ['Code','Allele','BP','Drop']
    for i in np.arange(len(classIIa)):
        rec_desc = classIIa[i].description
        seqs_pre = str(classIIa[i].seq)
        pre_df = rec_desc.split(' ')
        
        # There are far fewer sequences, so let's align them to
        # a structural example where we identified the platfrom
        # PDB ID: 4P4K
        test1 = 'LYS ALA ASP HIS VAL SER THR TYR ALA ALA PHE VAL GLN THR HIS ARG PRO THR GLY GLU PHE MET PHE GLU PHE ASP GLU ASP GLU MET PHE TYR VAL ASP LEU ASP LYS LYS GLU THR VAL TRP HIS LEU GLU GLU PHE GLY GLN ALA PHE SER PHE GLU ALA GLN GLY GLY LEU ALA ASN ILE ALA ILE LEU ASN ASN ASN LEU ASN THR LEU ILE GLN ARG SER ASN HIS THR'
        test1_f = test1.split(' ')
        dpalpha_struct = ''.join(convert_3Let(test1_f))
        
        # This combo REALLY punishes opening up gaps, which is what we want.
        aa = pairwise2.align.localms(dpalpha_struct,seqs_pre,0.5, -0.1, -5,-0.5)
        anchor = aa[0][0]
        seqs_align = aa[0][1]
        start_hold = anchor.find('KADH')
        end_hold = anchor.find('RSNHT') + 5
        seqs_f = seqs_align[start_hold:end_hold]
        
        if seqs_f.find('-') == -1:
            fin_vect = np.vstack((fin_vect,pre_df))
            concat_seq = concat_seq + [seqs_f]
        
    fin_vect_alpha = fin_vect
    concat_seq_alpha = concat_seq

    concat_seq = []
    fin_vect = ['Code','Allele','BP','Drop']
    for j in np.arange(len(classIIb)):
        rec_desc = classIIb[j].description
        seqs_pre = str(classIIb[j].seq)
        pre_df = rec_desc.split(' ')
        # There are far fewer sequences, so let's align them to
        # a structural example where we identified the platfrom
        # PDB ID: 4P4K
        test2 = 'GLN ALA PHE TRP ILE ASP LEU PHE GLU THR ILE GLY SER PRO GLU ASN TYR LEU PHE GLN GLY ARG GLN GLU CYS TYR ALA PHE ASN GLY THR GLN ARG PHE LEU GLU ARG TYR ILE TYR ASN ARG GLU GLU PHE VAL ARG PHE ASP SER ASP VAL GLY GLU PHE ARG ALA VAL THR GLU LEU GLY ARG PRO ASP GLU GLU TYR TRP ASN SER GLN LYS ASP ILE LEU GLU GLU GLU ARG ALA VAL PRO ASP ARG MET CYS ARG HIS ASN TYR GLU LEU GLY GLY PRO MET THR LEU GLN'
        test2_f = test2.split(' ')
        dpbeta_struct = ''.join(convert_3Let(test2_f))
        
        # This combo REALLY punishes opening up gaps, which is what we want.
        aa = pairwise2.align.localms(dpbeta_struct,seqs_pre,0.5, -0.1, -5,-0.5)
        anchor = aa[0][0]
        seqs_align = aa[0][1]
        start_hold = anchor.find('QAFW')
        end_hold = anchor.find('PMTLQ') + 5
        seqs_f = seqs_align[start_hold:end_hold]
        
        if seqs_f.find('-') == -1:
            fin_vect = np.vstack((fin_vect,pre_df))
            concat_seq = concat_seq + [seqs_f]
        
    fin_vect_beta = fin_vect
    concat_seq_beta = concat_seq

    alpha_df_pre = pandas.DataFrame(fin_vect_alpha[1:])
    alpha_df_pre.columns = fin_vect_alpha[0]
    alpha_seq_temp = pandas.DataFrame(concat_seq_alpha)
    alpha_df = pandas.concat([alpha_df_pre,alpha_seq_temp],axis=1)
    # remove the alleles "N", which is "null" according to IMGT
    alpha_noN = alpha_df[~alpha_df['Allele'].str.contains("N")]
    alpha_df = alpha_noN.drop_duplicates(0)

    beta_df_pre = pandas.DataFrame(fin_vect_beta[1:])
    beta_df_pre.columns = fin_vect_beta[0]
    beta_seq_temp = pandas.DataFrame(concat_seq_beta)
    beta_df = pandas.concat([beta_df_pre,beta_seq_temp],axis=1)
    # remove the alleles "N", which is "null" according to IMGT
    beta_noN = beta_df[~beta_df['Allele'].str.contains("N")]
    beta_df = beta_noN.drop_duplicates(0)

    return(alpha_df,beta_df)