# ADD A CATCH TO SUPPRESS WARNINGS JUST FOR PIP VERSION OF THE APP
import warnings
warnings.simplefilter("ignore")

from Bio import AlignIO
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import pairwise2
import numpy as np
import pandas

# More special stuff for the pip version of the script....
import aims_immune
# The -11 is because the filepath includes '__init__.py'
# So we need to remove that to get our data path.
datPath = aims_immune.__file__[:-11]

thing = True

def msa_loader(fastapath,label,drop_dups = False):
    thing = True; a = 0
    fin_seq = []; fin_id = []; fin_title=[]
    for seq_record in SeqIO.parse(fastapath,'fasta'):
        seqV=str(seq_record.seq)
        fasta_id = str(seq_record.description)
        # Keep this as a sort of index to match up
        titleV = label + '_' + str(a)
        fin_seq.append(seqV)
        fin_id.append(fasta_id)
        fin_title.append(titleV)
        a+=1

    fin_out = pandas.DataFrame([fin_seq,fin_id])
    fin_out.columns = fin_title
    return(fin_out)

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
def Ig_loader(fastapath,label,loops=6,drop_degens = False,return_index = False):
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
    if return_index:
        index_temp = total_abs1.index
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
    if return_index:
        index_temp_del=np.delete(index_temp,del_these,axis=0)

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
        if return_index:
            index_f = np.array(index_temp_del)[indices]

    else:
        f_Ig = final_Ig
        if return_index:
            index_f = index_temp_del

    final_title = [label + '_' + str(a) for a in np.arange(len(f_Ig))]
    final_Df = pandas.DataFrame(np.transpose(f_Ig),columns = final_title)

    if return_index:
        return(final_Df,index_f)
    else:
        return(final_Df)
#####################################################################################
def pep_loader(fastapath,label, scrape=False, start_label=0,drop_degens=False,len_cutoff = 12,return_index=False):

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
        dataF = csv_file['search_hit']
        # I believe the MHC allele is ALWAYS the last column,
        # but I should probably make sure of that at some point
        allele = csv_file[headers[-1]]
    else:
        data_pre = pandas.read_csv(fastapath,sep=',',header=0)['sequence']
        data = pandas.DataFrame(data_pre[data_pre.str.len() <= len_cutoff])
        if drop_degens:
            dataF = data.drop_duplicates()
        else:
            dataF = data
    if return_index:
        temp_index = dataF.index
    for i in np.arange(len(dataF)):
        # Replace 
        titleV = label + '_' + str(a+start_label)
        
        if thing:
            final_title = [titleV]
            thing = False
        else:
            final_title = final_title + [titleV]
        a = a+1

    finalDF = np.transpose(pandas.DataFrame(np.array(dataF)))
    finalDF.columns=final_title

    if scrape:
        finalAllele = np.transpose(pandas.DataFrame(np.array(allele)))
        finalAllele.columns=final_title
        return(finalDF,finalAllele)
    elif return_index:
        return(finalDF,temp_index)
    else:
        return(finalDF)

# Special Loaders for MHC Germline Analysis
####################################################################################
# A Temporary Change Here: Load in HLA A, B, C platform domains for analysis
def get_HLA():
    from Bio import SeqIO
    mhc_dir = datPath+'app_data/germline_data/'
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
    mhc_dir = datPath+'app_data/germline_data/'
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

############## Currently non-functional loading scripts ################

#######################################################################
# So as a note, this script actually loads in ALL MHC. Class I
def load_multiOrgMHC():
    records = list(SeqIO.parse(datPath+'app_data/germline_data/MHC_prot.fasta','fasta'))
    mhc_codes = pandas.read_csv(datPath+'app_data/germline_data/MHCprot_orgCode.csv')

    first = True
    for tester in records:
        ID = tester.id
        for sin_code in mhc_codes['IMGTcode']:
            find_code = tester.description.find(sin_code)
            if find_code != -1:
                metaDat = mhc_codes[mhc_codes['IMGTcode'] == sin_code]
                break
        pre_df = np.transpose(pandas.DataFrame(np.hstack([str(tester.seq),ID,metaDat.values[0]])))
        if first:
            fin_df = pre_df
            first = False
        else:
            fin_df = pandas.concat([fin_df,pre_df],axis=0)
    fin_df.columns = ['seq','ID','broadOrg', 'species', 'IMGTcode', 'commonName', 'numCode']
    pre_final_df = fin_df.drop_duplicates('seq')
    pre_final_df.index = np.arange(len(pre_final_df))

    # Specifically load in the mouse data
    records2 = list(SeqIO.parse(datPath+'app_data/germline_data/mouse_classI.fasta','fasta'))
    mouse_seqF = []
    for tester in records2:
        ID = tester.id
        seq = str(tester.seq)
        mouse_seqF = mouse_seqF + [seq]
    mouse_ids = ['H2D','H2K','H2L','H2M','H2Q','H2T']
    broadOrg = ['Mus','Mus','Mus','Mus','Mus','Mus']
    ID = ['','','','','','']
    species = ['Mouse','Mouse','Mouse','Mouse','Mouse','Mouse']
    commonName = ['Mouse','Mouse','Mouse','Mouse','Mouse','Mouse']
    numCode = [0,0,0,0,0,0]
    IMGTcode = ['Mus','Mus','Mus','Mus','Mus','Mus']

    mouse_tempdf = np.transpose(pandas.DataFrame([mouse_seqF,mouse_ids,broadOrg,species,IMGTcode,commonName,numCode]))
    mouse_tempdf.columns = pre_final_df.columns

    # Need to ALSO load in the human molecules. I THINK that allowing us to randomly subsample the 
    # MHC molecules should actually suffice to control for the heavy skewing of the dataset in number towards human
    # Loaded them in earlier, just process them into the "finalDF"
    classI = get_HLA()
    classIIa_df,classIIb_df = get_classII()

    classIIa_human_len = len(classIIa_df[0])
    classIIb_human_len = len(classIIb_df[0])
    classI_human_len = len(classI[0])

    classI_ids = classI[0]; classIIa_ids = classIIa_df['Allele']; classIIb_ids = classIIb_df['Allele']
    classI_broadOrg = ['Human']*classI_human_len; classIIa_broadOrg = ['Human']*classIIa_human_len; classIIb_broadOrg = ['Human']*classIIb_human_len 
    classI_ID = ['']*classI_human_len; classIIa_ID = ['']*classIIa_human_len; classIIb_ID = ['']*classIIb_human_len
    classI_numCode = [0]*classI_human_len; classIIa_numCode = [1]*classIIa_human_len; classIIb_numCode = [2]*classIIb_human_len;
    classI_IMGTcode = ['Hosa']*classI_human_len; classIIa_IMGTcode = ['Hosa']*classIIa_human_len; classIIb_IMGTcode = ['Hosa']*classIIb_human_len;

    classI_tempdf = np.transpose(pandas.DataFrame([classI[1].values,classI_ids,classI_broadOrg,classI_broadOrg,classI_IMGTcode,classI_broadOrg,classI_numCode]))
    classIIa_tempdf = np.transpose(pandas.DataFrame([classIIa_df[0].values,classIIa_ids,classIIa_broadOrg,classIIa_broadOrg,classIIa_IMGTcode,classIIa_broadOrg,classIIa_numCode]))
    classIIb_tempdf = np.transpose(pandas.DataFrame([classIIb_df[0].values,classIIb_ids,classIIb_broadOrg,classIIb_broadOrg,classIIb_IMGTcode,classIIb_broadOrg,classIIb_numCode]))

    classI_tempdf.columns = pre_final_df.columns
    classIIa_tempdf.columns = pre_final_df.columns
    classIIb_tempdf.columns = pre_final_df.columns

    final_df = pandas.concat([pre_final_df,mouse_tempdf,classI_tempdf,classIIa_tempdf,classIIb_tempdf])
    final_df.index = np.arange(len(final_df))
    
    return(final_df)

def load_multiOrgTCR():
    tcr_orgs = ['human','mouse','rhesus','sheeps','maMonkey','cow','goat']
    mhc_orgs = ['Human','Mouse','rhesus','Sheep','_Ma_','Cattle','Goat']

    # Define how many CDR loops are in the files you are loading
    num_loop=3
    for i in np.arange(len(tcr_orgs)):
        trav = Ig_loader(datPath+'app_data/germline_data/Ig_displays/trav_'+tcr_orgs[i]+'_cdrs.csv',label=mhc_orgs[i],loops=num_loop,drop_degens = True)
        trbv = Ig_loader(datPath+'app_data/germline_data/Ig_displays/trbv_'+tcr_orgs[i]+'_cdrs.csv',label=mhc_orgs[i],loops=num_loop,drop_degens = True)
        if tcr_orgs[i] == 'human':
            fin_trav = trav.loc[0:1]
            fin_trbv = trbv.loc[0:1]
        else:
            fin_trav = pandas.concat([fin_trav,trav],axis=1)
            fin_trbv = pandas.concat([fin_trbv,trbv],axis=1)

    return(fin_trav,fin_trbv)

def get_KIR(KIR):
    # Then go on to the loading script
    kir_dir = 'IPDKIR-Latest/fasta/'
    records = list(SeqIO.parse(kir_dir+KIR+'_prot.fasta','fasta'))
    # KIR2DL1,2,3 have D1-D2 arrangements
    # Other KIRs have D0-D2, or the 3-domain KIRs
        
    full_id = []; d1_l0_concat = []; d2_l0_concat = []
    d1_l1_concat = []; d2_l1_concat = []; 
    d1_l2_concat = []; d2_l2_concat = []; 
    connect_concat = []
    for i in np.arange(len(records)):
        rec_desc = records[i].description
        seqs_pre = str(records[i].seq)
        hla_loc = rec_desc.find('*')
        partial_id = rec_desc[hla_loc-1:]
        id_end = partial_id.find(' ')

        # So far, all of the truncated sequences end with an "N"
        # Search for it, and then shitcan these sequences
        if rec_desc.find('N') != -1:
            continue
        #############################################################################
        # This is the start of the D1 domain
        start_plat = seqs_pre.find('HRKP')
        if start_plat == -1:
            continue
        # This is the end of the D1 domain
        end_plat = seqs_pre.find('DIVI')
        if end_plat == -1:
            # Alternate possibility for this sequence
            end_plat = seqs_pre.find('DVVI')
            if end_plat == -1:
                continue
        else:
            # Otherwise need to add in the remainder of the sqeuence
            end_plat = end_plat + 4
            
        d1_seqs = seqs_pre[start_plat:end_plat]
        #############################################################################
        # This is the start of the D2 domain
        start_plat = seqs_pre.find('GLYEKP')
        if start_plat == -1:
            # Alternate possibility for this sequence
            start_plat = seqs_pre.find('GLYQKP')
            if start_plat == -1:
                continue
        # This is the end of the D2 domain
        end_plat = seqs_pre.find('VSVT')
        if end_plat == -1:
            end_plat = seqs_pre.find('VSVI')
            if end_plat == -1:
                continue
        else:
            # Otherwise need to add in the remainder of the sqeuence
            end_plat = end_plat + 4
        d2_seqs = seqs_pre[start_plat:end_plat]
        
        # Based on a few papers: Moradi and Berry, JBC 2015 // Moradi et al. Nat Comm 2021 // Fan, Long, and Wiley Nat Immuno 2001
        # We can pick out a few key regions for interacting regions of the KIRS
        # VERY specifically for the KIR D1 and D2 domains
        D1_loop0 = d1_seqs[12:22]
        D1_loop1 = d1_seqs[34:44]
        D1_loop2 = d1_seqs[63:73]
        connect = d2_seqs[0:10]
        D2_loop0 = d2_seqs[26:36]
        # Structure doesn't show any contact with this region, but appears highly variable.
        # Keep it in the analysis for now
        D2_loop1 = d2_seqs[48:58]
        D2_loop2 = d2_seqs[78:88]

        full_id = full_id +  [partial_id[:id_end]]
        d1_l0_concat = d1_l0_concat + [D1_loop0]
        d2_l0_concat = d2_l0_concat + [D2_loop0]
        d1_l1_concat = d1_l1_concat + [D1_loop1]
        d2_l1_concat = d2_l1_concat + [D2_loop1]
        d1_l2_concat = d1_l2_concat + [D1_loop2]
        d2_l2_concat = d2_l2_concat + [D2_loop2]
        connect_concat = connect_concat + [connect]
    sequence_frame = np.transpose(pandas.DataFrame(np.vstack((full_id,d1_l0_concat,d1_l1_concat,
                    d1_l2_concat,connect_concat,d2_l0_concat,d2_l1_concat,d2_l2_concat))))
    # remove duplicate entries
    nonDup_locs = sequence_frame[[1,2,3,4,5,6,7]].drop_duplicates().index
    seqFrameF = sequence_frame.loc[nonDup_locs]

    return(seqFrameF)