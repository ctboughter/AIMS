import numpy as np
import pandas
import matplotlib.pyplot as pl
import math
import matplotlib as mpl
from matplotlib import cm
from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import accuracy_score
from sklearn.utils import resample
from aims_immune import aims_loader as aimsLoad
from aims_immune import aims_classification as classy

# More special stuff for the pip version of the script....
import aims_immune
# The -11 is because the filepath includes '__init__.py'
# So we need to remove that to get our data path.
datPath = aims_immune.__file__[:-11]

# Define some initial stuff and import analysis functions:
#AA_key_old=['A','G','L','M','F','W','K','Q','E','S','P','V','I','C','Y','H','R','N','D','T']
AA_key=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
AA_key_dash = ['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','-']

# So we've got 46 orthogonal (or at least not super correlated)
# dimensions. Add them in to the matrix
# From "Hot spot prediction in protein-protein interactions by an ensemble system"
# Liu et. al. BMC Systems Biology
newnew=pandas.read_csv(datPath+'app_data/new_props')
oldold=pandas.read_csv(datPath+'app_data/old_props')

# These next 8 lines re-order the properties to guarantee that everything is in the correct order
# Would be a HUGE issue if alanine properties were given to arginine, etc.
prop_parsed = np.zeros((len(oldold),20))
for i in np.arange(len(AA_key)):
    prop_parsed[0:16,i]=oldold[AA_key[i]]

properties=np.zeros((len(newnew)+len(oldold),20))
for i in np.arange(len(AA_key)):
    properties[0:16,i]=oldold[AA_key[i]]
    properties[16:,i]=newnew[AA_key[i]]

AA_num_key_new=properties[1]
AA_num_key=np.arange(20)+1

def get_sequence_dimension(re_poly):
    num_loops,num_clones=np.shape(re_poly)
    seqlenF = []
    for i in np.arange(num_loops):
        if i == 0:
            max_len=len(re_poly[i,0])
        elif i <= num_loops:
            max_len=np.vstack((max_len,len(re_poly[i,0])))
        for j in np.arange(num_clones):
            seqlen = len(re_poly[i,j])
            if i == 0:
                if seqlen > max_len:
                    max_len = seqlen
            else:
                if seqlen > max_len[i]:
                    max_len[i] = seqlen
            if type(seqlenF) == list:
                seqlenF = seqlen
            else:
                seqlenF = np.vstack((seqlenF,seqlen))
    max_len=max_len+3 # Add 3 here so there's a more pronounced space between loops
    ## SO NOW MAX LEN SHOULD HAVE THE MAXIMUM LENGTH OF EACH CDR LOOP ##
    if num_loops == 1:
        sequence_dim = max_len
        seqlens = seqlenF.reshape(num_clones,num_loops)
    else:
        sequence_dim = int(sum(max_len))
        seqlens = seqlenF.reshape(num_clones,num_loops)
    return(max_len,sequence_dim,seqlens)

# NOTE, manuscript_arrange=False MUST be selected to run MHC analysis
# I used to re-arrange the CDR loops for a more position-accurate
# representation. In the current analysis, this isn't necessary.
def gen_tcr_matrix(pre_poly,AA_key=AA_key,key=AA_num_key_new,binary=False,
pre_mono=[],giveSize='',return_Size=False,manuscript_arrange=False,
alignment = 'center',bulge_pad = 8):
    if int(bulge_pad) != 8 and int(bulge_pad) != 6 and int(bulge_pad) != 4 and int(bulge_pad) != 2:
        # We can only accept the above options for bulge
        print('Can only accept bulge_pad of 8, 6, 4, or 2.')
        print('Default to bulge_pad = 8')
        bulge_pad = 8
    # Do this so that 
    if type(giveSize) == int:
        giveSize = [giveSize]
    if len(giveSize) == 0:
        if binary:
            # NEW ADDITION TO CLEAN THINGS UP A BIT #
            max_len1 = get_sequence_dimension(pre_poly)[0]
            max_len2 = get_sequence_dimension(pre_mono)[0]
            max_lenp=np.zeros(len(max_len1))
            for i in np.arange(len(max_len1)):
                max_lenp[i]=int(max(max_len1[i],max_len2[i]))
            if type(max_lenp) == int:
                sequence_dim = max_lenp
            else:
                sequence_dim = int(sum(max_lenp))
        else:
            max_lenp,sequence_dim,seqlens=get_sequence_dimension(pre_poly)
    else:
        max_lenp = giveSize
        if type(max_lenp) == int:
            sequence_dim = max_lenp
        else:
            sequence_dim = int(sum(max_lenp))
    final_poly=[] # initialize a variable
    # RE-ORGANIZE EVERYTHING SO IT LOOKS NICE IN MATRIX
    # But presumably, if you give a size, it's the size you want...
    if manuscript_arrange and len(giveSize) == 0:
        max_lenp = [max_lenp[0],max_lenp[1],max_lenp[2],max_lenp[5],max_lenp[4],max_lenp[3]]
    max_len = max_lenp
    
    for re_poly in [pre_poly,pre_mono]:
        # RE-ORGANIZE EVERYTHING SO IT LOOKS NICE IN MATRIX
        if manuscript_arrange:
            re_poly = [re_poly[0],re_poly[1],re_poly[2],re_poly[5],re_poly[4],re_poly[3]]

        numLoop,numClone = np.shape(re_poly)

        poly_PCA=np.zeros([numClone,sequence_dim])
        for i in range(numClone): # For all of our polyreactive sequences...
            loop=0
            for k in range(numLoop): # Scroll through all of the loops within a clone                
                leng=len(re_poly[k][i]) #should be this if re-sampling
                # the int clause below is how we center align

                if type(max_len) == int:
                    LoopMax = max_len
                    loopBuff = 0
                    if alignment.lower() == 'center':
                        count = int((max_len-leng)/2.0)
                    elif alignment.lower() == 'left':
                        count = 0
                    elif alignment.lower() == 'right':
                        count = int(max_len - leng)
                    elif alignment.lower() == 'bulge':
                        count = 0
                        b_count = 0 # b_count is bulge count
                        bulge = leng - bulge_pad
                        start_end = True # Have a trigger to start the end of the peptide
                        bulge_cen = int((max_len - bulge)/2.0)
                else:
                    # Can go ahead and replace these below... at some point
                    # For now define them for the purpose of the bulge alignment below
                    LoopMax = max_len[k]
                    loopBuff = int(sum(max_len[:k]))
                    if alignment.lower() == 'center':
                        count=int((max_len[k]-leng)/2.0)+int(sum(max_len[:k]))
                    elif alignment.lower() == 'left':
                        count = int(0) + int(sum(max_len[:k]))
                    elif alignment.lower() == 'right':
                        count = int(max_len[k] - leng) + int(sum(max_len[:k]))
                    elif alignment.lower() =='bulge':
                        # If the loop is too short, just do center alignment
                        if int(LoopMax) < int(bulge_pad+2):
                            count=int((max_len[k]-leng)/2.0)+int(sum(max_len[:k]))
                        else:
                            count = int(0) + int(sum(max_len[:k]))
                            b_count = 0
                            bulge = leng - bulge_pad
                            start_end = True
                            bulge_cen = int((max_len[k] - bulge)/2.0) + loopBuff
                for m in re_poly[k][i]: # SO IS THIS ONE for bootstrapping
                    for j in range(len(key)):
                        if m==AA_key[j]:
                            # Can't be doing this for short sequences
                            if alignment.lower() == 'bulge' and int(LoopMax) > int(bulge_pad+2):
                                if bulge < 0:
                                    if count > 3 + loopBuff:
                                        if start_end:
                                            start_end = False
                                            # So here is where we skip to the very end
                                            count = LoopMax-int(bulge_pad/2) + loopBuff
                                else:
                                    # Clearly don't need to do any of this stuff
                                    # if we DONT have a bulge, so make exception
                                    if count > int(bulge_pad/2)-1 + loopBuff and count < LoopMax-int(bulge_pad/2) +loopBuff:
                                        if b_count == 0:
                                            count = bulge_cen
                                        b_count += 1
                                    # This should then jump us over to the
                                    # end once we get past the bulk
                                    # ie, gaps should go around the bulk
                                    if b_count > bulge:
                                        if start_end:
                                            count = LoopMax-int(bulge_pad/2) + loopBuff
                                            start_end = False
                            poly_PCA[i][int(count)]=key[j]
                            count += 1
                loop=loop+1
        if binary:
            # Unfortunate naming here for the binary case
            # but it is what it is...
            if len(final_poly) == 0:
                final_poly = poly_PCA
            else:
                final_mono = poly_PCA
        else:
            break
    
    if binary and return_Size:
        return(final_poly,final_mono,max_lenp)
    elif binary:
        return(final_poly,final_mono)
    elif return_Size:
        return(poly_PCA,max_lenp)
    else:
        return(poly_PCA)

def calculate_shannon(poly_PCA):
    clones,aas = np.shape(poly_PCA)
    prob_poly_full=np.zeros((clones,aas,21))
    coverage = np.zeros(aas)
    # We technically have 21 entries, 20 AAs (1-20) and spaces 0... Take a look at all of that
    # Actually 22 now that we are allowing for '-' in MSA
    AAs=np.arange(0,22)
    #Actually for the time being, there is no "22". Whenever we see '-' we wouldn't see a space
    # so we can let them be equal for now

    for i in np.arange(clones):
        for j in np.arange((aas)):
            if poly_PCA[i,j] == 0:
                coverage[j] += 1
            for k in AAs:
                if poly_PCA[i,j]==k:
                    prob_poly_full[i,j,k]=prob_poly_full[i,j,k]+1
                    
    poly_count=np.sum(prob_poly_full,axis=0)/clones

    shannon_poly=np.zeros(len(poly_count))
    for i in np.arange(len(poly_count)):
        for j in np.arange(len(poly_count[0])):
            if poly_count[i,j]==0:
                continue
            shannon_poly[i]=shannon_poly[i]+(-poly_count[i,j]*math.log(poly_count[i,j],2))
    return(shannon_poly,poly_count,coverage/clones)

def calculate_MI(poly_PCA):
    shannon_poly,poly_count,coverage=calculate_shannon(poly_PCA)
    # Need to start removing explicit "126" sized matrices
    clones,aas = np.shape(poly_PCA)
    AAs=np.arange(0,21)
    MI_final_poly=np.zeros((aas,aas))
    poly_count_cond=np.zeros((21,aas,aas,21))
    save_count=np.zeros((21,aas))
    for location in np.arange(aas):
        # Copy same stuff as above, but now also calculate a conditional (cond) probability
        cond_count_poly=0
        #MI_final_poly=np.zeros((20,126))
        conditional_entropy_poly=np.zeros((21,aas))
        for res in AAs:
            prob_poly_cond=np.zeros((clones,aas,21))
            cond_count_poly=0
            for i in np.arange(clones): # Cycle through all clones
                if poly_PCA[i,location]==res: # Of those clones, only look at the ones that have AA 15 at position 73
                    cond_count_poly=cond_count_poly+1
                    for j in np.arange((aas)): # Then scroll through all of the positions to get conditional prob
                        for k in AAs:
                            if poly_PCA[i,j]==k:
                                prob_poly_cond[i,j,k]=prob_poly_cond[i,j,k]+1
                
        #So the matrix is built, now let's do our entropy...
        # for every position, we have a count... Can check and averaging over each position does add to 1.
            #if location == 75 and res == 5:
            #    return(prob_poly_cond,cond_count_poly)
            if cond_count_poly == 0:
                continue
            poly_count_cond[res,location,:,:]=np.sum(prob_poly_cond,axis=0)/cond_count_poly
            save_count[res,location]=cond_count_poly
    
            MI_poly=np.zeros(len(poly_count_cond[res,location,:,:])) # Note, really should change this... MI poly is a misnomer here. This is one of the conditional entropy terms...
            for i in np.arange(len(poly_count_cond[res,location,:,:])):
                for j in np.arange(len(poly_count_cond[res,location,:,:][0])):
                    if poly_count_cond[res,location,:,:][i,j]==0:
                        continue
                    MI_poly[i]=MI_poly[i]+(-poly_count_cond[res,location,:,:][i,j]*math.log(poly_count_cond[res,location,:,:][i,j],2))                
            
            condition_poly=np.zeros(len(shannon_poly))
            condition_poly[location]=shannon_poly[location]
            #if cond_count_poly==0:
            #    continue
                                
        ############################## NOTE, I THINK THIS IS WHERE MY ISSUE IS #####################
        # I need to sum up my MI_poly (which again, is actually conditional entropy) and THEN subtract from the standard entropy...
            conditional_entropy_poly[res]=MI_poly
        # Land at Shannon Entropy - global probability of a given residue * the conditional shannon entropy
        MI_final_poly[location]=shannon_poly-np.matmul(poly_count[location],conditional_entropy_poly)#-condition_poly
    return(MI_final_poly,poly_count_cond,poly_count)


def joint_prob(poly_PCA):
    #Joint Probability Distribution Time!
    clones,aas = np.shape(poly_PCA)
    # Might be a bit of a beast to calculate...
    joint_probs = np.zeros((aas,21,aas,21)) # Should be position x number amino acids x 
    AAs=np.arange(0,21)
    for cloneX in np.arange(clones): # Go through each clone in the matrix
        for posX in np.arange(aas): # GIVEN this first position in the matrix
            for resX in AAs: # cycle through all of the amino acids
                if poly_PCA[cloneX,posX] == resX: # test for existence to speed things up a bit
                    for posY in np.arange(aas): # cycle through every position 
                        for resY in AAs: #cycle through every residue at these positions
                            if poly_PCA[cloneX,posY] == resY:
                                joint_probs[posX,resX,posY,resY] = joint_probs[posX,resX,posY,resY] + 1
    joint=joint_probs/clones # divide out by the number of clones and number of Y positions
    return(joint)



def gen_dset_props(poly_PCA,props=properties,stdev=False):
    AA_num_key=np.arange(20)+1
    clones,aas = np.shape(poly_PCA)
    poly_prop=np.zeros([len(props),aas]) 
    for i in np.arange(len(props)): # For all of our properties...
        for j in np.arange(clones): # for every clone
            for k in np.arange(aas): # for every position
                for m in AA_num_key:
                    if poly_PCA[j,k]==m:
                        poly_prop[i,k]=poly_prop[i,k]+props[i][m-1]
                        
    poly_prop=poly_prop/clones

    if stdev:
        poly_prop_stdev=np.zeros([len(props),aas]) 
        for i in np.arange(len(props)): # For all of our properties...
            for j in np.arange(clones): # for every clone
                for k in np.arange(aas): # for every position
                    for m in AA_num_key:
                        if poly_PCA[j,k]==m:
                            #stdev = sqrt([sum(x-xavg)**2]/[N-1])... So below term is computing that sum
                            poly_prop_stdev[i,k]=poly_prop_stdev[i,k]+(props[i][m-1]-poly_prop[i,k])**2
                            
        poly_prop_stdev=np.sqrt(poly_prop_stdev/(clones-1))

        return(poly_prop,poly_prop_stdev)
    else:
        return(poly_prop)

def gen_clone_props(poly_PCA):
    # Difference between this and the above function is how the 
    # average is being taken... Either average over clones or average
    # over positions within the clones.
    clones,aas = np.shape(poly_PCA)
    
    # This whole code block used to be outside the function, but kept running into a very weird issue
    # where it seemed like the "props" were changing in a strange way.
    properties=np.zeros((len(newnew)+len(oldold),20))
    for i in np.arange(len(AA_key)):
        properties[0:16,i]=oldold[AA_key[i]]
        properties[16:,i]=newnew[AA_key[i]]
    props = properties[1:]

    # Re-normalize the properties for use in the matrix...
    for i in np.arange(len(props)):
        props[i] = props[i]/np.linalg.norm(props[i])

    poly_prop_pca=np.zeros([len(props),clones])
 
    for i in np.arange(len(props)): # For all of our properties...
        for j in np.arange(clones): # for every clone
            for k in np.arange(aas): # for every position
                for m in AA_num_key:
                    if poly_PCA[j,k]==m:
                        poly_prop_pca[i,j]=poly_prop_pca[i,j]+props[i,m-1]
    poly_prop_pca=poly_prop_pca/aas

    # 30 is a placeholder for now (30 positions tested)
    #single_pos_charge=np.zeros([30,clones]) 
    #single_pos_phob=np.zeros([30,clones]) 
    #for j in np.arange(clones): # for every clone
    #    a=0 # a is here as a placeholder
    #    for k in np.arange(68,98): # for every position
    #        for m in AA_num_key:
    #            if poly_PCA[j,k]==m:
                    # So I think that 0 and 1 should correspond to charge and hydrophobicity...
    #                single_pos_charge[a,j]=props[1,m-1]
    #                single_pos_phob[a,j]=AA_num_key_new[m-1]
    #        a=a+1
    return(poly_prop_pca)#,single_pos_charge,single_pos_phob) 

def get_props():
    # Identical to what we've got above, but it lets us also load it in to
    # the current variable space...
    AA_key=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

    # So we've got 46 orthogonal (or at least not super correlated)
    # dimensions. Add them in to the matrix
    # From "Hot spot prediction in protein-protein interactions by an ensemble system"
    # Liu et. al. BMC Systems Biology (these are "new props")

    newnew=pandas.read_csv(datPath+'app_data/new_props')
    oldold=pandas.read_csv(datPath+'app_data/old_props')

    props=np.zeros((len(newnew)+len(oldold),20))
    for i in np.arange(len(AA_key)):
        props[0:16,i]=oldold[AA_key[i]]
        props[16:,i]=newnew[AA_key[i]]

    AA_num_key=np.arange(20)+1
    AA_num_key_new=props[1]

    return(AA_key,AA_num_key,AA_num_key_new,props)

def prop_patterning(mono_PCA,poly_PCA,mat_size=100,props=properties[1:],ridZero=False,win_size = 3,returnBig=False):
    # Try to maximize differences across the properties by looking at patterning...

    # Re-normalize the properties for use in the matrix...
    for i in np.arange(len(props)):
        props[i] = props[i]-np.average(props[i])
        props[i] = props[i]/np.linalg.norm(props[i])

    # Since we'll be averaging shit, let's get rid of all the zeros...
    # However, that's going to result in bleed-over across the loops... Is this good or bad?
    # This is also going to have a strange effect based upon loop length...
    # WHATEVER, Try both
    
    if ridZero:
        poly_pca_NEW=np.transpose(poly_PCA)[~np.all(np.transpose(poly_PCA) == 0,axis=1)]
        poly_pca_NEW=np.transpose(poly_pca_NEW)
        mono_pca_NEW=np.transpose(mono_PCA)[~np.all(np.transpose(mono_PCA) == 0,axis=1)]
        mono_pca_NEW=np.transpose(mono_pca_NEW)
    else:
        poly_pca_NEW = poly_PCA
        mono_pca_NEW = mono_PCA

    poly_dim1,poly_dim2=np.shape(poly_pca_NEW)
    mono_dim1,mono_dim2=np.shape(mono_pca_NEW)

    # So this is where we should be able to do the averaging
    poly_prop_masks=np.zeros([len(props),poly_dim1,int(poly_dim2/win_size)])
    mono_prop_masks=np.zeros([len(props),mono_dim1,int(mono_dim2/win_size)])

    for i in np.arange(len(props)): # For all of our properties...
        for j in np.arange(poly_dim1): # for every clone
            for k in np.arange(int(poly_dim2/win_size)): # for every position
                for win_slide in np.arange(win_size):
                    # Make sure we don't go too far...
                    if k*win_slide+win_slide > int(poly_dim2/win_size):
                        continue
                    for m in AA_num_key:
                        if poly_pca_NEW[j,k*win_size+win_slide]==m:
                            poly_prop_masks[i,j,k]=poly_prop_masks[i,j,k]+props[i,m-1]

    for i in np.arange(len(props)): # For all of our properties...
        for j in np.arange(mono_dim1): # for every clone
            for k in np.arange(int(mono_dim2/win_size)): # for every position
                for win_slide in np.arange(win_size):
                    # Make sure we don't go too far...
                    if k*win_slide+win_slide > int(mono_dim2/win_size):
                        continue
                    for m in AA_num_key:
                        if mono_pca_NEW[j,k*win_size+win_slide]==m:
                            mono_prop_masks[i,j,k]=mono_prop_masks[i,j,k]+props[i,m-1]

    # And now we actually have to do some sort of intelligent analysis on these masks...
    mono_prop_line = np.average(mono_prop_masks/win_size,axis=1)
    poly_prop_line = np.average(poly_prop_masks/win_size,axis=1)

    line_diff = poly_prop_line - mono_prop_line

    # Pull out the biggest differences in these properties...
    max_diffs = np.zeros((mat_size,3))
    for i in np.arange(len(line_diff)):
        for j in np.arange(len(line_diff[0])):
            if line_diff[i,j] == 0:
                continue
            elif abs(line_diff[i,j]) > min(abs(max_diffs[:,0])):
                # Find the absolute value minimum.
                new_pos=np.where(abs(max_diffs[:,0]) == min(abs(max_diffs[:,0])))
                # Put_it_here should get the location of the min...
                # It should always be in max_diffs[X,0]
                put_it_here=int(new_pos[0][0])
                max_diffs[put_it_here,0] = line_diff[i,j]
                max_diffs[put_it_here,1] = i
                # Need the *win_size to revert back to physical location on the CDR loops...
                # Actually NOT true if we get rid of the zeros... little trickier there, come back to it
                max_diffs[put_it_here,2] = j*win_size

    # Now bring on back these properties for classification...
    new_mat_poly=np.zeros((mat_size,len(poly_pca_NEW))) 
    new_mat_mono=np.zeros((mat_size,len(mono_pca_NEW))) 

    for j in np.arange(len(poly_pca_NEW)): # for every clone
        for k in np.arange(len(max_diffs)): 
            for win_slide in np.arange(win_size):
                # This really shouldn't happen, but just in case...
                if max_diffs[k,2]+win_slide > len(poly_pca_NEW[0]):
                    continue
                for m in AA_num_key:
                    if poly_pca_NEW[j,int(max_diffs[k,2]+win_slide)]==m:
                        # So I think that 0 and 1 should correspond to charge and hydrophobicity...
                        new_mat_poly[k,j]=new_mat_poly[k,j] + props[int(max_diffs[k,1]),m-1]

    for j in np.arange(len(mono_pca_NEW)): # for every clone
        for k in np.arange(len(max_diffs)): 
            for win_slide in np.arange(win_size):
                # This really shouldn't happen, but just in case...
                if max_diffs[k,2]+win_slide > len(mono_pca_NEW[0]):
                    continue
                for m in AA_num_key:
                    if mono_pca_NEW[j,int(max_diffs[k,2]+win_slide)]==m:
                        # So I think that 0 and 1 should correspond to charge and hydrophobicity...
                        new_mat_mono[k,j]=new_mat_mono[k,j] + props[int(max_diffs[k,1]),m-1]
    if returnBig:
        return(new_mat_mono/win_size,new_mat_poly/win_size,max_diffs,poly_prop_masks,mono_prop_masks)
    else:
        return(new_mat_mono/win_size,new_mat_poly/win_size,max_diffs)

# Alright so this below section is identical to the above code... For now. Need to rationalize how to pair/score

def prop_pairing(ALL_mono,ALL_poly,mat_size=100,props=properties[1:],win_size = 3):
    # Try to maximize differences across the properties by looking at patterning...
    max_len1=get_sequence_dimension(ALL_poly)[0]
    max_len2=get_sequence_dimension(ALL_mono)[0]
    max_lenp=np.zeros(6)
    for i in np.arange(6):
        max_lenp[i]=max(max_len1[i],max_len2[i])
    max_size=int(max(max_lenp))

    # TEMPORARY
    #max_size=21

    poly_PCA=gen_tcr_matrixOLD(ALL_poly,max_size,key=AA_num_key_new)
    mono_PCA=gen_tcr_matrixOLD(ALL_mono,max_size,key=AA_num_key_new)

    # For this one (and maybe the one above) we should reshape the matrix
    poly_dim1,poly_dim2=np.shape(poly_PCA)
    mono_dim1,mono_dim2=np.shape(mono_PCA)
    # This reshapes the matrices in to #clones by CDR loop by arbitrary position
    poly_pca_re = poly_PCA.reshape(int(poly_dim1),6,int(poly_dim2/6))
    mono_pca_re = mono_PCA.reshape(int(mono_dim1),6,int(mono_dim2/6))

    # FOR THIS SPECIFIC APPLICATION, I THINK ONLY PROPERTIES 1,2 apply (charge and hydrophobicity)
    props=props[1:3,:]
    # So this is where we should be able to do the averaging... So now win_size can't be over 21
    poly_prop_masks=np.zeros([len(props),int(poly_dim1),6,int(poly_dim2/6/win_size)])
    mono_prop_masks=np.zeros([len(props),int(mono_dim1),6,int(mono_dim2/6/win_size)])

    # Make the property mask for the polyreactive matrix
    for i in np.arange(len(props)): # For all of our properties...
        for j in np.arange(poly_dim1): # for every clone
            for k in np.arange(int(poly_dim2/6/win_size)): # for every position
                for loop in np.arange(6): # for all 6 loops
                    for win_slide in np.arange(win_size):
                        # Make sure we don't go too far...
                        if k*win_slide+win_slide > int(poly_dim2/6/win_size):
                            continue
                        for m in AA_num_key:
                            if poly_pca_re[j,loop,k*win_size+win_slide]==m:
                                poly_prop_masks[i,j,loop,k]=poly_prop_masks[i,j,loop,k]+props[i,m-1]
    # Make the property mask for the monoreactive matrix
    for i in np.arange(len(props)): # For all of our properties...
        for j in np.arange(mono_dim1): # for every clone
            for k in np.arange(int(mono_dim2/6/win_size)): # for every position
                for loop in np.arange(6): # for all 6 loops
                    for win_slide in np.arange(win_size):
                        # Make sure we don't go too far...
                        if k*win_slide+win_slide > int(mono_dim2/6/win_size):
                            continue
                        for m in AA_num_key:
                            if mono_pca_re[j,loop,k*win_size+win_slide]==m:
                                mono_prop_masks[i,j,loop,k]=mono_prop_masks[i,j,loop,k]+props[i,m-1]

    # Alright now I've actually got to go through and determine if there are positive interactions between loops
    # FOR NOW, make the assumption that only AAs at the same position "see" each other... This is less egregious when averaging over windows
    Mmask_dim1,Mmask_dim2,Mmask_dim3,Mmask_dim4 = np.shape(mono_prop_masks)
    Pmask_dim1,Pmask_dim2,Pmask_dim3,Pmask_dim4 = np.shape(poly_prop_masks)

    poly_pair_score=np.zeros([6,6,Pmask_dim1,Pmask_dim2,Pmask_dim4])
    mono_pair_score=np.zeros([6,6,Mmask_dim1,Mmask_dim2,Mmask_dim4])

    for i in np.arange(Mmask_dim1): # For all of our properties
        for j in np.arange(Mmask_dim2): # for every clone
            for k in np.arange(Mmask_dim4): # for every position
                for loop1 in np.arange(6): # This and the bottom loop should be enough to find out loop-loop scoring
                    for loop2 in np.arange(loop1+1): # The plus 1 here also includes self-loop
                        if i == 0: # i=0 is charge as a property, so like-like is bad...
                            if (mono_prop_masks[i,j,loop1,k] > 0.0) & (mono_prop_masks[i,j,loop2,k] > 0.0):
                                mono_pair_score[loop1,loop2,i,j,k] = -1
                                mono_pair_score[loop2,loop1,i,j,k] = -1
                            elif (mono_prop_masks[i,j,loop1,k] < 0.0) & (mono_prop_masks[i,j,loop2,k] < 0.0):
                                mono_pair_score[loop1,loop2,i,j,k] = -1
                                mono_pair_score[loop2,loop1,i,j,k] = -1
                            elif (mono_prop_masks[i,j,loop1,k] < 0.0) & (mono_prop_masks[i,j,loop2,k] > 0.0):
                                mono_pair_score[loop1,loop2,i,j,k] = 1
                                mono_pair_score[loop2,loop1,i,j,k] = 1
                            elif (mono_prop_masks[i,j,loop1,k] > 0.0) & (mono_prop_masks[i,j,loop2,k] < 0.0):
                                mono_pair_score[loop1,loop2,i,j,k] = 1
                                mono_pair_score[loop2,loop1,i,j,k] = 1
                        elif i == 1: # i=1 is hydrophobicity as a property, so like-like is good!
                            if (mono_prop_masks[i,j,loop1,k] > 0.0) & (mono_prop_masks[i,j,loop2,k] > 0.0):
                                mono_pair_score[loop1,loop2,i,j,k] = 1
                                mono_pair_score[loop2,loop1,i,j,k] = 1
                            elif (mono_prop_masks[i,j,loop1,k] < 0.0) & (mono_prop_masks[i,j,loop2,k] < 0.0):
                                mono_pair_score[loop1,loop2,i,j,k] = 1
                                mono_pair_score[loop2,loop1,i,j,k] = 1
                            elif (mono_prop_masks[i,j,loop1,k] < 0.0) & (mono_prop_masks[i,j,loop2,k] > 0.0):
                                mono_pair_score[loop1,loop2,i,j,k] = -1
                                mono_pair_score[loop2,loop1,i,j,k] = -1
                            elif (mono_prop_masks[i,j,loop1,k] > 0.0) & (mono_prop_masks[i,j,loop2,k] < 0.0):
                                mono_pair_score[loop1,loop2,i,j,k] = -1
                                mono_pair_score[loop2,loop1,i,j,k] = -1

    for i in np.arange(Pmask_dim1): # For all of our properties
        for j in np.arange(Pmask_dim2): # for every clone
            for k in np.arange(Pmask_dim4): # for every position
                for loop1 in np.arange(6): # This and the bottom loop should be enough to find out loop-loop scoring
                    for loop2 in np.arange(loop1+1): # The plus 1 here also includes self-loop
                        if i == 0: # i=0 is charge as a property, so like-like is bad...
                            if (poly_prop_masks[i,j,loop1,k] > 0.0) & (poly_prop_masks[i,j,loop2,k] > 0.0):
                                poly_pair_score[loop1,loop2,i,j,k] = -1
                                poly_pair_score[loop2,loop1,i,j,k] = -1
                            elif (poly_prop_masks[i,j,loop1,k] < 0.0) & (poly_prop_masks[i,j,loop2,k] < 0.0):
                                poly_pair_score[loop1,loop2,i,j,k] = -1
                                poly_pair_score[loop2,loop1,i,j,k] = -1
                            elif (poly_prop_masks[i,j,loop1,k] < 0.0) & (poly_prop_masks[i,j,loop2,k] > 0.0):
                                poly_pair_score[loop1,loop2,i,j,k] = 1
                                poly_pair_score[loop2,loop1,i,j,k] = 1
                            elif (poly_prop_masks[i,j,loop1,k] > 0.0) & (poly_prop_masks[i,j,loop2,k] < 0.0):
                                poly_pair_score[loop1,loop2,i,j,k] = 1
                                poly_pair_score[loop2,loop1,i,j,k] = 1
                        elif i == 1: # i=1 is hydrophobicity as a property, so like-like is good!
                            if (poly_prop_masks[i,j,loop1,k] > 0.0) & (poly_prop_masks[i,j,loop2,k] > 0.0):
                                poly_pair_score[loop1,loop2,i,j,k] = 1
                                poly_pair_score[loop2,loop1,i,j,k] = 1
                            elif (poly_prop_masks[i,j,loop1,k] < 0.0) & (poly_prop_masks[i,j,loop2,k] < 0.0):
                                poly_pair_score[loop1,loop2,i,j,k] = 1
                                poly_pair_score[loop2,loop1,i,j,k] = 1
                            elif (poly_prop_masks[i,j,loop1,k] < 0.0) & (poly_prop_masks[i,j,loop2,k] > 0.0):
                                poly_pair_score[loop1,loop2,i,j,k] = -1
                                poly_pair_score[loop2,loop1,i,j,k] = -1
                            elif (poly_prop_masks[i,j,loop1,k] > 0.0) & (poly_prop_masks[i,j,loop2,k] < 0.0):
                                poly_pair_score[loop1,loop2,i,j,k] = -1
                                poly_pair_score[loop2,loop1,i,j,k] = -1

    # And now we actually have to do some sort of intelligent analysis on these masks...
    #mono_prop_line = np.average(mono_prop_masks/win_size,axis=1)
    #poly_prop_line = np.average(poly_prop_masks/win_size,axis=1)

    return(poly_pair_score,mono_pair_score,poly_prop_masks,mono_prop_masks)

# OK, SO RATHER THAN change up this entire script, lemme just bring back in the OLD
# gen_tcr_matrix script and use the old prop pairing...
# What do I mean by this... It means one of the scripts requires each loop to be
# the same number of entries. New way of making tcr_matrix doesn't do that.
# So here, we just make 6 loop entries that are of length = max_len
def gen_tcr_matrixOLD(re_poly,max_len,AA_key=AA_key,key=AA_num_key_new):
    poly_PCA=np.zeros([len(re_poly[0]),6*max_len])
    for i in range(len(re_poly[:][0])): # For all of our polyreactive sequences...
        loop=0
        for k in [0,1,2,5,4,3]: # Scroll through all of the loops within a clone
            #if poly[k][i][0]=='':
            #    continue
            #leng=len(re_poly[k][i][0]) # THIS LINE IS AN ISSUE WHEN RESAMPLING
            leng=len(re_poly[k][i]) #should be this if re-sampling
            count=0+max_len*loop
            if leng<max_len:
                count=int((max_len-leng)/2.0)+max_len*loop
            for m in re_poly[k][i]: # SO IS THIS ONE for bootstrapping
            #for m in re_poly[k][i][0]: # Need this version of the script when not doing bootstrap
                # The below commented out section is to go by individual AAs...
                for j in range(len(key)):
                    if m==AA_key[j]:
                        poly_PCA[i][count]=key[j]
                        count=count+1
            loop=loop+1
    return(poly_PCA)

###############################################################################################################

# OK SO NOW YOU GOTTA BE ABLE TO DO ALL OF THIS ANALYSIS ON A SINGLE CHAIN
# Everything else should be able to work downstream of this now... probably.
def gen_1Chain_matrix(pre_poly,AA_key=AA_key,key=AA_num_key_new,binary=False,pre_mono=[],giveSize='',return_Size=False):
    # Do this so that 
    if type(giveSize) == int:
        giveSize = [giveSize]
    if len(giveSize) == 0:
        if binary:
            # NEW ADDITION TO CLEAN THINGS UP A BIT #
            max_len1=get_sequence_dimension(pre_poly)[0]
            max_len2=get_sequence_dimension(pre_mono)[0]
            max_lenp=np.zeros(3)
            for i in np.arange(3):
                max_lenp[i]=max(max_len1[i],max_len2[i])

            sequence_dim = int(sum(max_lenp))
        else:
            max_lenp,sequence_dim,seqlens=get_sequence_dimension(pre_poly)
    else:
        max_lenp = giveSize
        sequence_dim = int(sum(max_lenp))

    max_len = max_lenp
    final_poly=[] # initialize a variable
    for re_poly in [pre_poly,pre_mono]:

        poly_PCA=np.zeros([len(re_poly[0]),sequence_dim])
        for i in range(len(re_poly[:][0])): # For all of our polyreactive sequences...
            loop=0
            for k in range(3): # Scroll through all of the loops within a clone
                leng=len(re_poly[k][i]) #should be this if re-sampling
                # the int clause below is how we center align
                count=int((max_len[k]-leng)/2.0)+int(sum(max_len[:k]))
                for m in re_poly[k][i]: # SO IS THIS ONE for bootstrapping
                    for j in range(len(key)):
                        if m==AA_key[j]:
                            poly_PCA[i][count]=key[j]
                            count=count+1
                loop=loop+1
        if binary:
            # Unfortunate naming here for the binary case
            # but it is what it is...
            if len(final_poly) == 0:
                final_poly = poly_PCA
            else:
                final_mono = poly_PCA
        else:
            break
    
    if binary and return_Size:
        return(final_poly,final_mono,max_lenp)
    elif binary:
        return(final_poly,final_mono)
    elif return_Size:
        return(poly_PCA,max_lenp)
    else:
        return(poly_PCA)

#### K.I.S.S. just make a new script to get the big matrix:
def getBig(mono_PCA,AA_key=AA_key, norm = 'msuv',prop_parse=False):
    # Try to maximize differences across the properties by looking at patterning...
    # Redifine "properties" because I was getting some weird errors...
        # Prop_parse removes the "hotspot" variables for better physical
    # interpretability (a key moving forward)
    if prop_parse:
        properties=np.zeros((len(oldold),20))
        for i in np.arange(len(AA_key)):
            properties[0:16,i]=oldold[AA_key[i]]
    else:
        properties=np.zeros((len(newnew)+len(oldold),20))
        for i in np.arange(len(AA_key)):
            properties[0:16,i]=oldold[AA_key[i]]
            properties[16:,i]=newnew[AA_key[i]]

    # Reminder, we skip two here because the old properties have homemade
    # amino acid keys that may not be physically meaningful.
    props = properties[2:]

    # Re-normalize the properties for use in the matrix...
    # msuv = Mean-subtracted unit vector. Not sure why I did it this way
    if norm=='msuv':
        for i in np.arange(len(props)):
            props[i] = props[i]-np.average(props[i])
            props[i] = props[i]/np.linalg.norm(props[i])
    elif norm=='zscore':
        for i in np.arange(len(props)):
            props[i] = props[i] - np.average(props[i])
            props[i] = props[i]/np.std(props[i])
    elif norm =='0to1':
        for i in np.arange(len(props)):
            props[i] = props[i] - np.min(props[i])
            props[i] = props[i]/(np.max(props[i])-np.min(props[i]))

    mono_pca_NEW = mono_PCA

    mono_dim1,mono_dim2=np.shape(mono_pca_NEW)

    # So this is where we should be able to do the averaging
    mono_prop_masks=np.zeros([len(props),mono_dim1,int(mono_dim2)])

    for i in np.arange(len(props)): # For all of our properties...
        for j in np.arange(mono_dim1): # for every clone
            for k in np.arange(int(mono_dim2)): # for every position
                # Hopefully this should speed things up a *tiny* bit
                if mono_pca_NEW[j,k]==0:
                    continue
                for m in AA_num_key:
                    if mono_pca_NEW[j,k]==m:
                        mono_prop_masks[i,j,k]=mono_prop_masks[i,j,k]+props[i,m-1]

    return(mono_prop_masks)

# So basically this and the above code are the "prop_patterning"
# code broken up in two so it makes more sense in a classifier.
def parse_props(X_train,y_train,mat_size=100):
    # ACTUALLY HAD TO CHANGE THIS FROM THE ORI VERSION, BECAUSE WE NOW HAVE A DIFFERENT
    # SHAPE TO THE MATRIX. IT WAS PROPxCLONExmatsize
    # IS NOW CLONExBIGMAT

    # ALSO, NEED TO TAKE AS INPUT X_train, y_train
    a = True
    b = True
    for i in np.arange(np.shape(y_train)[0]):
        if y_train[i] == 1:
            if a:
                mono_prop_masks = X_train[:,i]
                a = False
            else:
                mono_prop_masks = np.vstack((mono_prop_masks,X_train[:,i]))
        elif y_train[i] == 2:
            if b:
                poly_prop_masks = X_train[:,i]
                b = False
            else:
                poly_prop_masks = np.vstack((poly_prop_masks,X_train[:,i]))

    mono_prop_line = np.average(mono_prop_masks,axis=0)
    poly_prop_line = np.average(poly_prop_masks,axis=0)

    line_diff = poly_prop_line - mono_prop_line

    # need this to be size mat_size,2 because the 2 entries for each max discriminating
    # give the magnitude of diff and the location of diff
    max_diffs = np.zeros((mat_size,3))
    # Pull out the biggest differences in these properties...
    for i in np.arange(len(line_diff)):
        if line_diff[i] == 0:
            continue
        elif abs(line_diff[i]) > min(abs(max_diffs[:,0])):
            # Find the absolute value minimum.
            new_pos=np.where(abs(max_diffs[:,0]) == min(abs(max_diffs[:,0])))
            # Put_it_here should get the location of the min...
            # It should always be in max_diffs[X,0]
            put_it_here=int(new_pos[0][0])
            max_diffs[put_it_here,0] = line_diff[i]
            max_diffs[put_it_here,1] = i
            # WE NO LONGER CARE ABOUT REVERTING TO PHYSICAL LOCATION (for this version)
            # Need the *win_size to revert back to physical location on the CDR loops...
            #max_diffs[put_it_here,2] = i*win_size
    
    #Think I can just delete all of the other shit and return max_diffs
    # Then I can just pull out those locations.
    return(max_diffs)

###################################################
# Peptide stuff:
def gen_peptide_matrix(pre_pep1,AA_key=AA_key,key=AA_num_key_new,binary=False,pre_pep2=[]):
    # How many sequences do we have?
    numClone = len(pre_pep1)
    final_pep1=[] # initialize a variable
    # Allow for the possibility that you are doing a binary comparison
    for re_pep in [pre_pep1, pre_pep2]:
        numClone = len(re_pep[0])
        ### FOR NOW, HARD CODE A PEPTIDE SEQUENCE LEN MAX OF 18
        ### MIGHT NEED TO GET CREATIVE WITH THIS, MIGHT NEED TO
        ### ADAPT IT AS WE GO, JUST IN CASE.
        sequence_dim = 14
        pep_PCA=np.zeros([numClone,sequence_dim])

        for i in range(numClone): # For all of our polyreactive sequences...

            # Check the work of Guillame et al. [PNAS 2018]
            # for some ideas on how to encode this matrix

            # This loop is where we put these rules for encoding
            # For now, let's make simple assumptions
            # (pos 2 and last pos are anchors)
            # BUT code up a "bulge" region that is centrally aligned
            tot_len = len(re_pep[0][i])

            # NOTE BULGE SHOULD ALWAYS BE > 0!!! 
            # IF PEPTIDES OF LENGTH 8 EXIST, CHANGE TO -7
            bulge = tot_len - 8
            bulge_cen = int((sequence_dim - bulge)/2.0)

            count = 0
            b_count = 0
            start_end = True
            for m in re_pep[0][i]:
                if count > 3 and count < sequence_dim-4:
                    if b_count == 0:
                        count = bulge_cen
                        b_count = b_count + 1
                    else:
                        b_count = b_count + 1
                # This should then jump us over to the
                # end once we get past the bulk
                # ie, gaps should go around the bulk
                if b_count > bulge:
                    if start_end:
                        count = sequence_dim-4
                        start_end = False

                for j in range(len(key)):
                    if m==AA_key[j]:
                        pep_PCA[i][count]=key[j]
                        count=count+1
        if binary:
            # Unfortunate naming here for the binary case
            # but it is what it is...
            if len(final_pep1) == 0:
                final_pep1 = pep_PCA
            else:
                final_pep2 = pep_PCA
        else:
            break
    
    if binary:
        return(final_pep1,final_pep2)
    else:
        return(pep_PCA)

###### NEW PARALLEL PROCESSING SECTION! ##############
# Not for running parallel processing (that's a pain)
# Instead, do it for input/output for multiprocessing.
def gen_splits(splitMat, splitSize = 100):
    s1 = np.arange(0,np.shape(splitMat)[1],splitSize)
    s2 = np.arange(splitSize,np.shape(splitMat)[1],splitSize)
    if len(s1) != len(s2):
        if len(s1) > len(s2):
            s2 = np.hstack((s2,np.shape(splitMat)[1]))
        else:
            s1 = np.hstack((s1,np.shape(splitMat)[1]))
    final = np.transpose(np.vstack((s1,s2)))
    return(final)

def compile_MP(bigass_pre, pg1, pg2, final_size = 10, prop_parse=False, cat = True):
    # Alright this concatenate function should put the final product in the correct form
    if cat:
        total_mat = np.concatenate(bigass_pre, axis = 0)
    else:
        total_mat = bigass_pre
    prop_list_old = ['Phobic1','Charge','Phobic2','Bulk','Flex','Kid1','Kid2','Kid3','Kid4','Kid5','Kid6','Kid7','Kid8','Kid9','Kid10']
    prop_list_new = ['Hot'+str(b+1) for b in range(46)]

    if prop_parse:
        prop_names = prop_list_old
    else:
        prop_names = prop_list_old + prop_list_new
    num_locs = int(np.shape(bigass_pre)[1]/len(prop_names))
    Bigass_names = []
    for i in prop_names:
        for j in np.arange(num_locs):
            Bigass_names = Bigass_names + [ i + '-' + str(j) ]

    # FROM HERE DOWN, JUST DO_LINEAR SPLIT SCRIPT
    full_big = pandas.DataFrame(total_mat,columns = Bigass_names)
    drop_zeros = [column for column in full_big.columns if all(full_big[column] == 0 )]
    y = full_big.drop(full_big[drop_zeros], axis=1)
    #z = y.corr().abs()
    z_pre = np.abs(np.corrcoef(np.transpose(y)))
    z = pandas.DataFrame(z_pre,columns=y.columns,index=y.columns)
    # Select upper triangle of correlation matrix
    upper = z.where(np.triu(np.ones(z.shape), k=1).astype(bool))

    to_drop = [column for column in upper.columns if ( any(upper[column] > 0.75) ) ]

    parsed_mat = y.drop(y[to_drop], axis=1)

    Y_train = np.hstack((np.ones(np.shape(pg1)[1]),2*np.ones(np.shape(pg2)[1])))
    ID_big = pandas.concat([full_big,pandas.DataFrame(Y_train,columns=['ID'])],axis=1)

    #####################################################
    dframe_IDed = pandas.concat([parsed_mat,pandas.DataFrame(Y_train,columns=['ID'])],axis=1)
    mono_prop_masks = dframe_IDed[dframe_IDed['ID'] == 1.0]
    poly_prop_masks = dframe_IDed[dframe_IDed['ID'] == 2.0]
    mono_prop_line = np.average(mono_prop_masks,axis = 0)
    poly_prop_line = np.average(poly_prop_masks,axis = 0)
    # remove that one extra ID column here
    line_diff = (poly_prop_line - mono_prop_line)[:-1]
    # Take the absolute value of the differences
    parsed_vect_len = np.shape(parsed_mat)[1]
    diff_dframe = pandas.DataFrame(np.abs(line_diff).reshape(1,parsed_vect_len),columns = parsed_mat.columns)
    sort_diff = diff_dframe.sort_values(0,axis = 1)
    top_diffs = sort_diff.values[:,-final_size:]
    top_names = sort_diff.columns[-final_size:]
    ######################################################

    train_mat = np.array(parsed_mat[top_names])

    clf_all = LinearDiscriminantAnalysis(n_components=1,solver='svd')    
    mda_all=clf_all.fit_transform(train_mat,Y_train)

    # NOTE, THIS IS LIKE A "BEST CASE-SCENARIO" Accuracy
    # BUT, this does tell you how to best discriminate two classes.
    p_all = clf_all.predict(train_mat)
    acc_all = accuracy_score(Y_train.flatten(),p_all)
    # Give me the coefficients
    weights=clf_all.coef_
    return(ID_big, weights, acc_all, mda_all, parsed_mat, top_names)

def split_reshape(ID_big, matShape, total_props = 61):
    seq1_bigProps = np.array(ID_big[ID_big['ID'] == 1.0])[:,:-1]
    seq2_bigProps = np.array(ID_big[ID_big['ID'] == 2.0])[:,:-1]
    cloneNum1 = len(seq1_bigProps)
    cloneNum2 = len(seq2_bigProps)
    seq1_bigReshape = seq1_bigProps.reshape(cloneNum1,total_props,matShape)
    seq2_bigReshape = seq2_bigProps.reshape(cloneNum2,total_props,matShape)
    return(seq1_bigReshape,seq2_bigReshape)

###### Dpr-DIP Matrix #######
def gen_MSA_matrix(pre_poly,AA_key_dash=AA_key_dash,key=AA_num_key_new,binary=False,
pre_mono=[],giveSize='',return_Size=False):
    # Do this so that 
    if type(giveSize) == int:
        giveSize = [giveSize]
    if len(giveSize) == 0:
        if binary:
            # NEW ADDITION TO CLEAN THINGS UP A BIT #
            max_len1 = get_sequence_dimension(pre_poly)[0]
            max_len2 = get_sequence_dimension(pre_mono)[0]
            max_lenp=np.zeros(len(max_len1))
            for i in np.arange(len(max_len1)):
                max_lenp[i]=int(max(max_len1[i],max_len2[i]))
            if type(max_lenp) == int:
                sequence_dim = max_lenp
            else:
                sequence_dim = int(sum(max_lenp))
        else:
            max_lenp,sequence_dim,seqlens=get_sequence_dimension(pre_poly)
    else:
        max_lenp = giveSize
        if type(max_lenp) == int:
            sequence_dim = max_lenp
        else:
            sequence_dim = int(sum(max_lenp))
    final_poly=[] # initialize a variable

    max_len = max_lenp
    
    for re_poly in [pre_poly,pre_mono]:

        numLoop,numClone = np.shape(re_poly)

        poly_PCA=np.zeros([numClone,sequence_dim])
        for i in range(numClone): # For all of our polyreactive sequences...
            loop=0
            for k in range(numLoop): # Scroll through all of the loops within a clone                
                leng=len(re_poly[k][i]) #should be this if re-sampling
                # the int clause below is how we center align

                if type(max_len) == int:
                    count = int((max_len-leng)/2.0)
                else:
                    # Can go ahead and replace these below... at some point
                    # For now define them for the purpose of the bulge alignment below
                    count = int((max_len[k]-leng)/2.0) + int(sum(max_len[:k]))
                for m in re_poly[k][i]: # SO IS THIS ONE for bootstrapping
                    for j in range(len(key)):
                        if m==AA_key_dash[j]:
                            poly_PCA[i][int(count)]=key[j]
                            count += 1
                loop=loop+1
        if binary:
            # Unfortunate naming here for the binary case
            # but it is what it is...
            if len(final_poly) == 0:
                final_poly = poly_PCA
            else:
                final_mono = poly_PCA
        else:
            break
    
    if binary and return_Size:
        return(final_poly,final_mono,max_lenp)
    elif binary:
        return(final_poly,final_mono)
    elif return_Size:
        return(poly_PCA,max_lenp)
    else:
        return(poly_PCA)

def get_interaction_score(seq,MSA=True,scorMat='v1'):
    if scorMat == 'v1':
        scoring_mat = pandas.read_csv(datPath+'app_data/AA_interactionV1.csv',sep=',')
    elif scorMat == 'v2':
        scoring_mat = pandas.read_csv(datPath+'app_data/AA_interactionV2.csv',sep=',')
    elif scorMat == 'v3':
        scoring_mat = pandas.read_csv(datPath+'app_data/AA_interactionV3.csv',sep=',')
    scoring_mat.index = scoring_mat['Residue'].values
    scoring_pre = scoring_mat.drop(labels='Residue',axis=1)
    dash1 = np.transpose(pandas.DataFrame(np.zeros(20)))
    dash1.columns = scoring_mat['Residue'].values; dash1.index = ["-"]
    dash2 = pandas.DataFrame(np.zeros(20))
    dash2.index = scoring_mat['Residue'].values; dash2.columns= ["-"]
    scoreF = pandas.concat([scoring_pre,dash1],axis=0)
    scoring_final = pandas.concat([scoreF,dash2],axis=1)
    scoring_final['-']['-'] = 0.0

    # So seq input HAS to be 2 columns x N rows for each sample.
    # 2 being the interacting proteins, so let's get back out the length
    # Original version of the script assumed MSA input, but would like to expand that
    if MSA:
        len1 = len(seq.values[0,0])
        len2 = len(seq.values[1,0])
    else:
        # Find the maxlens
        len1 = 0
        for i in np.arange(len(seq.values[0])):
            if len(seq.values[0,i]) > len1:
                len1 = len(seq.values[0,i])
        len2 = 0
        for j in np.arange(len(seq.values[1])):
            if len(seq.values[1,j]) > len2:
                len2 = len(seq.values[j,i])

    for i in np.arange(np.shape(seq)[1]):
        cognate_scores = np.zeros([len1,len2])
        dpr_temp = seq.values[0,i]
        dip_temp = seq.values[1,i]
        # All of this just to center align sequences
        # of different lengths
        if len(dpr_temp) < len1:
            difflen = len1 - len(dpr_temp)
            split = difflen/2
            if split < 1:
                dpr_temp = dpr_temp + '-'
            else:
                pad = int(np.floor(split))
                dpr_temp = pad*'-' + dpr_temp + pad*'-'
                dpr_temp = dpr_temp + '-'*int(np.ceil(split)-pad)
        if len(dip_temp) < len2:
            difflen = len2 - len(dip_temp)
            split = difflen/2
            if split < 1:
                dip_temp = dip_temp + '-'
            else:
                pad = int(np.floor(split))
                dip_temp = pad*'-' + dip_temp + pad*'-'
                dip_temp = dip_temp + '-'*int(np.ceil(split)-pad)

        for a in np.arange(len(dpr_temp)):
            for b in np.arange(len(dip_temp)):
                res1 = dpr_temp[a]
                res2 = dip_temp[b]
                cognate_scores[a,b] = scoring_final[res1][res2]
        if i == 0:
            cognate_final = cognate_scores
        else:
            cognate_final = np.vstack((cognate_final,cognate_scores))
        
    cognate_scores = cognate_final.reshape(i+1,len1,len2)
    cognate_average = np.average(cognate_scores,axis=0)
    cognate_std = np.std(cognate_scores,axis=0)
    return(cognate_scores,cognate_average,cognate_std)

def create_msa_pairs(file1, file2, pair_list, names, label = 'Pair'):
    csv1 = pandas.read_csv(file1,sep=','); csv2 = pandas.read_csv(file2,sep=',')
    csv1_seq = np.transpose(pandas.DataFrame(csv1['Sequence']))
    csv2_seq = np.transpose(pandas.DataFrame(csv2['Sequence']))
    
    csv1_seq.columns = csv1[names[0]]; csv2_seq.columns = csv2[names[1]]
    pairs_frame = pandas.read_csv(pair_list,sep=',',header=None)
    pairs = pairs_frame.values
    
    for i in np.arange(len(pairs)):
        seq1 = csv1_seq[pairs[i,0]]
        seq2 = csv2_seq[pairs[i,1]]
        new_frame = pandas.DataFrame([seq1,seq2])
        new_frame.index = [names[0],names[1]]
        new_frame.columns = [ label + ' ' + str(i)]
        if i == 0:
            final_pairs = new_frame
        else:
            final_pairs = pandas.concat([final_pairs,new_frame],axis=1)
    return(final_pairs)

########## So this is a first draft of loading metadata into AIMS ##########
# For now, this is specifically for loading in distance metadata for the
# Dpr-DIP application. Need to create a standard formatting for future use
def load_metadata(file):
    dist_load = pandas.read_csv(file,sep = '|',header=1)
    dist_load.columns = ['Query Chain', 'Interacting Chains', 'Dist', 'AtomClasses']
    first = True
    for i in dist_load['Interacting Chains'].values:
        # So it seems like there are always 3 spaces before the numbers
        skip = False
        hold = i.find('   ')
        count = 3
        zz = i[hold+count:]
        hold2 = zz.find(' ')
        while hold2 == 0:
            count += 1
            zz = i[hold+count:]
            hold2 = zz.find(' ')
        if first:
            interact_num = int(zz[:hold2])
            first = False
        else:
            interact_num = np.hstack((interact_num, int(zz[:hold2])))
    first = True
    for i in dist_load['Query Chain'].values:
        # So it seems like there are always 3 spaces before the numbers
        skip = False
        hold = i.find('   ')
        count = 3
        zz = i[hold+count:]
        hold2 = zz.find(' ')
        while hold2 == 0:
            count += 1
            zz = i[hold+count:]
            hold2 = zz.find(' ')
        if first:
            query_num = int(zz[:hold2])
            first = False
        else:
            query_num = np.hstack((query_num, int(zz[:hold2])))
            
    dist_load['Query Chain'] = query_num
    dist_load['Interacting Chains'] = interact_num

    first = True
    for i in dist_load[['Query Chain','Interacting Chains']].drop_duplicates().values:
        res1 = dist_load[dist_load['Query Chain'] == i[0]]
        resf = res1[res1['Interacting Chains'] == i[1]]
        avg_dist = np.average(resf['Dist'])
        if first:
            # WE CAN CHANGE THIS TO EITHER MIN, MAX, OR AVG DIST
            input_vects = np.hstack((i,avg_dist))
            new_frame = pandas.DataFrame(input_vects)
            first = False
        else:
            # WE CAN CHANGE THIS TO EITHER MIN, MAX, OR AVG DIST
            input_vects = np.hstack((i,avg_dist))
            new_frame = pandas.concat([new_frame,pandas.DataFrame(input_vects)],axis=1)
    meta_dist = np.transpose(new_frame)
    meta_dist.columns = ['Dpr','DIP','Dist']

    # Need to define the different alignments for visually matching to alignment figure AND for incorporating
    # in some distance metadata for these...
    # Dpr first
    dpr_pdb_numbering = [33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,
                        0,0,0,0,0,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107]
    dpr_aln_numbering = [40,41,47,48,49,50,51,52,53,54,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,
                        0,0,0,0,0,113,114,115,116,117,118,120,121,122,123,124,125,126,127,128,129,130,131,132]
    # Then DIP
    dip_pdb_numbering = [140,141,142,143,144,145,146,147,148,149,150,151,152,153,154,155,156,157,158,159,160,
                        161,162,163,164,165,0,0,0,0,0,184,185,186,187,188,189,190,191,192,193,194,195,196,197,198,199,200,201,202,203,204]
    dip_aln_numbering = [33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,
                        0,0,0,0,0,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100]

    # And export it all:
    return(meta_dist,dpr_pdb_numbering,dpr_aln_numbering,dip_pdb_numbering,dip_aln_numbering)

def generate_score_df(seq,scores,dist,pairs,mols):
    first = True
    for i in np.arange(np.shape(seq)[1]):
        for j in pairs:
            res1 = np.transpose(seq)[mols[0]].values[i][int(j[0])]
            res2 = np.transpose(seq)[mols[1]].values[i][int(j[1])]
            
            score_temp1 = scores[i][int(j[0]),int(j[1])] * dist.values[int(j[0]),int(j[1])]
            score_temp2 = scores[i][int(j[0]),int(j[1])]
            
            vect = [res1,res2,score_temp1,score_temp2,int(j[0]),int(j[1]),i]
            if first:
                important_res = vect
                first = False
            else:
                important_res = np.vstack((important_res,vect))

    res_df = pandas.DataFrame(important_res); res_df.columns = ['Mol1_res','Mol2_res','Pair_Score','Unweight_Score','Res1','Res2','ID']
    return(res_df)

def reshape_scores(dataframe):
    first = True
    for i in dataframe['ID'].drop_duplicates():
        sin_prot_score = dataframe[dataframe['ID'] == i]['Pair_Score']
        if first:
            score_reshape = [float(a) for a in sin_prot_score.values]
            first = False
        else:
            score_reshape = np.vstack((score_reshape, [float(a) for a in sin_prot_score.values]))
    return(score_reshape)
    
def msa_maxdiff(frame, pairs, scores, top_X):
    avg_score_df = pandas.DataFrame(scores)
    pos_index = pandas.DataFrame(abs(avg_score_df[0] - avg_score_df[1])).sort_values(0).index[-top_X:]
    first = True
    for i in pairs[pos_index]:
        sub_df = frame[frame['Res1'] == str(int(i[0]))]
        fin_df = sub_df[sub_df['Res2'] == str(int(i[1]))]
        if first:
            TOP = fin_df
            first = False
        else:
            TOP = pandas.concat((TOP,fin_df))

    return(TOP)

# So the score_reshape vars are #pairs x #vectors, we want to change that to #pairs x 6
# Need to program in options to either count or accumulate some kind of average score
def general_score(dataframe,option = 'count'):
    first=True; scores = [-2.0, -1.0, -0.5, 0.0, 1.0, 2.0]
    for i in dataframe['ID'].drop_duplicates():
        sin_prot_score = dataframe[dataframe['ID'] == i]['Unweight_Score']
        sin_prot_weighted = dataframe[dataframe['ID'] == i]['Pair_Score']
        if first:
            first = False
            counter = [0,0,0,0,0,0]
            for j in np.arange(len(sin_prot_score)):
                for k in np.arange(len(scores)):
                    if float(sin_prot_score.values[j]) == scores[k]:
                        if option == 'count':
                            counter[k] = counter[k] + 1
                        elif option == 'sum':
                            counter[k] = counter[k] + float(sin_prot_weighted.values[j])
            counterF = counter
        else:
            counter = [0,0,0,0,0,0]
            for j in np.arange(len(sin_prot_score)):
                for k in np.arange(len(scores)):
                    if float(sin_prot_score.values[j]) == scores[k]:
                        if option == 'count':
                            counter[k] = counter[k]+1
                        elif option == 'sum':
                            counter[k] = counter[k] + float(sin_prot_weighted.values[j])
            counterF = np.vstack((counterF,counter))
    return(counterF)

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

def full_AA_freq(seq,norm='num_AA'):
    AA_freq_all = np.zeros((20))
    digram_all = np.zeros((20,20))
    AAs = 0; num_seq = 0
    datlen,looplen = np.shape(seq.values)
    for i in np.arange(datlen):
        for j in np.arange(looplen):
            if len(seq.values[i][j]) == 1:
                temp_seq = seq.values[i][j][0]
            else:
                temp_seq = seq.values[i][j]
            num_seq = num_seq + 1
            AAs = AAs + len(temp_seq)
            for loc in np.arange(len(temp_seq)):
                res1 = temp_seq[loc]
                if loc + 1 < len(temp_seq):
                    res2 = temp_seq[loc+1]
                else:
                    res2 = -1
                matched = False
                for mat_loc1 in np.arange(len(AA_key)):
                    if AA_key[mat_loc1] == res1:
                        AA_freq_all[mat_loc1] = AA_freq_all[mat_loc1] + 1
                        matched = True
                        if res2 != -1:
                            for mat_loc2 in np.arange(len(AA_key)):
                                if AA_key[mat_loc2] == res2:
                                    digram_all[mat_loc1,mat_loc2] = digram_all[mat_loc1,mat_loc2] + 1
                                    break
                    if matched:
                        break
    if norm == 'num_AA':
        AA_freq_all = AA_freq_all/AAs
        digram_all = digram_all/AAs
    elif norm == 'num_seq':
        AA_freq_all = AA_freq_all/num_seq
        digram_all = digram_all/num_seq
    return(AA_freq_all,digram_all)

#### Begin Changes Made Explicitly for MHC Germline Analysis Manuscript##############
####################################################################################
def labelIT(dlen,labels):
    for j in np.arange(dlen):
        if j ==0:
            tab = [labels+ str(j)];
            organ = [labels]
        else:
            tab = tab + [labels + str(j)]; 
            organ = organ + [labels]
    return(tab,organ)

# Get out specific regions of the MHC Helices
def get_mhcSub(mhc,seq_choice,multiOrg=False,input_loc=[]):
    # This variable is only important for multi-org analysis
    # Key for matching sequences to metadata down the line
    loc_choice = input_loc
    if mhc == 'classI':
        alpha1 = [55,56,59,60,63,66,67,70,74,77,80]
        pep_contact = [3,5,7,20,22,24,43,57,61,64,65,68,71,72,75,78,79,82,93,95,97,112,114,141,145,150,154,157,169]
        alpha2 = [143,144,147,148,149,152,153,156,160,161,164,165,167,168]
    elif mhc == 'classIIa':
        II_alpha_contacts = [51, 53, 55, 56, 58, 59, 62, 63, 66, 69, 70, 73, 74]
        alpha1 = II_alpha_contacts
        alpha2 = []
        pep_contact =[]
    elif mhc == 'classIIb':
        II_beta_contacts = [63, 66, 67, 70, 71, 73, 74, 76, 77, 80, 83, 84, 87, 88, 91, 92]
        alpha1 = II_beta_contacts
        alpha2 = []
        pep_contact = []


    first = True
    for i in np.arange(len(seq_choice)):
        sub_seq1 = ''
        for j in alpha1:
            let = seq_choice[i][j]
            sub_seq1 = sub_seq1 + let
        sub_seq2 = ''
        for j in alpha2:
            let = seq_choice[i][j]
            sub_seq2 = sub_seq2 + let
        sub_seq3 = ''
        for j in pep_contact:
            let = seq_choice[i][j]
            sub_seq3 = sub_seq3 + let
        # Remove incomplete sequences
        if sub_seq1.find('-') != -1 or sub_seq2.find('-') != -1 or sub_seq3.find('-') != -1:
            continue
        # Make sure we dont have empty entries for classII
        if mhc == 'classI':
            pre_seq = np.array([sub_seq1,sub_seq2,sub_seq3])
        else:
            pre_seq = sub_seq1
        # Make our final output
        if first:
            plot_seq = pre_seq
            first = False
            if multiOrg:
                seq_loc = loc_choice[i]
        else:
            plot_seq = np.vstack((plot_seq,pre_seq))
            if multiOrg:
                seq_loc = np.hstack((seq_loc,loc_choice[i]))
    if multiOrg:
        return(plot_seq,seq_loc)
    else:
        return(plot_seq)

# A fuction used for randomizing TCR-MHC pairing to calculate information theoertic quantitites
def randomize_tcr_mhc_pair(trv_mat_df,mhc_mat_df,multiOrg = False,mhc_orgs=['Human']):
    first = True
            
    for org in mhc_orgs:
        if multiOrg:
            trv_sub = trv_mat_df[trv_mat_df.index.str.contains(org)]
            mhc_sub = mhc_mat_df[mhc_mat_df.index.str.contains(org)]
        else:
            trv_sub = trv_mat_df
            mhc_sub = mhc_mat_df

        # Some of our organism subsets dont have representatives in every subset
        if len(mhc_sub) == 0:
            continue
        if len(trv_sub) == 0:
            continue

        if len(mhc_sub) > len(trv_sub):
            mhc_re = []
            trv_re = trv_sub
            mhc_pre=mhc_sub.sample(frac=1).reset_index(drop=True)
            smallLen = len(trv_re)
            # We take advantage of the fact that "sample" randomly shuffles
            # the sequences in the subset of interest
            mhc_re = mhc_pre[:smallLen]
            
        elif len(mhc_sub) < len(trv_sub):
            trv_re = []
            mhc_re = mhc_sub
            trv_pre = trv_sub.sample(frac=1).reset_index(drop=True)
            smallLen = len(mhc_re)
            # We take advantage of the fact that "sample" randomly shuffles
            # the sequences in the subset of interest
            trv_re = trv_pre[:smallLen]

        # So at this point, they should already be equal values
        MHCseq = mhc_re.values
        TRVseq = trv_re.values
        # So in this version of the script, we take the SAME structural subfeatures and
        # apply a single sequence to each one. This feels the most rigorous.
        for i in np.arange(smallLen):
            tcr = TRVseq[i,:]
            mhc = MHCseq[i]

            if first:
                matF = np.hstack([tcr,mhc])
                first = False
            else:
                mat_pre = np.hstack([tcr,mhc])
                matF = np.vstack((matF,mat_pre))
    return(matF)

# This script is for actively breaking down the interaction scores for TCR-MHC interactions
def get_byAllele_scores(SCORES,tcr_names,hla_names,mhc_type = 'classI',ScoreAlpha1=True,ScoreAlpha2 = False,
                           len_weight=False,score_weight=False,clash=False,BADclash = False):
    # In this, we want to take the scores and identify WHERE
    # The positive interactions are coming from.
    # This is for if you want a subset
    if mhc_type == 'classI':
        alpha1 = [55,56,59,60,63,66,67,70,74,77,80]
        alpha2 = [143,144,147,148,149,152,153,156,160,161,164,165,167,168]
    elif mhc_type == 'classII_alpha':
        II_alpha_contacts = [51, 53, 55, 56, 58, 59, 62, 63, 66, 69, 70, 73, 74]
        alpha1 = II_alpha_contacts
        alpha2 = []
    elif mhc_type == 'classII_beta':
        II_beta_contacts = [63, 66, 67, 70, 71, 73, 74, 76, 77, 80, 83, 84, 87, 88, 91, 92]
        alpha1 = II_beta_contacts
        alpha2 = []
    
    output_hist = np.zeros((len(tcr_names),len(hla_names)))
    output_df = pandas.DataFrame(output_hist)
    output_df.columns = [i[0] for i in hla_names]
    output_df.index = [i for i in tcr_names]
    First = True
    for i in hla_names:
        for j in tcr_names:
            if First:
                fin_name = [i[0],j]
                First = False
            else:
                fin_name = np.vstack((fin_name,[i[0],j]))
    
    noMatch = []
    first = True
    for seq in np.arange(len(SCORES)):
        lochla = fin_name[seq][0]
        loctcr = fin_name[seq][1]
        test_score = SCORES[seq]
        len1, len2 = np.shape(test_score)
        save_coords = []
        for i in np.arange(len1):
            # Break when you are reaching past the end
            if i + 2 >= len1:
                break
            if ScoreAlpha1:
                #do it for the alpha1 loop
                for j in np.arange(len(alpha1)):
                    # Break when you are past the end
                    if j + 2 >= len(alpha1):
                        break
                    # Looking for contiguous regions
                    # of some positive interactions
                    # can change criteria later
                    test1 = test_score[i,alpha1[j]]
                    test2 = test_score[i+1,alpha1[j+1]]
                    test3 = test_score[i+2,alpha1[j+2]]
                    if clash:
                        tf1 = (test1 < 0)
                        tf2 = (test2 < 0)
                        tf3 = (test3 < 0)
                    elif BADclash:
                        # only test tf1
                        tf1 = (test1 == -2)
                        tf2 = True
                        tf3 = True
                    else:
                        tf1 = (test1 > 0)
                        tf2 = (test2 > 0)
                        tf3 = (test3 > 0)
                    if tf1 & tf2 & tf3:
                        # Save only the middle coords
                        save_coords = save_coords + [[i+1,alpha1[j],alpha1[j+1],alpha1[j+2]]]

                        # Add to a counter using this complicated line
                        if score_weight:
                            output_df.loc[loctcr,lochla] = output_df.loc[loctcr,lochla] + (test1+test2+test3)/3
                        elif BADclash:
                            output_df.loc[loctcr,lochla] = output_df.loc[loctcr,lochla] + test1
                        else:
                            output_df.loc[loctcr,lochla] = output_df.loc[loctcr,lochla] + 1
            if ScoreAlpha2:
                # Do it again for the alpha2 loop
                for j in np.arange(len(alpha2)):
                    # Break when you are past the end
                    if j + 2 >= len(alpha2):
                        break
                    # Looking for contiguous regions
                    # of some positive interactions
                    # can change criteria later
                    test1 = test_score[i,alpha2[j]]
                    test2 = test_score[i+1,alpha2[j+1]]
                    test3 = test_score[i+2,alpha2[j+2]]
                    if clash:
                        tf1 = (test1 < 0)
                        tf2 = (test2 < 0)
                        tf3 = (test3 < 0)
                    elif BADclash:
                        # only test tf1
                        tf1 = (test1 == -2)
                        tf2 = True
                        tf3 = True
                    else:
                        tf1 = (test1 > 0)
                        tf2 = (test2 > 0)
                        tf3 = (test3 > 0)
                    if tf1 & tf2 & tf3:
                        # Save only the middle coords
                        save_coords = save_coords + [[i+1,alpha2[j],alpha2[j+1],alpha2[j+2]]]

                        # Add to the counter above
                        if score_weight:
                            output_df.loc[loctcr,lochla] = output_df.loc[loctcr,lochla] + (test1+test2+test3)/3
                        elif BADclash:
                            output_df.loc[loctcr,lochla] = output_df.loc[loctcr,lochla] + test1
                        else:
                            output_df.loc[loctcr,lochla] = output_df.loc[loctcr,lochla] + 1
                            
        # Before we move on, if you want to, scale all these by the CDR lengths
        # Obviously based on this simple counting metric, longer CDRs will have
        # higher interaction scores by default
        if len_weight:
            # len1 should give the length of the CDR of interest
            output_df.loc[loctcr,lochla] = output_df.loc[loctcr,lochla]/len1
        if len(save_coords) == 0:
            # There were no contiguous matches for a bunch of these sequences!
            # Save the list of these sequences to be sure
            noMatch = noMatch + [seq]
            continue
            
        if first:
            final_coords = save_coords
            first = False
        else:
            final_coords = np.vstack((final_coords,save_coords))
    if first:
        print('no matches at all')
        return(1,2,3)
    else:
        return(output_df,final_coords,noMatch)

# Make a script to identify regions of contiguous positive interactions
# Start out first looking for ONLY strings of positive interactions
# Given the curvature of TCRs, lets only look in strings of 3

def get_byRes_scores(SCORES,mhc_type = 'classI',scoreWeight = False,lenWeight=False):
    # This is for if you want FULL interact scores
    #test_score = SCORES[0]
    #len1, len2 = np.shape(test_score)
    #bind_loc = np.zeros(len2)
    # This is for if you want a subset
    if mhc_type == 'classI':
        alpha1 = [55,56,59,60,63,66,67,70,74,77,80]
        alpha2 = [143,144,147,148,149,152,153,156,160,161,164,165,167,168]
    elif mhc_type == 'classII_alpha':
        II_alpha_contacts = [51, 53, 55, 56, 58, 59, 62, 63, 66, 69, 70, 73, 74]
        alpha1 = II_alpha_contacts
        alpha2 = []
    elif mhc_type == 'classII_beta':
        II_beta_contacts = [63, 66, 67, 70, 71, 73, 74, 76, 77, 80, 83, 84, 87, 88, 91, 92]
        alpha1 = II_beta_contacts
        alpha2 = []
    
    hist_alpha1 = np.zeros(len(alpha1))
    hist_alpha2 = np.zeros(len(alpha2))

    noMatch = []
    first = True
    for seq in np.arange(len(SCORES)):
        test_score = SCORES[seq]
        len1, len2 = np.shape(test_score)
        save_coords = []
        for i in np.arange(len1):
            # Break when you are reaching past the end
            if i + 2 >= len1:
                break
            #do it for the alpha1 loop
            for j in np.arange(len(alpha1)):
                # Break when you are past the end
                if j + 2 >= len(alpha1):
                    break
                # Looking for contiguous regions
                # of some positive interactions
                # can change criteria later
                test1 = (test_score[i,alpha1[j]] > 0)
                test2 = (test_score[i+1,alpha1[j+1]] > 0)
                test3 = (test_score[i+2,alpha1[j+2]] > 0)
                if test1 & test2 & test3:
                    # Save only the middle coords
                    save_coords = save_coords + [[i+1,alpha1[j],alpha1[j+1],alpha1[j+2]]]
                    
                    # Add to a histogram
                    if scoreWeight:
                        mult1 = test_score[i,alpha1[j]]
                        mult2 = test_score[i+1,alpha1[j+1]]
                        mult3 = test_score[i+2,alpha1[j+2]]
                    else:
                        mult1 = 1; mult2 = 1; mult3 = 1
                        
                    hist_alpha1[j] = hist_alpha1[j] + 1 * mult1
                    hist_alpha1[j+1] = hist_alpha1[j+1] + 1 * mult2
                    hist_alpha1[j+2] = hist_alpha1[j+2] + 1 * mult3
                    
            # Do it again for the alpha2 loop
            for j in np.arange(len(alpha2)):
                # Break when you are past the end
                if j + 2 >= len(alpha2):
                    break
                # Looking for contiguous regions
                # of some positive interactions
                # can change criteria later
                test1 = (test_score[i,alpha2[j]] > 0)
                test2 = (test_score[i+1,alpha2[j+1]] > 0)
                test3 = (test_score[i+2,alpha2[j+2]] > 0)
                if test1 & test2 & test3:
                    # Save only the middle coords
                    save_coords = save_coords + [[i+1,alpha2[j],alpha2[j+1],alpha2[j+2]]]
                    
                    # Add to a histogram
                    if scoreWeight:
                        mult1 = test_score[i,alpha2[j]]
                        mult2 = test_score[i+1,alpha2[j+1]]
                        mult3 = test_score[i+2,alpha2[j+2]]
                    else:
                        mult1 = 1; mult2 = 1; mult3 = 1
                        
                    hist_alpha2[j] = hist_alpha2[j] + 1 * mult1
                    hist_alpha2[j+1] = hist_alpha2[j+1] + 1 * mult2
                    hist_alpha2[j+2] = hist_alpha2[j+2] + 1 * mult3
                    
        if len(save_coords) == 0:
            # There were no contiguous matches for a bunch of these sequences!
            # Save the list of these sequences to be sure
            noMatch = noMatch + [seq]
            continue
            
        if first:
            final_coords = save_coords
            first = False
        else:
            final_coords = np.vstack((final_coords,save_coords))
    if lenWeight:
        fin_data = np.hstack((hist_alpha1,hist_alpha2))/len(SCORES)
    else:
        fin_data = np.hstack((hist_alpha1,hist_alpha2))
    return(fin_data)

# convert matF back into sequences...
# So the input to matF should be a single numeric row of the AIMS matrix
def decode_mat(matF,num_key_AA,key_AA):
    hold_str = []
    for i in matF:
        if i == 0:
            hold_str = hold_str + ['-']
        else:
            for j in num_key_AA:
                if i == j:
                    hold_str = hold_str + [key_AA[j-1]]
    final = ''.join(hold_str)
    return(final)

# This script is necessary to pull out the CDR 1 and 2 loop from gene names
def pull_cdr_1_2(gene_list,chain='trav',organism='Human'):
    org = organism.lower()
    trv = chain.lower()
    fin_trv, trv_name_pre = aimsLoad.Ig_loader(datPath+'app_data/germline_data/'+trv+'_'+org+'_cdrs.csv','tcr',loops=3,return_index = True)    
    fin_trv.columns = trv_name_pre
    
    fin_cdrs = []
    for i in gene_list.values:
        # Should probably include some sort of failsafe in case 
        try:
            cdrs = [fin_trv[i].values]
        except:
            # So in dealing with MiSeq data, these may be pseudogenes frequently.
            #print("Gene " +i+" not found")
            cdrs = [['','','']]
        
        fin_cdrs = fin_cdrs + cdrs
    
    return(pandas.DataFrame(fin_cdrs))

# Turn a single column dataframe of metadata into a numeric dataframe
# Currently only for tags/integers, but will be edited in the future for numeric float data
def encode_meta(metadat):
    temp_df = metadat
    uniq_labels = temp_df.drop_duplicates().values
    dict_f = {}
    for i in np.arange(len(uniq_labels)):
        dict_temp = {uniq_labels[i][0]:i}
        dict_f.update(dict_temp)
    encoded=temp_df.replace({0:dict_f})
    return(encoded)

# For returning cluster purity of biophysical clusters with metadata
# Do it a slightly different way than prev. version and return a matrix of fractions
# Way more thorough and less complicated than trying to quantify purity.
def calc_cluster_purity(final_breakdown,meta_name,do_stats = True,num_reps=1000):
    clusters = np.sort(final_breakdown['cluster'].drop_duplicates())
    metas = final_breakdown[meta_name].drop_duplicates()
    temp_mat = np.zeros((len(clusters),len(metas)))
    pure_df = pandas.DataFrame(temp_mat)
    pure_df.index = clusters
    pure_df.columns = metas.values
    cluster_size = []
    for i in np.arange(len(clusters)):
        cluster_order = np.sort(final_breakdown['cluster'].drop_duplicates())[i]
        sub_clust = final_breakdown[final_breakdown['cluster'] == cluster_order]
        cluster_size = cluster_size + [len(sub_clust)]
        for j in np.arange(len(metas)):
            num_ant = len(sub_clust[sub_clust[meta_name]==metas.values[j]])
            pure_df.iloc[i,j] = num_ant/len(sub_clust)
    if do_stats:
        # Try to do some stats off the pure_df
        # We're going to stick with the permutation test, 
        # But should likely look into Fisher's exact test at some point
        test_mat = pandas.DataFrame(np.zeros((len(clusters),len(metas))))
        for rep in np.arange(num_reps):
            permute_temp = final_breakdown[meta_name].values
            permute_meta = resample(permute_temp,replace=False)
            
            # randomly assign "clusters":
            prev_clustered = 0 
            for ri in np.arange(len(cluster_size)):
                permute_clust = permute_meta[prev_clustered:cluster_size[ri]]
                for rj in np.arange(len(metas)):
                    num_ant = len(permute_clust[permute_clust==metas.values[rj]])
                    test_mat.iloc[ri,rj] = num_ant/len(permute_clust)

            test_mat.index = pure_df.index
            test_mat.columns = pure_df.columns

            z_temp = test_mat >= pure_df
            z_check = z_temp.astype(int)

            if rep == 0:
                num_sig = z_check
            else:
                num_sig += z_check

        p_mat = (1+num_sig)/num_reps
        return(pure_df,p_mat)
    else:
        return(pure_df)

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

# This function is really only in here to make the notebook look a little cleaner
# Changing this so that 2D always comes first...
def get_plotdefs(clust_show,proj_show,chosen_map1,chosen_map2,leg1,leg2):
    if clust_show.lower() == 'both':
        if proj_show.lower() == 'both':
            fig3d = pl.figure(figsize = (25, 20))
            plotloc = [221,222,223,224]
            plottype = ['2d','3d','2d','3d']
            dattype = ['clust','clust','meta','meta']
            plotem=[chosen_map1,chosen_map1,chosen_map2,chosen_map2]
            legends = [leg1,leg1,leg2,leg2]
        elif proj_show.lower() == '2d':
            fig3d = pl.figure(figsize = (20, 10))
            plotloc = [121,122]
            plottype=['2d','2d']
            dattype = ['clust','meta']
            plotem = [chosen_map1, chosen_map2]
            legends = [leg1, leg2]
        elif proj_show.lower() == '3d':
            fig3d = pl.figure(figsize = (20, 10))
            plotloc = [121,122]
            plottype=['3d','3d']
            dattype = ['clust','meta']
            plotem = [chosen_map1, chosen_map2]
            legends = [leg1, leg2]
    elif clust_show.lower() == 'clusters':
        if proj_show.lower() == 'both':
            fig3d = pl.figure(figsize = (20, 10))
            plotloc = [121,122]
            plottype = ['2d','3d']
            dattype = ['clust','clust']
            plotem=[chosen_map1,chosen_map1]
            legends = [leg1,leg1]
        elif proj_show.lower() == '2d':
            fig3d = pl.figure(figsize = (10, 10))
            plotloc = [111]
            plottype=['2d']
            dattype = ['clust']
            plotem = [chosen_map1]
            legends = [leg1]
        elif proj_show.lower() == '3d':
            fig3d = pl.figure(figsize = (10, 10))
            plotloc = [111]
            plottype=['3d']
            dattype = ['clust']
            plotem = [chosen_map1]
            legends = [leg1]
    elif clust_show.lower() == 'metadata':
        if proj_show.lower() == 'both':
            fig3d = pl.figure(figsize = (20, 10))
            plotloc = [121,122]
            plottype = ['2d','3d']
            dattype = ['meta','meta']
            plotem=[chosen_map2,chosen_map2]
            legends = [leg2,leg2]
        elif proj_show.lower() == '2d':
            fig3d = pl.figure(figsize = (10, 10))
            plotloc = [111]
            plottype=['2d']
            dattype = ['meta']
            plotem = [chosen_map2]
            legends = [leg2]
        elif proj_show.lower() == '3d':
            fig3d = pl.figure(figsize = (10, 10))
            plotloc = [111]
            dattype = ['meta']
            plottype=['3d']
            plotem = [chosen_map2]
            legends = [leg2]

    return(fig3d,plotloc,plottype,plotem,legends,dattype)

# Finally, a hard-coded in way to calculate the AIMS distances for direct comparison to TCRdist
def calc_AIMSdist(seqSet1, seqSet2='',matrix='',align='center',normalize='msuv',special='',pad=6,form = ''):
    # If the user provides two sequence sets, concatenate them to start
    if len(seqSet2) == 0:
        dsetF = seqSet1
    else:
        if seqSet1.equals(seqSet2):
            print("ERROR: If trying to compare internal distances, just do not provide a second dataframe")
            return()
        dsetF = pandas.concat([seqSet1,seqSet2],axis=0)
    
    if len(matrix) == 0:
        bigass = classy.get_bigass_matrix(np.transpose(dsetF.values), alignment = align, norm = normalize,special=special,bulge_pad=pad )
        fin_mat = pandas.DataFrame(bigass,index=dsetF.index)
    else:
        # This assumes a user-provided matrix that has matched entries
        # plus, an index which matches the subset incdices.
        
        # You can use this user-supplied matrix to look at different distances
        # like from the parsed matrix. Cant generate parsed_mat within this
        # function because it would introduce too much variability
        fin_mat = matrix.loc[dsetF.index]
    
    if len(seqSet2) == 0:
        seqs1 = seqSet1
        seqs2 = seqSet1
    else:
        seqs1 = seqSet1
        seqs2 = seqSet2
        
    numSeq1 = len(seqs1); numSeq2 = len(seqs2)
    dist_calc = np.zeros((numSeq1,numSeq2))
    for i in np.arange(numSeq1):
        for j in np.arange(numSeq2):
            vect1 = fin_mat.loc[seqs1.index[i]].values
            vect2 = fin_mat.loc[seqs2.index[j]].values
            dist_calc[i,j] = np.sqrt(sum((vect1 - vect2)**2))

    if form == 'tri':
        if len(seqs1) != len(seqs2):
            print("ERROR: Cannont return upper triangle of matrix if matrix is not symmetric")
            return(dist_calc)
        return(dist_calc[np.triu_indices(len(dist_calc))])
    else:
        return(dist_calc)

# Finally, a way to calculate statistics for the AIMS analysis
# Also note, I was able to very nicely test that the p-value converges as num_rep->inf
# Typically pretty fast. I think it's supposed to converge as num_rep -> N
def do_statistics(data1,data2,num_reps = 1000,test='median',multi_test='none',alpha=0.05,test_func = [],func_val=0):
    # Should probably have a failsafe to make sure that we're looking across the proper axes
    # prop_axis should include the number of samples. We assume it should be axis 0
    # Ideally, we wont have many situations where the number of samples is exactly equal...
    if len(np.shape(data1)) == 1:
        prop_axis = 0; test_axis = 0
    elif np.shape(data1)[1] == np.shape(data2)[1]:
        prop_axis = 0; test_axis = 1
    elif np.shape(data1)[0] == np.shape(data2)[0]:
        prop_axis = 1; test_axis = 0
    elif np.shape(data1) == np.shape(data2):
        prop_axis = 1; test_axis = 0


    # Alright lets finally get around to calculating stats.
    len1 = np.shape(data1)[prop_axis]
    preDat = np.concatenate((data1,data2),axis=prop_axis)

    if test.lower()=='average':
        z0=np.average(data1,axis=prop_axis)-np.average(data2,axis=prop_axis)
    elif test.lower()=='median':
        z0=np.median(data1,axis=prop_axis)-np.median(data2,axis=prop_axis)
    elif test.lower()=='diff':
        # In some cases (namely MI and Entropy) we might just want to look at a simple difference
        z0 = data1-data2
    elif test.lower()=='function':
        # Allow for custom functions (or MI/Shannon entropy) to calc sigFigs
        if test_func == full_AA_freq:
            temp1 = test_func(data1)[func_val]
            temp2 = test_func(data2)[func_val]
        else:
            temp1 = test_func(np.transpose(np.array(data1)))[func_val]
            temp2 = test_func(np.transpose(np.array(data2)))[func_val]
        z0 = temp1 - temp2
    # add more tests in the future?

    if type(z0) == np.float64:
        num_sig = 0

    for rep in np.arange(num_reps):
        # Turns out that np.random.shuffle mixes up the entries in place
        tempAll = preDat
        if prop_axis == 0:
            allDat = resample(tempAll,replace=False)
            re_dat1 = allDat[0:len1]
            re_dat2 = allDat[len1:]
        else:
            allDat = np.transpose(resample(np.transpose(tempAll),replace=False))
            re_dat1 = allDat[:,0:len1]
            re_dat2 = allDat[:,len1:]

        if test.lower() == 'average':
            z = np.average(re_dat1,axis=prop_axis) - np.average(re_dat2,axis=prop_axis)
        elif test.lower() == 'median':
            z = np.median(re_dat1,axis=prop_axis) - np.median(re_dat2,axis=prop_axis)
        elif test.lower() == 'diff':
            z = re_dat1 - re_dat2
        elif test.lower() == 'function':
            if test_func == full_AA_freq:
                temp1 = test_func(pandas.DataFrame(re_dat1))[func_val]
                temp2 = test_func(pandas.DataFrame(re_dat2))[func_val]
            else:
                temp1 = test_func(np.transpose(np.array(re_dat1)))[func_val]
                temp2 = test_func(np.transpose(np.array(re_dat2)))[func_val]
            z = temp1 - temp2
        
        if type(z0) == np.float64:
            if z**2 >= z0**2:
                num_sig += 1
        else:
            # It turns out this works pretty nicely, doing
            # element-wise comparisons of the data
            z_temp = z**2 >= z0**2
            z_check = z_temp.astype(int)

            if rep == 0:
                num_sig = z_check
            else:
                num_sig += z_check
    
    p = (num_sig+1)/(num_reps+1)

    # Add in a module for multiple test correction.
    # Multi-test correction isn't about correcting p-values
    # it is about correcting what is considered stat-sig.
    if multi_test.lower() == 'none':
        return(p)
    # Have a failsafe in case someone mistakenly tries to multi-test correct a float
    if type(p) == np.float64:
        return(p)
    else:
        stat_sig = []
        if multi_test.lower() == 'bonferroni' or multi_test == 1:
            for tester in p:
                if tester < alpha/np.shape(data1)[test_axis]:
                    stat_sig = stat_sig + ['*']
                else:
                    stat_sig = stat_sig + ['ns']
            return(p,stat_sig)
        elif multi_test.lower() == 'benjamini-hochberg' or multi_test == 2:
            # first, sort the p-values in ascending order:
            # Need to be a little bit smarter since we're ordering the p-values
            # Make sure we relate back to the actual order of the data...
            pre_sort = pandas.DataFrame(p).sort_values(0)
            p_ordered = pre_sort.sort_valeus(0)
            p_loc = np.array(pre_sort.sort_values(0).index)

            stat_sig = ['']*len(p_ordered)
            for k in np.arange(len(p_ordered)):
                if p_ordered[k] < k/np.shape(data1)[test_axis]*alpha:
                    stat_sig[p_loc[k]] = ['*']
                else:
                    stat_sig[p_loc[k]] = ['ns']
            return(p,stat_sig)
        else:
            return('ERROR: Bad Test Correction Variable')

# Generate a simulated repertoire for response to reviewers...
# order of operations here:
# First: Select a V- and J- pairing
# Second: Determine #AAs left
# Third: Determine constrained combinatorics of added and removed AAs (all possible)
# Fourth: Only select those options that satisfy length distribution.
# Fifth: Select added amino acids based upon supplied probabilities, stich sequences
#### Simulation Assumptions ################################
#- 3 or 6 nucleotide deletions (remove 1-2 amino acids) randomly from the end of V- and start of J-
#- So IMGT ALWAYS removes these p-nucleotides from their protein displays... so go off of direct deletion values from Murugan et al.
#- 3 or 6 nucleotide insertions (add 1-2 amino acids) randomly to the end of V- and start of J-
#- LASTLY, maybe another 2-3 amino acids added for the D-segment... Although this is also a bit wonky...
#- Basically the D segments encode 4 AAs, but also get deletions from BOTH ends.
#- I will NOT be going off of the distribution of D-segments, although you could conceivably do so
def gen_sim_repertoire(numSeq=100,adic={},lens=14):
    from sklearn.utils import resample
    import random
    import itertools

    trav = pandas.read_csv(datPath+'app_data/germline_data/trav_human_cdrs.csv')
    trbv = pandas.read_csv(datPath+'app_data/germline_data/trbv_human_cdrs.csv')
    traj = pandas.read_csv(datPath+'app_data/germline_data/traj_human_cdrs.csv')
    trbj = pandas.read_csv(datPath+'app_data/germline_data/trbj_human_cdrs.csv')

    for rep in np.arange(numSeq):
        selV = random.randrange(0,len(trbv))
        selJ = random.randrange(0,len(trbj))

        myV = trbv.iloc[selV]['cdr3']
        myJ = trbj.iloc[selJ]['cdr3']

        del_opts = list(itertools.product([0,1,2],repeat=2))
        del_sel = del_opts[random.randrange(0,9)]

        if del_sel[0] == 0:
            finV = myV
        else:
            finV = myV[:-del_sel[0]]
        if del_sel[1] == 0:
            finJ = myJ
        else:
            finJ = myJ[del_sel[1]:]

        germ_len = len(finV) + len(finJ)

        # get the number of AAs needed left depending
        # on whether you want a range of lengths or just one
        if type(lens) == int:
            addAA = lens - germ_len 
        else:
            addAA_opts = np.array(lens)-germ_len
            sel_add = random.randrange(0,len(addAA_opts))
            addAA = addAA_opts[sel_add]

        # Alright now we need to decide which amino acids we're going to add...
        # Ok how I should probably do this is have my input be a dictionary, use that to find the 
        # location of the given weight (i.e. W is position 0), and then use the dict itself to change the
        # value at this position...
        rev_altered = 'WFMLIVPYHAGSTDECNQRK'
        aa_prob_weight = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]

        for i in adic.keys():
            aa_loc = rev_altered.find(i)
            aa_prob_weight[aa_loc] = adic[i]

        fin_AAs = []
        for i in np.arange(len(aa_prob_weight)):
            fin_AAs = fin_AAs + [rev_altered[i]*aa_prob_weight[i]]

        added = []
        for i in np.arange(addAA):
            added = added + [random.choice(''.join(fin_AAs))]

        final_seq = finV+''.join(added)+finJ

        if rep == 0:
            repertoire = final_seq
        else:
            repertoire = np.vstack((repertoire,final_seq)) 

    return(repertoire)

def prep_distCalc(sorted_seqs,splitSize=1000):
    import itertools
    zz = np.arange(0,len(sorted_seqs)+1,1000)

    s1 = np.arange(0,len(sorted_seqs),splitSize)
    s2 = np.arange(splitSize,len(sorted_seqs),splitSize)
    if len(s1) != len(s2):
        if len(s1) > len(s2):
            s2 = np.hstack((s2,len(sorted_seqs)))
        else:
            s1 = np.hstack((s1,len(sorted_seqs)))
    final = np.transpose(np.vstack((s1,s2)))

    xx = list(itertools.combinations_with_replacement(final,2))
    return(xx)

def get_distClusts(dists,metadat,max_d=5):
    # max_d determines where to draw your cutoff for what defines a "cluster"
    # Might need to tweak a little bit depending on if you're using AIMS or TCRdist
    from scipy.cluster.hierarchy import dendrogram, linkage
    from scipy.spatial.distance import squareform
    from scipy.cluster.hierarchy import fcluster
    
    mat = dists
    square_dist = squareform(mat)
    linkage_matrix = linkage(square_dist, "single")
    # This is only for *showing* the last N clusters
    #dend2 = dendrogram(linkage_matrix,truncate_mode='lastp',p=50)  # show only the last p merged clusters)

    clusters = fcluster(linkage_matrix, max_d, criterion='distance')
    aimsDist_df = pandas.DataFrame(clusters)
    aimsDist_final = pandas.concat([metadat,aimsDist_df],axis=1)
    aimsDist_final.columns = [['meta','Distcluster']]

    return(aimsDist_final)