import numpy as np
import pandas
import matplotlib.pyplot as pl
import math
import matplotlib as mpl
from matplotlib import cm

# Define some initial stuff and import analysis functions:
#AA_key_old=['A','G','L','M','F','W','K','Q','E','S','P','V','I','C','Y','H','R','N','D','T']
AA_key=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

# So we've got 46 orthogonal (or at least not super correlated)
# dimensions. Add them in to the matrix
# From "Hot spot prediction in protein-protein interactions by an ensemble system"
# Liu et. al. BMC Systems Biology
newnew=pandas.read_csv('app_data/new_props')
oldold=pandas.read_csv('app_data/old_props')

# Again, ugly to hard code in the number of properties (62) but 
# For now no harm no foul
properties=np.zeros((62,20))
for i in np.arange(len(AA_key)):
    properties[0:16,i]=oldold[AA_key[i]]
    properties[16:,i]=newnew[AA_key[i]]

AA_num_key_new=properties[1]
AA_num_key=np.arange(20)+1

def get_sequence_dimension(re_poly):
    num_loops,num_clones=np.shape(re_poly)
    for i in np.arange(num_loops):
        if i == 0:
            max_len=len(re_poly[i,0])
        elif i <= 5:
            max_len=np.vstack((max_len,len(re_poly[i,0])))
        for j in np.arange(num_clones):
            if i == 0:
                if len(re_poly[i,j]) > max_len:
                    max_len = len(re_poly[i,j])
            else:
                if len(re_poly[i,j]) > max_len[i]:
                    max_len[i] = len(re_poly[i,j])
    max_len=max_len+3 # Add 3 here so there's a more pronounced space between loops
    ## SO NOW MAX LEN SHOULD HAVE THE MAXIMUM LENGTH OF EACH CDR LOOP ##
    if num_loops == 1:
        sequence_dim = max_len
    else:
        sequence_dim = int(sum(max_len))
    return(max_len,sequence_dim)

# NOTE, manuscript_arrange=False MUST be selected to run MHC analysis
# I used to re-arrange the CDR loops for a more position-accurate
# representation. In the current analysis, this isn't necessary.
def gen_tcr_matrix(pre_poly,key=AA_num_key_new,binary=False,
pre_mono=[],giveSize=[],return_Size=False,manuscript_arrange=False,
alignment = 'center'):
    # Do this so that 
    if giveSize == []:
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
            max_lenp,sequence_dim=get_sequence_dimension(pre_poly)
    else:
        max_lenp = giveSize
        if type(max_lenp) == int:
            sequence_dim = max_lenp
        else:
            sequence_dim = int(sum(max_lenp))
    final_poly=[] # initialize a variable
    # RE-ORGANIZE EVERYTHING SO IT LOOKS NICE IN MATRIX
    # But presumably, if you give a size, it's the size you want...
    if manuscript_arrange and giveSize == []:
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
                        bulge = leng - 8 # Minus 8 because we are conserving 4 AA on each end
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
                        if int(LoopMax) < int(10):
                            count=int((max_len[k]-leng)/2.0)+int(sum(max_len[:k]))
                        else:
                            count = int(0) + int(sum(max_len[:k]))
                            b_count = 0
                            bulge = leng - 8
                            start_end = True
                            bulge_cen = int((max_len[k] - bulge)/2.0) + loopBuff
                for m in re_poly[k][i]: # SO IS THIS ONE for bootstrapping
                    for j in range(len(key)):
                        if m==AA_key[j]:
                            # Can't be doing this for short sequences
                            if alignment.lower() == 'bulge' and int(LoopMax) > int(10):
                                if bulge < 0:
                                    if count > 3 + loopBuff:
                                        if start_end:
                                            start_end = False
                                            # So here is where we skip to the very end
                                            count = LoopMax-4 + loopBuff
                                else:
                                    # Clearly don't need to do any of this stuff
                                    # if we DONT have a bulge, so make exception
                                    if count > 3 + loopBuff and count < LoopMax-4 +loopBuff:
                                        if b_count == 0:
                                            count = bulge_cen
                                        b_count += 1
                                    # This should then jump us over to the
                                    # end once we get past the bulk
                                    # ie, gaps should go around the bulk
                                    if b_count > bulge:
                                        if start_end:
                                            count = LoopMax-4 + loopBuff
                                            start_end = False
                            poly_PCA[i][int(count)]=key[j]
                            count += 1
                loop=loop+1
        if binary:
            # Unfortunate naming here for the binary case
            # but it is what it is...
            if final_poly == []:
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
    # We technically have 21 entries, 20 AAs (1-20) and spaces 0... Take a look at all of that
    AAs=np.arange(0,21)
    #print(AAs)

    for i in np.arange(clones):
        for j in np.arange((aas)):
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
    return(shannon_poly,poly_count)

def calculate_MI(poly_PCA):
    shannon_poly,poly_count=calculate_shannon(poly_PCA)
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
    properties=np.zeros((62,20))
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

    newnew=pandas.read_csv('app_data/new_props')
    oldold=pandas.read_csv('app_data/old_props')

    props=np.zeros((62,20))
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
def gen_tcr_matrixOLD(re_poly,max_len,key=AA_num_key_new):
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
def gen_1Chain_matrix(pre_poly,key=AA_num_key_new,binary=False,pre_mono=[],giveSize=[],return_Size=False):
    # Do this so that 
    if giveSize == []:
        if binary:
            # NEW ADDITION TO CLEAN THINGS UP A BIT #
            max_len1=get_sequence_dimension(pre_poly)[0]
            max_len2=get_sequence_dimension(pre_mono)[0]
            max_lenp=np.zeros(3)
            for i in np.arange(3):
                max_lenp[i]=max(max_len1[i],max_len2[i])

            sequence_dim = int(sum(max_lenp))
        else:
            max_lenp,sequence_dim=get_sequence_dimension(pre_poly)
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
            if final_poly == []:
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
def getBig(mono_PCA):
    # Try to maximize differences across the properties by looking at patterning...
    # Redifine "properties" because I was getting some weird errors...
    properties=np.zeros((62,20))
    for i in np.arange(len(AA_key)):
        properties[0:16,i]=oldold[AA_key[i]]
        properties[16:,i]=newnew[AA_key[i]]

    props = properties[2:]

    # Re-normalize the properties for use in the matrix...
    for i in np.arange(len(props)):
        props[i] = props[i]-np.average(props[i])
        props[i] = props[i]/np.linalg.norm(props[i])

    mono_pca_NEW = mono_PCA

    mono_dim1,mono_dim2=np.shape(mono_pca_NEW)

    # So this is where we should be able to do the averaging
    mono_prop_masks=np.zeros([len(props),mono_dim1,int(mono_dim2)])

    for i in np.arange(len(props)): # For all of our properties...
        for j in np.arange(mono_dim1): # for every clone
            for k in np.arange(int(mono_dim2)): # for every position
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

def gen_peptide_matrix(pre_pep1,key=AA_num_key_new,binary=False,pre_pep2=[]):
    # How many sequences do we have?
    numClone = len(pre_pep1)
    final_pep1=[] # initialize a variable
    # Allow for the possibility that you are doing a binary comparison
    for re_pep in [pre_pep1, pre_pep2]:
        numClone = len(re_pep[0])
        ### FOR NOW, HARD CODE A PEPTIDE SEQUENCE LEN MAX OF 18
        ### MIGHT NEED TO GET CREATIVE WITH THIS, MIGHT NEED TO
        ### ADAPT IT AS WE GO, JUST IN CASE.
        sequence_dim = 18
        pep_PCA=np.zeros([numClone,sequence_dim])

        for i in range(numClone): # For all of our polyreactive sequences...
            # this 'count' variable is how we center align... probably 
            # NOT what we want to do for the peptide analysis... except maybe
            # for the 'betwixt anchors' section
            count = int((sequence_dim - len(re_pep[0][i]))/2.0)

            # Check the work of Guillame et al. [PNAS 2018]
            # for some ideas on how to encode this matrix

            # This loop is where we put these rules for encoding
            # For now, let's make simple assumptions
            # (pos 2 and last pos are anchors)
            # BUT code up a "bulge" region that is centrally aligned
            for m in re_pep[0][i]:
                for j in range(len(key)):
                    if m==AA_key[j]:
                        pep_PCA[i][count]=key[j]
                        count=count+1
        if binary:
            # Unfortunate naming here for the binary case
            # but it is what it is...
            if final_pep1 == []:
                final_pep1 = pep_PCA
            else:
                final_pep2 = pep_PCA
        else:
            break
    
    if binary:
        return(final_pep1,final_pep2)
    else:
        return(pep_PCA)
###################################################
# Peptide stuff:
def gen_peptide_matrix(pre_pep1,key=AA_num_key_new,binary=False,pre_pep2=[]):
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
            if final_pep1 == []:
                final_pep1 = pep_PCA
            else:
                final_pep2 = pep_PCA
        else:
            break
    
    if binary:
        return(final_pep1,final_pep2)
    else:
        return(pep_PCA)
