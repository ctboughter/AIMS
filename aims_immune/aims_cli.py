#!/usr/bin/env python
# coding: utf-8

# # Welcome to the AIMS CLI
# # Section 0: Loading in Modules and Defining Figure Formatting
# This first cell is just loading in all of the necessary python modules (which you should have already installed) and defining figure font, size, etc.

# In[1]:
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as pl
from matplotlib import rcParams
from matplotlib import rc
from matplotlib.lines import Line2D
from mpl_toolkits import mplot3d
import pandas
import os
from aims_immune import aims_loader as aimsLoad
from aims_immune import aims_analysis as aims
from aims_immune import aims_classification as classy
import matplotlib.gridspec as gridspec
from sklearn.utils import resample
import argparse
import distutils

# This bit is for that figure formatting. Change font and font size if desired
font = {'family' : 'Arial',
        'weight' : 'bold',
        'size'   : 20}
COLOR = 'black'
rcParams['text.color'] = 'black'
rcParams['axes.labelcolor'] = COLOR
rcParams['xtick.color'] = COLOR
rcParams['ytick.color'] = COLOR
rc('font', **font)
# Lastly this custom colormap is for 
import matplotlib as mpl
upper = mpl.cm.jet(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
    lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
cmap = mpl.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

# Taking a crack at writing a parser
###############################################################################################
def main():
    parser = argparse.ArgumentParser(description = 'Run the AIMS Pipeline in a Single Command')
    parser.add_argument("-dd", "--datDir",help="Data Directory [string]",required=False,default='./',type=str)
    parser.add_argument("-od", "--outputDir",help="Output Directory [string]",required=False, default='AIMS_out',type=str)
    parser.add_argument("-m","--molecule",help="Molecule of interest (Ig, MSA, or Peptide) [string]",required=True,type=str)
    parser.add_argument("-f","--fileNames",help="Name of files [list of strings]",required=True,type=str,nargs='*')
    # No default for datName, need to account for this in the "main" instance
    parser.add_argument("-n","--datNames",help="Name of the datasets [list of strings]",required=False,default=[],nargs='*',type=str)
    parser.add_argument("-a","--align",help="Sequence alignment scheme [string]",required=False,default='center',type=str)
    # Need to fix this so the default is the shape of the input file
    parser.add_argument("-nl","--numLoop",help="Number of loops for Ig [int]",required=False,default=1,type=int)
    parser.add_argument("-dp","--dropDup",help="Drop duplicate sequences? [T/F]",required=False,default=False,type=distutils.util.strtobool)
    parser.add_argument("-p","--parallel",help="Turn on Parallel Processing? [T/F]",required=False,default=False,type=distutils.util.strtobool)
    parser.add_argument("-s","--subset",help="Take a subset of input data? [T/F]",required=False,default=False,type=distutils.util.strtobool)
    parser.add_argument("-ss","--subStart",help="Start points for data subset [list of ints]",required=False,default=[],type=int,nargs='*')
    parser.add_argument("-se","--subEnd",help="End points for data subset [list of ints]",required=False,default=[],type=int,nargs='*')
    parser.add_argument("-bp","--bulgePad",help="Padding for bulge format [int]",required=False,default=8,type=int)
    parser.add_argument("-np","--normProp",help="Normalize biophysical properties? [msuv,zscore, or 0to1]",required=False,default='msuv',type=str)
    parser.add_argument("-rn","--REnorm",help="Renormalize BPHYS mat by entropy? [T/F]",required=False,default=True,type=distutils.util.strtobool)
    parser.add_argument("-cd","--clustData",help="Data format for dim. red. and clustering [string]",required=False,default='parse',type=str)
    parser.add_argument("-pa","--projAlg",help="Algorithm for data projection [string]",required=False,default='pca',type=str)
    parser.add_argument("-us","--umapSeed",help="Random seed for UMAP projection [int]",required=False,default=[],type=int)
    parser.add_argument("-c","--clustAlg",help="Clustering algorithm of choice [string]",required=False,default='optics',type=str)
    parser.add_argument("-cs","--clustSize",help="Set cluster size for clustering algorithm [string]",required=False,default=10,type=int)
    # Need to come back to the metadata section soon
    parser.add_argument("-mf","--metaForm",help="Format of metadata [string]",required=False,default='category',type=str)
    parser.add_argument("-mn","--metaName",help="Name of the incorporated metadata [string]",required=False,default='meta',type=str)
    
    ##################################################
    # NEW FEATURES!!!! Woohoo
    parser.add_argument("-ds","--DOstats",help="Run statistics for relevant analysis? [T/F]",required=False,default=False,type=distutils.util.strtobool)
    parser.add_argument("-db","--DOboot",help="Run bootstrapping for relevant analysis? [T/F]",required=False,default=False,type=distutils.util.strtobool)
    parser.add_argument("-gd","--GETdist",help="Run AIMSdist calculations? [T/F]",required=False,default=False,type=distutils.util.strtobool)
    parser.add_argument("-pd","--PARdist",help="Make AIMSdist calc parallel? [T/F] (must specify --GETdist True)",required=False,default=False,type=distutils.util.strtobool)
    parser.add_argument("-pp","--Plotprops",help="Plot biophysical properties for clusters? [T/F]",required=False,default=False,type=distutils.util.strtobool)
    parser.add_argument("-bt","--boots",help="How many bootstrap replicase should run? [int]",required=False,default=1000,type=int)
    parser.add_argument("-mi","--MIboots",help="How many MI bootstrap replicas should run? [int]",required=False,default=10,type=int)
    parser.add_argument("-aa","--AAorder",help="Order of amino acids for figures [single string of 20]",required=False,default='',type=str)
    parser.add_argument("-cc","--colors",help="Color selection for scatter/line plots [list of strings]",required = False,default=['purple','orange'],type=str,nargs='*')
    ##################################################
    # # # # # # # # # # # # 
    parser.add_argument("-sp","--showProj",help="Show 2D or 3D data projections [string]",required=False,default='both',type=str)
    parser.add_argument("-sc","--showClust",help="Show metadata or clustering [string]",required=False,default='both',type=str)
    parser.add_argument("-nb","--normBar",help="Normalize cluster purity bar plots? [T/F]",required=False,default=True,type=distutils.util.strtobool, nargs=1)
    parser.add_argument("-as","--analysisSel",help="Select the data subset type to further analyze [cluster or metadata]",required=False,default='cluster',type=str)
    parser.add_argument("-lo","--seqlogo",help = "Want to create seqlogo plots? [T/F]",required=False,default=False,type=distutils.util.strtobool)
    parser.add_argument("-ln","--logoNum",help = "What sequence length do you want to generate seqlogos for? [int]",required=False,default=14,type=int)
    parser.add_argument("-sv","--saveSeqs",help = "Want to save clustered sequences in separate files? [T/F]",required=False,default=False,type=distutils.util.strtobool)
    parser.add_argument("-sd","--selDat",help="Specify which precisely which data subsets you want to further analyze [list of ints]",
                        required=False, type=int, nargs='*', default=[0,1])
    parser.add_argument("-p1","--prop1",help="Biophysical property of choice to analyze [int]", required=False,default=1,type=int)
    parser.add_argument("-p2","--prop2",help="Biophysical property of choice to analyze [int]", required=False,default=2,type=int)
    parser.add_argument("-ms","--matSize",help="Matrix size for linear discriminant analysis [int]", required=False,default=10,type=int)
    # Debugging functions:
    parser.add_argument("-sl","--showLabel",help="Show on projections? [T/F]",required=False,default=False,type=distutils.util.strtobool)

    args = parser.parse_args()
    return(args)
###############################################################################################

# All of the 
def run():
    #try:
    args = main()
    #except:
    #    print("Error, check yourself (more error handling to come)")

    # A few of the argparse options don't have default entries
    # because we are dependent on them being 
    ###############################################################################################
    #### Definition of ALL variables used in this notebook: Need to let users add them in with flags
    # For now the defaults are commented out so we can compare/contrast any issues
    datDir = args.datDir #"app_data"
    outputDir = args.outputDir #"AIMS_out"

    # Look to the notebooks for the other examples [ig vs msa vs peptide]
    molecule = args.molecule #'Ig'
    fileName = args.fileNames  #["siv_tl8.csv","siv_cm9.csv"]
    datName = args.datNames #["TL8","CM9"]
    if len(datName) != len(fileName):
        datName = []
        for i in np.arange(len(fileName)):
            datName.append("dat"+str(i))
    num_loop = args.numLoop #1
    drop_duplicates = args.dropDup  #True
    parallel_process = args.parallel  #False

    # Note, you can't really take a subset of peptide data
    subset = args.subset #True
    subset_starts = args.subStart #[164,214,275,327]
    subset_ends = args.subEnd #[214,275,327,376]

    # The "align" option only changes alignments for Ig, not MSA encodings
    align = args.align #'bulge'
    # This "pad" value only changes things if your alignemnt is "bulge"
    # Otherwise keep it defined as 8, again another dummy variable that you need to keep
    pad = args.bulgePad #6

    # Specifically normalization for the bphys property matrix
    normalize = args.normProp # True
    renormalize = args.REnorm # True
    # For the clustering
    dchoice = args.clustData #'parse'
    reduce = args.projAlg #'pca'
    # set an optional umap_seed
    umap_seed = args.umapSeed #69
    # Clustering specifics
    clust = args.clustAlg #'optics'
    # This will be an all in one variable for whichever algorithm
    # minimum cluster size for optics, EPS for DBSCAN, and NClust for kmeans
    clust_size = args.clustSize #10
    # Incorporating metadata
    meta_form = args.metaForm # 'category'

    ############################################
    # New features! Woohoo
    DOstats = args.DOstats # False
    bootstrap = args.DOboot # False
    boots = args.boots # 1000
    MIboots = args.MIboots # 10
    AAorder = args.AAorder # 'WFMLIVPYHAGSTDECNQRK'
    if len(AAorder) == 20:
        custom_key = True
    elif len(AAorder) == 0:
        custom_key = False
    else:
        print("ERROR: Custom Key Wrong Length!")
        # Break out early with this error
        quit()
    colors = args.colors
    GETdist= args.GETdist
    parallel_dist = args.PARdist
    plot_props = args.Plotprops
    ############################################

    ######### Do you want to show 2D projection, 3D, or both?##################
    proj_show = args.showProj # 'both'
    clust_show = args.showClust # 'both'
    # This is more of a debug function
    show_labels = args.showLabel #False
    meta_name = args.metaName # 'Dset'
    # Do you want to noramlize cluster_purity barplots
    norm = args.normBar #True
    # Decide if you want to analyze the data by cluster or by metadata.
    subset_sel = args.analysisSel #'cluster' # other option is 'metadata'
    # Do you want to save each individual cluster of sequences as an AIMS-compatible file?
    seqlogo = args.seqlogo #False
    save_subSeqs = args.saveSeqs #False
    seqlogo_size = args.logoNum #14
    # Which data subsets do you want to use
    sub_sels = args.selDat #[0,1]
    # Biophysical properties shown for positional average
    prop1 = args.prop1 # 1
    prop2 = args.prop2 # 2
    # Matsize for the linear discriminant analysis
    matSize= args.matSize #10

    # THEN, we need a whole bunch of functions so users don't have to sit through this forever.
    # Specifically, ways to only run portions of the analysis repeatedly or to save progress
    # Big one will be an option to save the "bigass matrix" and then load it back in to skip that step
    ###############################################################################################

    # # Section 1: Analysis Mode Definition and Pre-Processing
    # This section is critical for loading in data and either utilizing the AIMS Ig (TCR and Antibody), AIMS Peptide (Eluted from MHC or otherwise), or AIMS MSA (all other molecules) analysis modules

    # Trying to make it as easy as possible for new users to get a grasp of how AIMS works, while also reducing some clutter on the GitHub page. If you preferred older versions of AIMS, where each type of analysis was a completely different Jupyter notebook, you can click "tags" on the GitHub and download version 0.7.5 or earlier. 
    # Note that this has important implications for the way that your data is read in. Only three options, "Ig","MSA", or "Peptide". MSA should work for every type of molecule, even Ig molecules if you really wanted to. Would be a little messy though due to poor conservation of CDR3 loops
    # Even if you only have a single file to analyze, it is required to define your path in a list (should be clearer below in the code)
    # The below cell is where we finally require users to input their own information for their files. How much info depends on the molecules of interest

    # In[2]:
    #####################
    """
    Throughout this script, I will try to include comments like this where explanations are important
    Use the below information to alter this code for your own personal analysis
    - molecule is either "Ig","MSA", or "peptide"
    - datDir is the location of your data
    - outputDir is where figures will be saved
        -If you need help understanding how to define these "Dirs" go to the ReadTheDocs website
    - fileName is a list of the full dataset filenames
    - datName is a more human readable title for the data
    """
    #####################

    try: 
        os.listdir(outputDir)
    except:
        os.makedirs(outputDir)

    # # In the Below Cell We Convert Our Sequences to an AIMS-Readable Format
    # If there are downstream issues in the first section where figures are generated (Section 2) come back here and check formatting of the outputs here
    # # In this Section We Also Determine if We Want Subsets of the Data to be Taken
    # For instance, you may be interested in only a select region of the MSA or only a subset of the included CDR loops
    # In the example subset we provide, we are selecting out regions of the MHC MSA that correspond to consensus alignments of the MHC alpha-helices and beta-strands
    # In[3]:
    if len(fileName) != len(datName):
        print("A mistake! You don't have the proper number of labels for your files")
    elif molecule.lower() == 'ig':
        for i in np.arange(len(fileName)):
            seq_pre = aimsLoad.Ig_loader(datDir+'/'+fileName[i],label=datName[i],loops=num_loop,drop_degens = drop_duplicates)
            if i == 0:
                seqPRE = seq_pre
            else:
                seqPRE = pandas.concat([seqPRE,seq_pre],axis=1)
        seqF = seqPRE
    elif molecule.lower() == 'peptide':
        for i in np.arange(len(fileName)):
            seq_pre = aimsLoad.pep_loader(datDir+'/'+fileName[i],label=datName[i])
            if i == 0:
                seqPRE = seq_pre
            else:
                seqPRE = pandas.concat([seqPRE,seq_pre],axis=1)
        seqF = seqPRE
    elif molecule.lower() == 'msa':
        for i in np.arange(len(fileName)):
            seq_pre = aimsLoad.msa_loader(datDir+'/'+fileName[i],label=datName[i],drop_dups = drop_duplicates)
            if i == 0:
                seqAll = seq_pre
            else:
                seqAll = pandas.concat([seqAll,seq_pre],axis=1)
        # Have to reshape our sequences
        seqs = np.array(seqAll.loc[0].values).reshape(1,len(seqAll.loc[0].values))
        seqPRE = pandas.DataFrame(seqs)
        seqPRE.columns = seqAll.columns
        if subset:
            seqF = aims.get_msa_sub(seqPRE,subset_starts,subset_ends)
        else:
            seqF = seqPRE
        # Save our FASTA headers as metadata. May be useful downstream or not
        metaF = seqAll.loc[1]

    # In[4]:
    mat_size = aims.get_sequence_dimension(np.array(seqF))[0]
    # General changes that need to be done for every type of molecule
    AA_num_key = aims.get_props()[1]
    if num_loop != 1:
        for i in np.arange(len(mat_size)):
            if i == 0:
                xtick_loc = [mat_size[i]/2]
            else:
                pre_loc = sum(mat_size[:i])
                xtick_loc = xtick_loc + [mat_size[i]/2 + pre_loc]
        else:
            xtick_loc = mat_size/2


    # # Section 2: Sequence Visualization via AIMS Matrix Encoding
    # If looking at Ig molecules, you can decide if you would like to align to the center, left, or right of each sequence.
    # There is a 3rd option, "bulge" which aligns the germline regions of CDR3 (and other loops) and then center aligns what is left. 
    # Change the 'align' variable to one of these four. Pretty easy to visualize each time you do so in below matrix

    # New little added bit to create a custom amino acid order
    if custom_key:
        my_AA_key = [a for a in AAorder]
    else:
        # My AA key here is the "standard" AIMS key that has been used in previous papers
        # Note changing the key doesn't change anything BUT the ordering of Amino acids in some figures
        my_AA_key=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

    # Just in case you want to do MSA analysis, you need to use "my_AA_key_dash"
    my_AA_key_dash = my_AA_key + ['-']

    # In[5]:
    if molecule.lower() == 'ig':
        seq_MIpre = aims.gen_tcr_matrix(np.array(seqF),AA_key = my_AA_key,key = AA_num_key, giveSize = mat_size, alignment = align, bulge_pad=pad)
    elif molecule.lower() == 'peptide':
        seq_MIpre = aims.gen_tcr_matrix(np.array(seqF),AA_key = my_AA_key,key = AA_num_key, giveSize = mat_size, alignment = align, bulge_pad=pad)
    elif molecule.lower() == 'msa':
        AA_num_key_dash = np.hstack((AA_num_key,[0]))
        seq_MIpre = aims.gen_MSA_matrix(np.array(seqF),AA_key_dash = my_AA_key_dash,key = AA_num_key_dash, giveSize = mat_size)

    seq_MIf = pandas.DataFrame(np.transpose(seq_MIpre),columns = seqF.columns)
    fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(16,8))
    x = ax[0,0].imshow(np.transpose(seq_MIf), interpolation='nearest', aspect='auto',cmap=cmap)
    pl.colorbar(x)
    ax[0,0].set_ylabel('Sequence Number')
    ######
    # It will help to have vertical dashed black lines to guide the viewer
    seq1_len = np.shape(seqF)[1]
    Numclones = int(seq1_len)
    if type(mat_size) != int:
        for i in np.arange(len(mat_size)-1):
            ax[0,0].plot( (mat_size[i] + sum(mat_size[:i]) - 0.5) * np.ones(Numclones),np.arange(Numclones),'k--',linewidth = 3)
    #######
    ax[0,0].set_xlabel('Sequence Position')
    pl.savefig(outputDir+'/AIMS_mat.pdf',format='pdf')
    pl.close()
    # # Section 3: Calculate our Biophysical Property Matrices
    # Depending on what type of analysis you are doing, this will likely be the slowest step in this entire notebook
    # The only thing that you *might* need to change in the below code is the decision to normalize the biophysical properties or not. Default is to normalize

    # In[6]:
    # Process this new matrix and apply biophysical propery "masks"
    # This has to be changed from the binary case, because we aren't looking for differences
    dsetF = seqF.values

    if molecule.lower() == 'ig':
        special = ''
    elif molecule.lower() == 'peptide':
        special = ''
    elif molecule.lower() == 'msa':
        special = 'MSA'

    #################### PARALLEL PROCESSING TO CREATE BIG MATRIX #######################
    if parallel_process:
        import multiprocessing as mp
        def boot_it(data):
            bigass = classy.get_bigass_matrix(dsetF[:,data[0]:data[1]],AA_key=my_AA_key,AA_key_dash=my_AA_key_dash, giveSize = mat_size, alignment = align,special = special, norm=normalize,bulge_pad=pad)
            return(bigass)
        def do_boot(data):
            with mp.Pool() as pool:
                results = pool.map(boot_it, data)
                return(results)
        if __name__ == "__main__":
            # Probably a smarter way to calculate #seqs per node, but do 100 for now
            final = aims.gen_splits(splitMat = seq_MIf, splitSize = 100)
            big_pre = do_boot(final)
        total_mat = np.concatenate(big_pre, axis = 0)
    else:
        #################### Or, Don't Parallelize to CREATE BIG MATRIX #######################
        bigass = classy.get_bigass_matrix(dsetF,AA_key=my_AA_key,AA_key_dash=my_AA_key_dash, giveSize = mat_size, alignment = align, norm = normalize,special=special,bulge_pad=pad )
        total_mat = bigass
        
    # Generate a large list of property names and matrix positions so you can pinpoint strong
    # contributors to discrimating features between datasets or clusters
    prop_list_old = ['Phobic1','Charge','Phobic2','Bulk','Flex','Kid1','Kid2','Kid3','Kid4','Kid5','Kid6','Kid7','Kid8','Kid9','Kid10']
    prop_list_new = ['Hot'+str(b+1) for b in range(46)]
    prop_names = prop_list_old + prop_list_new
    num_locs = int(np.shape(total_mat)[1]/61)
    Bigass_names = []
    for i in prop_names:
        for j in np.arange(num_locs):
            Bigass_names = Bigass_names + [ i + '-' + str(j) ]

    ########################################################################################
    if renormalize:
        entropy_pre,freq_pre,cov_pre = aims.calculate_shannon(np.transpose(seq_MIf.values))

        repeat_ent = []
        for i in np.arange(len(prop_names)):
            repeat_ent = repeat_ent + [2**entropy_pre]

        refactor = np.array(repeat_ent).reshape(61*len(entropy_pre))

        pp_mat = total_mat*refactor
    else:
        pp_mat = total_mat
    ##########################################################################################

    # Drop Highly Correlated Vectors and Vectors where entry=0 for all entries
    ###### Currently drop vectors with over 0.75 corr. coef. ################
    full_big = pandas.DataFrame(pp_mat,columns = Bigass_names)
    drop_zeros = [column for column in full_big.columns if all(full_big[column] == 0 )]
    y = full_big.drop(full_big[drop_zeros], axis=1)
    z_pre = np.abs(np.corrcoef(np.transpose(y)))
    z = pandas.DataFrame(z_pre,columns=y.columns,index=y.columns)
    # Select upper triangle of correlation matrix
    upper = z.where(np.triu(np.ones(z.shape), k=1).astype(bool))
    # If you did want to change that corr. coef. cutoff, do so here
    to_drop = [column for column in upper.columns if ( any(upper[column] > 0.75) ) ]

    # Your final product of a parsed matrix
    parsed_mat = y.drop(y[to_drop], axis=1)

    # This is a new, important variable to account for the cases where renormalization
    # is used. We need non-renormed data for downstream repertoire characterization
    NonNorm_big = pandas.DataFrame(total_mat,columns = Bigass_names)

    # Let's have some default metadata we can pull from later
    tokenized_dset = []
    for i in np.arange(len(datName)):
        for j in seqF.columns:
            if str(j).find(datName[i]) != -1:
                tokenized_dset.append(i)
    token_df = pandas.DataFrame(tokenized_dset,columns=['ID'])
    IDed_full_big = pandas.concat([full_big,token_df],axis=1)

    # Lastly, create a good-ole traditional averaged bphys property matrix
    # i.e. each sequence gets a single value for averaged charge, averaged flexibility, etc...
    posLen,seqLen = np.shape(seq_MIf)
    # The 61 is hardcoded here because it is our number of properties. Eventually we will let users define which properties to use
    seq_bigReshape = np.array(full_big).reshape(seqLen,61,posLen)


    # # Section 4: Sequence Projection & Clustering
    # # The next few cells are particularly powerful for isolating interesting populations in the dataset using PCA, UMAP, and KMeans Clustering
    # This is really a section where you should take your time and toggle some of these settings. Look at your data using PCA and UMAP
    # Try to use "full", "parse", or "avg" biophysical properties for each entry as input into the dimensionality reduction 
    # The below cell is important to define, so we set it apart from the rest of the code

    # In[7]:
    # DEFINE WHICH FORM OF THE DATASET YOU WOULD LIKE TO ANALYZE HERE
    # Perform dim. red. on the whole dataset? enter "full"
    # Want to do it on a dataset with highly correlated vectors removed? enter "parse"
    # Lastly, can just do it on a matrix of the per-sequence average over 61 props? enter "avg"
    # # The below section then uses the above information to actually calculate these things
    # 
    # There is a lot of opportunity to go EVEN deeper into the code here to tweak your analysis. Important to do if you really care about the data. Alter things like:
    # - Nclust for Kmeans
    # - min_samples for OPTICS
    # - eps for DBSCAN
    # - Setting random seeds for UMAP (see the ReadTheDocs for the disclaimer for using UMAP)

    # In[8]:
    # Don't change these if statements, change the stuff below that
    if dchoice.lower() == 'full':
        chosen_dset = full_big
    elif dchoice.lower() == 'parse':
        chosen_dset = parsed_mat
    elif dchoice.lower() == 'avg':
        chosen_dset = np.average(seq_bigReshape,axis=2)

    if reduce == 'pca':
        from sklearn.decomposition import PCA
        pca = PCA(n_components=3, svd_solver='full')
        final=pca.fit_transform(chosen_dset)
        transform = pandas.DataFrame(np.transpose(final),columns = seq_MIf.columns)
        print("PCA Explained Variance Ratio:")
        print(pca.explained_variance_ratio_)
    elif reduce == 'umap':
        import umap
        if isinstance(umap_seed,int):
            reducer = umap.UMAP(n_components=3, n_neighbors = 25,random_state=umap_seed) 
        else:
            reducer = umap.UMAP(n_components=3, n_neighbors = 25)
        final = reducer.fit_transform(chosen_dset)
        transform = pandas.DataFrame(np.transpose(final),columns = seq_MIf.columns)

    # Cluster the results:
    import sklearn.cluster as cluster
    clust_input = np.array(np.transpose(transform))
    if clust == 'kmean':
        NClust = clust_size
        clusts = cluster.KMeans(n_clusters=NClust).fit_predict(clust_input)
    elif clust == 'optics':
        clusts = cluster.OPTICS(min_samples=clust_size).fit_predict(clust_input)
    elif clust == 'dbscan':
        clusts = cluster.DBSCAN(eps=clust_size).fit_predict(clust_input)

    cluster_dset = pandas.DataFrame(clusts,columns=['cluster'])
    # # SECTION 5: Metadata Incorporation
    # By default, we'll assume that the user DOES NOT have any metadata to add (and the code reflects this), but there are sections to show how metadata could be incorporated are included in the code. The ReadTheDocs will be updated soon to go into more detail here for you to add custom metadata

    # In[9]:
    # We can incorporate metadata either defining a categorical map or a quantitative map
    ###################################################################################
    # meta_form is either "category" or "quant". Check ReadTheDocs if more descriptions are needed
    # If using default metadata (i.e. loaded file), it should be "category"
    ###################################################################################

    if meta_form == 'category':
        # This is the easiest default metadata definition.
        # Just based on the files that were loaded in (useless if only 1 file)
        metapre = token_df
        # Convert the metadata from numbers to strings
        meta_conv = []
        for i in token_df.values:
            meta_conv.append(datName[i[0]])

        metadat = pandas.DataFrame(meta_conv)
        # However, you can also read in your own metadata, or build it from scratch
        # JUST MAKE SURE YOUR METADATA IS A PANDAS DATAFRAME IN THE END
    elif meta_form == 'quant':
        # This "quantitative" metadata will just count from 1 to the length of the dataset,
        # coloring the points on the plot in order.
        metaPRE = np.arange(len(token_df))
        metadat = pandas.DataFrame(metaPRE)
        # But, you could add things like MFI, binding affinity, or GEX data for a particular gene

    # From there, not much should change unless you want to give a unique name to your metadata
    meta_map = aims.encode_meta(metadat)
    meta_map.columns = [meta_name]
    meta_leg = metadat.drop_duplicates().values

    # We also need to define our clusters more clearly
    clust_map = cluster_dset
    clust_leg = [a[0] for a in cluster_dset.drop_duplicates().sort_values('cluster').values]
    clust_name = 'cluster'


    # # Section 6: Plotting Clustering and (Optional) Metadata
    # The previous section was just for calculating these things. We are now going to provide users with a few different options to visualize their data
    # 
    # By default we show 2D and 3D cluster projections with metadata and clustered data but give users control over which is shown

    # In[10]:
    # Could optionally plot other data or change legends if you would like
    chosen_map1 = clust_map; leg1 = clust_leg
    chosen_map2 = meta_map; leg2 = meta_leg

    # Define a colormap to color metadata. Can change this if you want, look at matplotlib colormap options
    cmapF = pl.get_cmap('tab20b')
    # So there's no reason to have a function for this, other than to make this notebook look a bit prettier
    # Really just defining a bunch of stuff repeatedly with if statements
    fig3d,plotloc,plottype,plotem,legends,dattype = aims.get_plotdefs(clust_show,proj_show,chosen_map1,chosen_map2,leg1,leg2)

    # Now plot the actual stuff, and save the object handles in a list
    colorhandle= []; ax=[]
    for i in np.arange(len(plotloc)):
        print(dattype[i])
        if dattype[i] == 'clust':
            if clust == 'kmean':
                cmap_use = pl.get_cmap('rainbow')
            else:
                cmap_use = cmap
        else:
            cmap_use = cmapF
        if plottype[i] == '3d':
            ax.append(fig3d.add_subplot(plotloc[i],projection=plottype[i]))
            # So the 3D plot has shading, which means instead we're better off splitting the data up 
            colorhandle.append(ax[i].scatter(clust_input[:,0],clust_input[:,1],clust_input[:,2],c = plotem[i].values.reshape(len(plotem[i]),), cmap=cmap_use))
            ax[i].set_xlabel('AX1',labelpad=20); ax[i].set_ylabel('AX2',labelpad=20); ax[i].set_zlabel('AX3',labelpad=15)
        else:
            ax.append(fig3d.add_subplot(plotloc[i]))
            colorhandle.append(ax[i].scatter(clust_input[:,0],clust_input[:,1],c = plotem[i].values.reshape(len(plotem[i]),), cmap=cmap_use))
            ax[i].set_xlabel('AX1'); ax[i].set_ylabel('AX1')

    # This code is somewhat wild just for getting a properly color-coded legend in there...
    fig3d.canvas.draw()
    # We use the metadata mapped colors down the line, so if you aren't plotting it then you need to save some other way
    need_meta = True
    for i in np.arange(len(legends)):
        cmap_pre = pandas.DataFrame(colorhandle[i].get_facecolors())
        if dattype[i] == 'clust':
            mapDF = pandas.concat([clust_map,cmap_pre],axis=1)
            mapped_colors= mapDF.sort_values('cluster').drop_duplicates('cluster').values[:,1:]
        else:
            need_meta = False
            mapDF = pandas.concat([meta_map,cmap_pre],axis=1)
            mapped_colors= mapDF.sort_values(meta_name).drop_duplicates(meta_name).values[:,1:]
            keep_map = mapped_colors
        # Do two things at once here. Also add in axis labels
        # Don't plot duplicate legends
        if i > 0:
            if dattype[i] == dattype[i-1]:
                continue
        if len(legends[i]) > 6:
            # Don't show exceedingly long legends
            continue
        else:
            legend_elements = []
            # This line won't do anything if the array is properly shaped, will do stuff it if isn't.
            legends[i] = np.array(legends[i]).reshape(len(legends[i]))
            for j in np.arange(len(legends[i])):
                element = [Line2D([0], [0], marker='o', color='w', label=legends[i][j],markerfacecolor=mapped_colors[j], markersize=10)]
                legend_elements = legend_elements+element
            ax[i].legend(bbox_to_anchor=(0.5, 1.1),handles=legend_elements,loc='upper center',ncol=len(legends[i]),title=meta_name)

    if need_meta:
        cmap_discrete = cmapF(np.linspace(0, 1, len(clust_input)))
        cmap_pre = pandas.DataFrame(cmap_discrete)
        mapDF = pandas.concat([meta_map,cmap_pre],axis=1)
        mapped_colors= mapDF.sort_values(meta_name).drop_duplicates(meta_name).values[:,1:]
        keep_map = mapped_colors


    if show_labels:
        # ONLY show this for a 2D plot
        for num in np.arange(len(ax)):
            if plottype[num] == '2d':
                break
        a = 0; plot1 = clust_input[:,0]; plot2 = clust_input[:,1]
        plot_labels = seqF.columns.values
        for i,j in zip(plot1,plot2):
            ax[num].annotate(str(plot_labels[a]),xy=(i,j),fontsize=14)
            a+=1

    pl.savefig(outputDir+'/AIMS_projections.png',format='png')
    pl.close()

    # # Section 7: Quantify Cluster Compositions
    # This step might be meaningless if you don't have ANY metadata to go off of, and/or are not comparing/contrasting two datasets. This is fine, you should still run this step to make sure that everything is properly defined
    # 

    # In[11]:
    cmap3 = pl.get_cmap('tab20b')
    fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(14,8))
    # Just in case your metadata has weird indices:
    meta_map.index = cluster_dset.index
    final_breakdown = pandas.concat([cluster_dset,meta_map],axis=1)
    a = 0
    for i in np.sort(final_breakdown['cluster'].drop_duplicates()):
        if i == -1:
            continue # Dont count the unclustered
        sub_clust = final_breakdown[final_breakdown['cluster'] == i]
        if len(sub_clust) == 0:
            continue
        bottom=0

        for j in sub_clust[meta_name].drop_duplicates().values:
            sub_sub = sub_clust[sub_clust[meta_name] == j]
            if norm:
                pl.bar(a,len(sub_sub)/len(sub_clust),bottom = bottom,color=keep_map[int(j)],edgecolor='black')
                bottom += len(sub_sub)/len(sub_clust)
            else:
                pl.bar(a,len(sub_sub),bottom = bottom,color=keep_map[int(j)],edgecolor='black')
                bottom += len(sub_sub)
        a = a+1

    if len(meta_leg) < 6:
        meta_legF = np.array(meta_leg).reshape(len(meta_leg))
        legend_elements= []
        for j in np.arange(len(meta_leg)):
            element = [Line2D([0], [0], marker='o', color='w', label=meta_legF[j],markerfacecolor=keep_map[j], markersize=10)]
            legend_elements = legend_elements+element
        pl.legend(handles=legend_elements,ncol=len(meta_leg))
    pl.xlabel('Cluster Number')
    if norm:
        pl.ylabel('Fraction Per Cluster')
    else:
        pl.ylabel('Count Per Cluster')

    cluster_purity = aims.calc_cluster_purity(final_breakdown,meta_name)
    print("Cluster Purities:")
    print(np.transpose(cluster_purity))

    # Really not sure why *just* this file is mad being saved as a pdf. 
    # Had to put this line in to try to fix it though.

    try:
        pl.savefig(outputDir+'/AIMS_clusterQuant.pdf',format='pdf')
    except:
        pl.savefig(outputDir+'/AIMS_clusterQuant.png',format='png')
    pl.close()

    ###### CALCULATE CLUSTER SIGNIFICANCE ###########
    fig, ax = pl.subplots(1, 2,squeeze=False,figsize=(16,6))

    cluster_purity,purity_pVal = aims.calc_cluster_purity(final_breakdown,meta_name)

    x = ax[0,0].imshow(cluster_purity,interpolation='nearest', aspect='auto',vmax = 1,cmap='plasma')
    y = ax[0,1].imshow(purity_pVal,interpolation='nearest', aspect='auto',vmax = 0.05,cmap='Greys_r')

    ax[0,0].set_xlabel('MetaData #'); ax[0,1].set_xlabel('MetaData #')
    ax[0,0].set_ylabel('Cluster #'); ax[0,1].set_ylabel('Cluster #')
    ax[0,0].set_title('Cluster Purity'); ax[0,1].set_title('Cluster Significance')
    pl.colorbar(x,ax=ax[0,0]); pl.colorbar(y,ax=ax[0,1])

    pl.savefig(outputDir+'/AIMS_clusterPurity.pdf',format='pdf')
    pl.close()

    # # Section 8: Defining Data Subsets of Interest to Further Characterize Repertoire
    # VERY IMPORTANT STEP here for all downstream analysis. Define whether you want to compare/contrast clustered sequences or analyze sequences based upon the associated metadata

    # In[12]:
    # # Then, Plot This Dataset of Interest and Visualize How Similar/Different Each Look

    # In[13]:

    # NEW FEATURE!!!
    # You can view biophysical properties as
    # a function of cluster or metadata here,
    # rather than just the sequence colors
    #########################################
    # Recommend changing colors if showing diff props
    showProp = 1 #1=Charge, 2=phob
    # Show lines separating major groups?
    show_lines = True
    #########################################

    if subset_sel.lower() == 'cluster':
        chosen_map = clust_map; chosen_name = clust_name
    elif subset_sel.lower() == 'metadata':
        chosen_map = meta_map; chosen_name = meta_name
    fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(16,8))
    for i in np.sort(chosen_map[chosen_name].drop_duplicates()):
        if i == -1:
            continue
        if plot_props:
            sub_props = seq_bigReshape[:,showProp]
            subbDF = pandas.DataFrame(np.transpose(sub_props))
            subbDF.columns = seq_MIf.columns
            pre_clust = subbDF[subbDF.columns[chosen_map[chosen_map[chosen_name] == i].index]]
        else:
            pre_clust = seq_MIf[seq_MIf.columns[chosen_map[chosen_map[chosen_name] == i].index]]
        clustID = np.transpose(pandas.DataFrame(i*np.ones(np.shape(pre_clust)[1])))
        clustID.columns = pre_clust.columns
        pre_clustF = pandas.concat([pre_clust,clustID],axis=0)
        if i == 0:
            clustered = pre_clustF
        else:
            clustered = pandas.concat([clustered, pre_clustF],axis = 1)
        if show_lines:
            ax[0,0].plot(np.arange(len(seq_MIf)),np.ones(len(seq_MIf))*(np.shape(clustered)[1]),'black',linewidth = 3)
    ax[0,0].set_ylabel('Sequence Number')
    if type(mat_size) != int:
        for i in np.arange(len(mat_size)-1):
            ax[0,0].plot( (mat_size[i] + sum(mat_size[:i]) - 0.5) * np.ones(np.shape(clustered)[1]),np.arange(np.shape(clustered)[1]),'k--',linewidth = 3)
    if plot_props:
        ttt = np.transpose(np.array(clustered))[:,:-1]
        scaledd = np.max([np.abs(np.min(ttt)),np.abs(np.max(ttt))])
        xyz = ax[0,0].imshow(np.transpose(np.array(clustered))[:,:-1], interpolation='nearest', aspect='auto',cmap='bwr',vmin=-scaledd,vmax=scaledd)
    else:
        xyz = ax[0,0].imshow(np.transpose(np.array(clustered))[:,:-1], interpolation='nearest', aspect='auto',cmap=cmap)
    pl.colorbar(xyz)
    pl.savefig(outputDir+'/AIMS_'+subset_sel+'_subViz.pdf',format='pdf')
    pl.close()


    ###############################################################
    # NEW: AIMSDIST!
    if GETdist:
        get_distClusts = True

        for i in chosen_map.sort_values(chosen_name).drop_duplicates().values:
            if i == -1:
                continue
            sub_MI_temp = seq_MIf[seq_MIf.columns[chosen_map[chosen_map[chosen_name] == i[0]].index]]
            sub_seqs_temp = np.transpose(seqF[sub_MI_temp.columns])
            if i == 0:
                sorted_seqs = sub_seqs_temp
            else:
                sorted_seqs = pandas.concat([sorted_seqs,sub_seqs_temp])
        ########################################################################################################3
        if parallel_dist:
            import multiprocessing as mp
            def boot_it(data):
                if data[0][0] == data[1][0]:
                    dist_temp = aims.calc_AIMSdist(sorted_seqs[data[0][0]:data[0][1]])
                else:
                    dist_temp = aims.calc_AIMSdist(sorted_seqs[data[0][0]:data[0][1]],sorted_seqs[data[1][0]:data[1][1]])
                return(data,dist_temp)
            def do_boot(data):
                with mp.Pool() as pool:
                    results = pool.map(boot_it, data)
                    return(results)
            if __name__ == "__main__":
                # Probably a smarter way to calculate #seqs per node, but do 100 for now
                xx = aims.prep_distCalc(sorted_seqs)
                dist_pre = do_boot(xx)
            dist_matF = np.zeros((len(sorted_seqs),len(sorted_seqs)))
            for i in np.arange(len(dist_pre)):
                # Set 1 will be our x-axis of the matrix
                # set 2 will be our y-axis of the matrix
                set1 = dist_pre[i][0][0]
                set2 = dist_pre[i][0][1]
                # Set 3 is then the data that goes in that space
                set3 = dist_pre[i][1]

                # Need to fill both the matrix entry and the 
                # transpose of that entry!!!
                dist_matF[set1[0]:set1[1],set2[0]:set2[1]] = set3
                dist_matF[set2[0]:set2[1],set1[0]:set1[1]] = np.transpose(set3)
            dists = dist_matF
        else:
            dists = aims.calc_AIMSdist(sorted_seqs)
        fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(10,8))
        x = pl.imshow(np.transpose(dists), interpolation='nearest', aspect='auto',vmin=0,vmax=11)

        # Optionally can get back distance clusters:
        if get_distClusts:
            distance_clusters = aims.get_distClusts(dists,metadat,max_d=5)
            distance_clusters.to_csv(outputDir+'/dist_clust.csv',index=False)
        pl.colorbar(x)
        pl.savefig(outputDir+'/AIMSdist.pdf',format='pdf')
        pl.close

    # # Section 9: Isolation of Individual Groups for Downstream Characterization
    # Run the next cell to visualize your clusters of choice, and then go through the remainder of the AIMS modules
    # 
    # Selecting the sub_sels [0,1] will select either the first two clusters or the first two metadata entries. Remember that python is 0-indexed, so if you want to look at a very specific metadata cluster then make sure you take that into account! You can *technically* visualize every cluster all at once, but that is REALLY not recommended

    # In[14]:
    # So now that we've done some clustering, pick out the most interesting or sections of the data:
    # OPTIONALLY YOU CAN VISUALIZE/ANALYZE EVERY CLUSTER. STRONGLY NOT RECOMMENDED IF YOU HAVE MANY CLUSTERS
    ##sub_sels = np.arange(len(chosen_map))

    # Redefine the raw sequences so we can manipulate them if needed
    dset = seqF
    fig, ax = pl.subplots(len(sub_sels), 1,squeeze=False,figsize=(16,4*len(sub_sels)))
    label=[]
    a = 0
    for i in sub_sels:
        # Look at umap dset or pca dset
        sub_MI = seq_MIf[seq_MIf.columns[chosen_map[chosen_map[chosen_name] == i].index]]
        sub_seqs = np.transpose(dset[sub_MI.columns])
        if subset_sel.lower() == 'metadata':
            label.append(meta_legF[i])
        else:
            label.append('cluster'+str(i))
        ax[a,0].imshow(np.transpose(sub_MI), interpolation='nearest', aspect='auto',cmap=cmap)

        datlen = np.shape(sub_MI)[1]
        datID = np.transpose(pandas.DataFrame(datlen*[a]))
        datID.columns = sub_MI.columns
        sub_matPRE = pandas.concat([sub_MI,datID],axis=0)
        if a == 0:
            sub_matF = sub_matPRE
            sub_seqF = sub_seqs
        else:
            sub_matF = pandas.concat([sub_matF,sub_matPRE],axis=1)
            sub_seqF = pandas.concat([sub_seqF,sub_seqs],axis=0)
        a+=1 

        if save_subSeqs:
            sub_seqs.to_csv(outputDir+'/'+label[i]+'_all.txt',header=None,index=None)
        if seqlogo:
            seqlogo1 = sub_seqs[sub_seqs[0].str.len() == seqlogo_size]
            seqlogo1.to_csv(outputDir+'/'+label[i]+'_logo.txt',header=None,index=None)

    print("Visualizing the subsets: "+str([str(a) for a in label]))
    pl.savefig(outputDir+'/AIMS_'+subset_sel+'_selected.pdf',format='pdf')
    pl.close()


    # # Leaving This Here in Case You Have Seqlogo Installed
    # Can be useful to visualize some sequences in this way. This code has not been tested in a while, and has not been tested at all for MSA analysis

    # In[15]:
    # # Section 10: Generate Subsets of Matrices We Have Calculated Previously
    # This will save some time compared to how things were done before.

    # In[16]:


    # define "RE" variables so we don't mess with any variables upstream (in case you want to re-run)
    full_big_re = NonNorm_big; full_big_re.index = seq_MIf.columns
    parsed_mat_re = parsed_mat; parsed_mat_re.index = seq_MIf.columns
    # use the transpose of the sub_mat to find the 
    ref_sub = np.transpose(sub_matF)
    ref_sub.columns = np.arange(len(sub_matF))

    sub_big = full_big_re.loc[sub_matF.columns]
    sub_parsed = parsed_mat_re.loc[sub_matF.columns]


    # # Section 11: Position Sensitive Biophysical Properties for Every Clone in the Dataset
    # Here, we can visualize the how similar or dissimilar the biophysical properties are for each sequence within a given repertoire.
    # 
    # Right now we show only the charge and the hydropathy, but you can change "prop1" in either cell to visualize a different property (see ReadTheDocs for more info)

    # In[17]:


    # Generate the position sensitive charge across all clones in the dataset
    # Which property will you want to look at down the line? 1 = charge, 2 = hydrophobicity... see full list in eLife paper (Boughter et al. 2020)

    fig, axs = pl.subplots(1, len(sub_sels),squeeze=False,figsize=(10*len(sub_sels),8))

    for dat in np.arange(len(sub_sels)):
        take_sub = ref_sub[ref_sub[len(sub_matF)-1] == dat].index
        take_big = sub_big.loc[take_sub]
        temp_bigReshape = np.array(take_big).reshape(len(take_big),61,posLen)
        min_temp = np.min(temp_bigReshape[:,prop1,:])
        max_temp = np.max(temp_bigReshape[:,prop1,:])
        if abs(min_temp) > abs(max_temp):
            propMin = min_temp
            propMax = -min_temp
        elif abs(max_temp) > abs(min_temp):
            propMin = -max_temp
            propMax = max_temp
        else:
            propMin = min_temp
            propMax = max_temp
        x = axs[0,dat].imshow(temp_bigReshape[:,prop1,:],interpolation='nearest', aspect='auto',cmap=cm.seismic, vmin = propMin, vmax = propMax)
        axs[0,dat].set_xlabel('Sequence Position'); axs[0,dat].set_ylabel('Sequence Number')
        axs[0,dat].set_title(str(label[sub_sels[dat]]) + ' - Charge')
        fig.colorbar(x, ax=axs[0,dat])

    pl.savefig(outputDir+'/AIMS_'+subset_sel+'_charge.pdf',format='pdf')
    pl.close()


    # In[18]:


    fig, axs = pl.subplots(1, len(sub_sels),squeeze=False,figsize=(10*len(sub_sels),8))

    for dat in np.arange(len(sub_sels)):
        take_sub = ref_sub[ref_sub[len(sub_matF)-1] == dat].index
        take_big = sub_big.loc[take_sub]
        temp_bigReshape = np.array(take_big).reshape(len(take_big),61,posLen)
        min_temp = np.min(temp_bigReshape[:,prop2,:])
        max_temp = np.max(temp_bigReshape[:,prop2,:])
        if abs(min_temp) > abs(max_temp):
            propMin = min_temp
            propMax = -min_temp
        elif abs(max_temp) > abs(min_temp):
            propMin = -max_temp
            propMax = max_temp
        else:
            propMin = min_temp
            propMax = max_temp
        x = axs[0,dat].imshow(temp_bigReshape[:,prop2,:],interpolation='nearest', aspect='auto',cmap=cm.BrBG, vmin = propMin, vmax = propMax)
        axs[0,dat].set_xlabel('Sequence Position'); axs[0,dat].set_ylabel('Sequence Number')
        axs[0,dat].set_title(str(label[sub_sels[dat]]) + ' - Hydropathy')
        fig.colorbar(x, ax=axs[0,dat])

    pl.savefig(outputDir+'/AIMS_'+subset_sel+'_hydropathy.pdf',format='pdf')
    pl.close()


    # # Section 12: Averaged Biophysical Properties for Each Group of Interest
    # Look at biophysical properties averaged over clones AND position
    # 
    # Note for old users of the software, you might get different looking results because originally I normalized vectors to unit length, but NOT 0 mean. I now do both.

    # In[19]:


    fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(16,8))
    x_axis = np.array([-0.2,0.9,2,3.1])
    # We want to exclude prop0 (the simple 1-21 AA representation entries)
    full_avg = []; full_std = []; plot_lab = []
    a=0
    for dat in np.arange(len(sub_sels)):
        sin_avg =[]; sin_std = []
        b=0
        for prop in np.arange(4):
            propF = prop+1
            take_sub = ref_sub[ref_sub[len(sub_matF)-1] == dat].index
            take_big = sub_big.loc[take_sub]
            temp_bigReshape = np.array(take_big).reshape(len(take_big),61,posLen)
            if bootstrap:
                prop_avg = []
                for i in np.arange(boots):
                    re_big = resample(temp_bigReshape)
                    prop_avg.append(np.average(re_big[:,propF,:]))
                fin_avg = np.average(prop_avg,axis=0)
                fin_std = np.std(prop_avg,axis=0)
                if b == 0:
                    plot_lab.append(pl.bar(x_axis[b]+a/len(sub_sels), fin_avg,yerr=fin_std,width=1/len(sub_sels),alpha=0.5,color=colors[dat]))
                else:
                    pl.bar(x_axis[b]+a/len(sub_sels), fin_avg,yerr=fin_std,width=1/len(sub_sels),alpha=0.5,color=colors[dat])
            
            else:
                # We don't want to plot prop1, we want to plot the rest of them
                plot_avg = np.average(np.average(temp_bigReshape[:,propF,:],axis=1))
                plot_std = np.std(np.std(temp_bigReshape[:,propF,:],axis=1))
                sin_avg.append(plot_avg)
                sin_std.append(plot_std)
            b+=1
        a+=1

        if bootstrap == False:
            full_avg.append(sin_avg)
            full_std.append(sin_std)

    if bootstrap == False:
        for i in np.arange(len(sub_sels)):
            ax[0,0].bar(x_axis+i/len(sub_sels), full_avg[i],yerr = full_std[i],alpha = 0.5, width = 1/len(sub_sels),color=colors[i])
            #ax[0,0].bar(x_axis+i/len(sub_sels), full_avg[i],alpha = 0.5, width = 1/len(sub_sels))
            ax[0,0].legend(label)
    else:
        ax[0,0].legend(plot_lab,label)
    ax[0,0].set_xticks([0.2,1.3,2.4,3.5])
    ax[0,0].set_xticklabels(['Charge','Hydrophobicity','Bulkiness','Flexibility'])
    ax[0,0].set_xlabel('Biophysical Property')
    ax[0,0].set_ylabel('Normalized Property Value')
    pl.savefig(outputDir+'/AIMS_'+subset_sel+'_netAvgProp.pdf',format='pdf')
    pl.close()

    ##########################################################################
    # Calculating statistical significance
    if DOstats:
        prop_names = ['Charge','Hydrophobicity','Bulkiness','Flexibility']
        take_sub1 = ref_sub[ref_sub[len(sub_matF)-1] == 0].index
        take_big1 = sub_big.loc[take_sub1]
        temp_bigReshape1 = np.array(take_big1).reshape(len(take_big1),61,posLen)
        take_sub2 = ref_sub[ref_sub[len(sub_matF)-1] == 1].index
        take_big2 = sub_big.loc[take_sub2]
        temp_bigReshape2 = np.array(take_big2).reshape(len(take_big2),61,posLen)

        # Significance for Bar plots
        bar_sig = []; bar_names = []
        for i in np.arange(4):
            propF = i+1
            data1 = np.average(temp_bigReshape1[:,propF,:],axis=1)
            data2 = np.average(temp_bigReshape2[:,propF,:],axis=1)
            p = aims.do_statistics(data1,data2,num_reps=10000,test='average')
            bar_sig = bar_sig + [p]
            bar_names = bar_names + [prop_names[i]]
            print('p-value for '+prop_names[i]+' bar plot: '+str(p))
        BARsigF = np.vstack((bar_names,bar_sig))
        pandas.DataFrame(BARsigF).to_csv(outputDir+'/bar_sig.csv',index=False)


    # # Section 13: Position-Sensitive Averaged Biophysical Properties
    # Here we effectively average the figures from section 11 over the y-axes, giving a general idea of trends in the biophysical properties of within-group sequences
    # 
    # NOTE: This figure frequently looks chaotic/hard to parse for MSA analysis. More useful for Ig analysis

    # In[20]:


    # Now get the position sensitive avarege biophysical properties
    fig, ax = pl.subplots(2, 1,squeeze=False,figsize=(14,10))
    prop_sel = [1,2]
    full_avg = []; full_std = []

    for dat in np.arange(len(sub_sels)):
        a=0
        for prop in prop_sel:
            take_sub = ref_sub[ref_sub[len(sub_matF)-1] == dat].index
            take_big = sub_big.loc[take_sub]
            temp_bigReshape = np.array(take_big).reshape(len(take_big),61,posLen)
            if bootstrap:
                prop_avg = []
                for i in np.arange(boots):
                    re_big = resample(temp_bigReshape)
                    prop_avg.append(np.average(re_big[:,prop,:],axis=0))
                fin_avg = np.average(prop_avg,axis=0)
                fin_std = np.std(prop_avg,axis=0)
                ax[a,0].plot(fin_avg,marker='o',linewidth=2.5,color=colors[dat])
                ax[a,0].fill_between(np.arange(len(fin_avg)),fin_avg+fin_std,fin_avg-fin_std,alpha=0.3,color=colors[dat])
            else:
                take_sub = ref_sub[ref_sub[len(sub_matF)-1] == dat].index
                take_big = sub_big.loc[take_sub]
                temp_bigReshape = np.array(take_big).reshape(len(take_big),61,posLen)
                plot_avg = np.average(temp_bigReshape[:,prop,:],axis=0)
                plot_std = np.std(temp_bigReshape[:,prop,:],axis=0)

                ax[a,0].plot(plot_avg,marker='o',linewidth=2.5,color=colors[dat])
                ax[a,0].fill_between(np.arange(len(plot_avg)),plot_avg+plot_std,plot_avg-plot_std,alpha=0.3,color=colors[dat])
            a+=1

    # Draw some nice lines to guide 
    y11, y12 = ax[0,0].get_ylim();y21, y22 = ax[1,0].get_ylim()
    if type(mat_size) != int:
        for i in np.arange(len(mat_size)-1):
            ax[0,0].plot( (mat_size[i] + sum(mat_size[:i]) - 0.5) * np.ones(100),np.linspace(y11,y12,100),'black',linewidth = 3)
            ax[1,0].plot( (mat_size[i] + sum(mat_size[:i]) - 0.5) * np.ones(100),np.linspace(y21,y22,100),'black',linewidth = 3)

    legend_elements=[]
    for j in np.arange(len(sub_sels)):
        element = [Line2D([0], [0], marker='o', color='w', label=label[j],markerfacecolor=colors[j], markersize=10)]
        legend_elements = legend_elements+element
    pl.legend(handles=legend_elements,ncol=len(sub_sels))

    ax[0,0].set_ylabel('Normalized Charge')
    ax[1,0].set_ylabel('Normalized Hydropathy')
    if type(mat_size)==int:
        ax[0,0].set_xlim([-0.5,mat_size-0.5])
        ax[1,0].set_xlim([-0.5,mat_size-0.5])
    else:
        ax[0,0].set_xlim([-0.5,sum(mat_size)-0.5])
        ax[1,0].set_xlim([-0.5,sum(mat_size)-0.5])

    pl.xlabel('Sequence Position')
    pl.savefig(outputDir+'/AIMS_posSensAvg.pdf',format='pdf')
    pl.close()

    ##########################################################################
    # Calculating statistical significance
    if DOstats:
        fig, ax = pl.subplots(2, 1,squeeze=False,figsize=(16,10))
        take_sub1 = ref_sub[ref_sub[len(sub_matF)-1] == 0].index
        take_big1 = sub_big.loc[take_sub1]
        temp_bigReshape1 = np.array(take_big1).reshape(len(take_big1),61,posLen)
        take_sub2 = ref_sub[ref_sub[len(sub_matF)-1] == 1].index
        take_big2 = sub_big.loc[take_sub2]
        temp_bigReshape2 = np.array(take_big2).reshape(len(take_big2),61,posLen)

        a=0
        for propF in prop_sel:
            #Significance for position-sensitive
            data1 = temp_bigReshape1[:,propF,:]
            data2 = temp_bigReshape2[:,propF,:]

            p = aims.do_statistics(data1,data2,num_reps=10000,test='average')
            ax[a,0].plot(p,color='black',marker='o')
            ax[a,0].plot(np.arange(np.shape(data1)[1]),np.ones(np.shape(data1)[1])*0.05,linewidth=3,color='red')
            a+=1

        y11, y12 = ax[0,0].get_ylim();y21, y22 = ax[1,0].get_ylim()
        if type(mat_size) != int:
            for i in np.arange(len(mat_size)-1):
                ax[0,0].plot( (mat_size[i] + sum(mat_size[:i]) - 0.5) * np.ones(100),np.linspace(y11,y12,100),'black',linewidth = 3,linestyle='--')
                ax[1,0].plot( (mat_size[i] + sum(mat_size[:i]) - 0.5) * np.ones(100),np.linspace(y21,y22,100),'black',linewidth = 3,linestyle='--')

        pl.xlabel('Sequence Position')
        ax[0,0].set_ylabel('Charge p-value')
        ax[1,0].set_ylabel('Hydropathy p-value')
        pl.savefig(outputDir+'/AIMS_PosSens_pval.pdf',format='pdf')
        pl.close()


    # # Section 14: Information Theoretic Calculations
    # Use the Shannon Entropy and Mutual Information to quantify the diversity and inter-relations between the amino acids used in each sequence within a group
    # NOTE: These metrics are more useful for large datasets, less so far small datasets from DBSCAN or OPTICS identified clustered
    # Also NOTE: Again, these plots can be a bit hard to interpret for longer MSA sequences. Just too many points to read.

    # In[21]:

    # Calculate the Shannon Entropy, a proxy for diversity
    fig = pl.figure(figsize=(16,8))
    gs = gridspec.GridSpec(2, 1,height_ratios=[1,4])

    ax1 = pl.subplot(gs[0])
    ax2 = pl.subplot(gs[1])

    poses = len(seq_MIf)
    entropy = []; frequencies = []; coverage=[]
    for dat in np.arange(len(sub_sels)):
        temp_MI = sub_matF[ref_sub[ref_sub[len(sub_matF)-1] == dat].index].iloc[0:-1]
        if bootstrap:
            boot_entropy = []; boot_frequencies = []; boot_cov = []
            for i in np.arange(boots):
                re_MI = resample(np.transpose(np.array(temp_MI)))
                entropy_pre,freq_pre,cov_pre = aims.calculate_shannon(re_MI)
                boot_entropy.append(entropy_pre); boot_frequencies.append(freq_pre)
                boot_cov.append(cov_pre)
            ent_avg = np.average(boot_entropy,axis=0)
            ent_std = np.std(boot_entropy,axis=0)
            freq_avg = np.average(boot_frequencies,axis=0)
            cov_avg = np.average(boot_cov,axis=0)
            ax2.plot(ent_avg,marker='o',linewidth=2.5,color=colors[dat])
            pl.fill_between(np.arange(len(ent_avg)),ent_avg+ent_std,ent_avg-ent_std,alpha=0.3,color=colors[dat])
            entropy.append(ent_avg); frequencies.append(freq_avg); coverage.append(1-cov_avg)
        else:
            entropy_pre,freq_pre,cov_pre = aims.calculate_shannon(np.transpose(np.array(temp_MI)))
            ax2.plot(entropy_pre,marker='o',linewidth=2.5,color=colors[dat])
            entropy.append(entropy_pre); frequencies.append(freq_pre); coverage.append(1-cov_pre)

    ax1.imshow(coverage,aspect='auto',interpolation='nearest',cmap='Greys')

    legend_elements=[]
    for j in np.arange(len(sub_sels)):
        element = [Line2D([0], [0], marker='o', color='w', label=label[j],markerfacecolor=colors[j], markersize=10)]
        legend_elements = legend_elements+element
    pl.legend(handles=legend_elements,ncol=len(sub_sels))

    pl.xlabel('Sequence Position'); pl.ylabel('Shannon Entropy (Bits)')

    if type(mat_size) != int:
        for i in np.arange(len(mat_size)-1):
            ax2.plot( (mat_size[i] + sum(mat_size[:i]) - 0.5) * np.ones(100),np.linspace(0,4.2,100),'black',linewidth = 3)
    pl.savefig(outputDir+'/AIMS_entropy.pdf',format='pdf')
    pl.close()

    ############# Do stats? #############
    if DOstats:
        data1 = sub_matF[ref_sub[ref_sub[len(sub_matF)-1] == 0].index].iloc[0:-1]
        data2 = sub_matF[ref_sub[ref_sub[len(sub_matF)-1] == 1].index].iloc[0:-1]
        p_ent = aims.do_statistics(data1,data2,num_reps=1000,test='function',test_func=aims.calculate_shannon,func_val = 0)

        fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(10,8))
        pl.plot(p_ent,color='black',marker='o',linewidth=2.5)
        pl.plot(np.arange(len(p_ent)),np.ones(len(p_ent))*0.05,color='red',linewidth=3)

        pl.xlabel('Sequence Position')
        pl.ylabel('p-value')

        pl.savefig(outputDir+'/aims_entropy_pval.pdf',format='pdf')
        pl.close()


    # # NOTE: The Mutual Information Calculation Can Be Quite Slow
    # Especially the case for MSAs with many amino acids. 

    # In[22]:


    # And then the mutual information:
    fig, ax = pl.subplots(1, len(sub_sels),squeeze=False,figsize=(14,5*len(sub_sels)))

    poses = len(seq_MIf)
    MI = []; ent_cond = []; count = []; MIMax = []
    for dat in np.arange(len(sub_sels)):
        temp_MI = sub_matF[ref_sub[ref_sub[len(sub_matF)-1] == dat].index]
        MI_pre,ent_cond_pre,count_pre = aims.calculate_MI(np.transpose(np.array(temp_MI)))
        MIMax.append(np.max(MI_pre))
        MI.append(MI_pre); ent_cond.append(ent_cond_pre); count.append(count_pre)
        ax[0,dat].set_title(label[dat])
        
    maxmax = np.max(MIMax)
    for dat in np.arange(len(sub_sels)):
        ax[0,dat].imshow(MI[dat],vmin=0,vmax=maxmax,cmap=cm.Greys)

    # Help Guide the eyes a bit
    if type(mat_size) != int:
        for i in np.arange(len(mat_size)-1):
            for j in np.arange(len(sub_sels)):
                ax[0,j].plot( (mat_size[i] + sum(mat_size[:i]) - 0.5) * np.ones(100),np.linspace(0,poses,100),'black',linewidth = 3)
                ax[0,j].plot( np.linspace(0,poses,100), (mat_size[i] + sum(mat_size[:i]) - 0.5) * np.ones(100) ,'black',linewidth = 3)

    pl.savefig(outputDir+'/AIMS_MI.pdf',format='pdf')
    pl.close()


    #####################################################
    # Amino acid frequency comparison for each individual dataset

    # Calculate the probabilities of seeing each amino acid at each position
    fig, ax = pl.subplots(1, 2,squeeze=False,figsize=(18,10))
    #pl.title(str(label[0])+ ' AA Frequency - ' + str(label[1]) + ' AA Frequency')

    AA_key=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

    freq1 = pandas.DataFrame(frequencies[0][:,1:])
    freq2 = pandas.DataFrame(frequencies[1][:,1:])
    # Remember that the "frequencies" are calculated with your custom
    # key in mind! So you need to carry that down here
    if custom_key:
        freq1.columns = my_AA_key
        freq2.columns = my_AA_key 
        fin_key = my_AA_key
    else:
        freq1.columns = AA_key
        freq2.columns = AA_key
        fin_key = AA_key

    x=ax[0,0].pcolormesh(freq1,vmin=0,vmax=0.25,cmap=cm.Greys)
    x=ax[0,1].pcolormesh(freq2,vmin=0,vmax=0.25,cmap=cm.Greys)

    #pl.colorbar(x); pl.ylabel('Sequence Position')
    xax=pl.setp(ax,xticks=np.arange(20)+0.5,xticklabels=fin_key)

    place=0
    if type(mat_size) == int:
        pl.plot(np.arange(21),place*np.ones(21),'black')
    else:
        for i in mat_size:
            place += i
            ax[0,0].plot(np.arange(21),place*np.ones(21),'black')
            ax[0,1].plot(np.arange(21),place*np.ones(21),'black')

    ax[0,0].set_xlabel("Amino Acid")
    ax[0,1].set_xlabel("Amino Acid")
    ax[0,0].set_ylabel("Sequence Position")
    ax[0,1].set_ylabel("Sequence Position")

    #pl.colorbar(x)
    pl.savefig(outputDir+'/AIMS_freq_sides.pdf',format='pdf')
    pl.close()

    # # Section 15: Binary Comparisons Between Datasets
    # Anything from here on requires a comparison between only two datasets. This can be two clusters, two loaded in files, two metadata subsets, whatever.
    # 
    # You can continue on from the previous sections even if you were analyzing multiple datasets. The analysis will just, by default, compare only the first two datasets

    # In[23]:


    # A bit easier to look at the DIFFERENCE in mutual information:
    fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(14,10))
    max_temp = np.max(MI[0]-MI[1]); min_temp = np.min(MI[0]-MI[1])
    if abs(min_temp) > abs(max_temp):
        propMin = min_temp
        propMax = -min_temp
    elif abs(max_temp) > abs(min_temp):
        propMin = -max_temp
        propMax = max_temp
    else:
        propMin = min_temp
        propMax = max_temp
    x = pl.imshow(MI[0] - MI[1], cmap=cm.PuOr, vmin = propMin, vmax = propMax)
    pl.colorbar(x); pl.title(str(label[0])+ ' MI - ' + str(label[1]) + ' MI')
    # Help Guide the eyes a bit
    if type(mat_size) != int:
        for i in np.arange(len(mat_size)-1):
            ax[0,0].plot((mat_size[i] + sum(mat_size[:i]) - 0.5) * np.ones(100),np.linspace(0,poses,100),'black',linewidth = 3)
            ax[0,0].plot( np.linspace(0,poses,100), (mat_size[i] + sum(mat_size[:i]) - 0.5) * np.ones(100) ,'black',linewidth = 3)

    pl.xlabel('Sequence Position'); pl.ylabel('Sequence Position')
    pl.savefig(outputDir+'/AIMS_MIdiff.pdf',format='pdf')
    pl.close()

    ###############################################################
    # Do stats?
    if DOstats:
        data1 = sub_matF[ref_sub[ref_sub[len(sub_matF)-1] == 0].index].iloc[0:-1]
        data2 = sub_matF[ref_sub[ref_sub[len(sub_matF)-1] == 1].index].iloc[0:-1]
        p_MI = aims.do_statistics(data1,data2,num_reps=MIboots,test='function',test_func=aims.calculate_MI,func_val = 0)

        dim1, dim2 = np.shape(p_MI)
        alpha = 0.05

        empty_mat_mi = np.zeros((dim1,dim2))
        for a in np.arange(dim1):
            for b in np.arange(dim2):
                if p_MI[a,b] < alpha:
                    empty_mat_mi[a,b] = 1

        fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(10,8))

        x= pl.imshow(empty_mat_mi,cmap=cm.Greys)
        pl.xlabel('Sequence Position')
        pl.ylabel('Sequence Position')

        pl.savefig(outputDir+'/AIMS_MIdiff_sig.pdf',format='pdf')
        pl.close()


    # In[24]:


    # Calculate the probabilities of seeing each amino acid at each position
    fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(16,8))
    pl.title(str(label[0])+ ' AA Frequency - ' + str(label[1]) + ' AA Frequency')
    freqMax = np.max(freq1.values-freq2.values); freqMin = np.min(freq1.values-freq2.values)
    freqBound = max(abs(freqMax),abs(freqMin))

    x=ax[0,0].pcolormesh(freq1-freq2,vmin=-freqBound,vmax=freqBound,cmap=cm.PuOr)
    AA_key=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    pl.colorbar(x); pl.ylabel('Sequence Position')
    xax=pl.setp(ax,xticks=np.arange(20)+0.5,xticklabels=fin_key)

    place=0
    if type(mat_size) == int:
        pl.plot(np.arange(21),place*np.ones(21),'black')
    else:
        for i in mat_size:
            place += i
            pl.plot(np.arange(21),place*np.ones(21),'black')

    pl.savefig(outputDir+'/AIMS_freqDiff.pdf',format='pdf')
    pl.close()

    if DOstats:
        data1 = sub_matF[ref_sub[ref_sub[len(sub_matF)-1] == 0].index].iloc[0:-1]
        data2 = sub_matF[ref_sub[ref_sub[len(sub_matF)-1] == 1].index].iloc[0:-1]
        p_freq = aims.do_statistics(data1,data2,num_reps=1000,test='function',test_func=aims.calculate_shannon,func_val = 1)

        dim1, dim2 = np.shape(p_freq)
        alpha = 0.05

        AA_key=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

        empty_mat = np.zeros((dim1,dim2))
        for a in np.arange(dim1):
            for b in np.arange(dim2):
                if p_freq[a,b] < alpha:
                    empty_mat[a,b] = 1

        # Need stuff to optionally shift around the x-axis:
        empty_frame = pandas.DataFrame(empty_mat[:,1:])
        if custom_key:
            empty_frame.columns = my_AA_key
        else:
            empty_frame.columns = AA_key

        fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(10,8))
        x = pl.pcolormesh(empty_frame,cmap=cm.Greys)
        xax=pl.setp(ax,xticks=np.arange(20)+0.5,xticklabels=fin_key)
        pl.ylabel('Sequence Position')
        pl.xlabel('Amino Acid')
        pl.savefig(outputDir+'/aa_diff_statSig.pdf',format='pdf')
        pl.close()


    # # Alright Now We Need to Bring Back in Linear Discriminant Analysis

    # In[25]:


    bigass1 = full_big_re.loc[ref_sub[ref_sub[len(sub_matF)-1] == 0].index]
    seq1_len = len(bigass1)
    bigass1.index = np.arange(seq1_len)
    bigass2 = full_big_re.loc[ref_sub[ref_sub[len(sub_matF)-1] == 1].index]
    seq2_len = len(bigass2)
    bigass2.index = np.arange(seq2_len)

    ######### Run the actual classification here ###############
    # Change the matSize variable to alter the number of vectors used to generate a classification
    bigF,weights,acc_all,mda_all,final,top_names = classy.do_linear_split(bigass1,bigass2,got_big=True,matSize=matSize)
    ############################################################

    import seaborn as sns
    fig = pl.figure(figsize = (12, 12))
    dset = ["Linear Discriminant Analysis" for x in range(seq1_len+seq2_len)]
    reacts = [label[0] for x in range(seq1_len)] + [label[1] for x in range(seq2_len)]

    d1 = {'Dataset': dset, 'Linear Discriminant 1': mda_all.reshape(len(mda_all)),
        'Reactivity' : reacts}
    df1 = pandas.DataFrame(data=d1)
    sns.set(style="white", color_codes=True,font_scale=1.5)
    sns.swarmplot(x="Dataset", y="Linear Discriminant 1", data=df1, hue = 'Reactivity', palette = "Dark2")
    print("Classification Accuracy")
    print(acc_all)
    pl.savefig(outputDir+'/AIMS_LDA.pdf',format='pdf')
    pl.close()


    # In[26]:


    # Show the top properties that differentiate the two populations
    # show_top = how many of these top values do you want to show? don't recommend more than ~5
    # solely due to how busy the figure gets
    # Again, see eLife paper for biophysical property definitions
    import seaborn as sns
    show_top = 5
    dset_parse = final[top_names[0:show_top]]
    dset_ID_pre1 = bigF['ID']
    dset_ID_pre2 = dset_ID_pre1.replace(1.0,label[0])
    dset_ID = dset_ID_pre2.replace(2.0,label[1])

    bigass_parse_dset = pandas.concat([dset_parse,dset_ID],axis = 1)
    try:
        sns.pairplot(bigass_parse_dset,hue = 'ID')
        pl.savefig(outputDir+'/AIMS_pairplot.pdf',format='pdf')
        pl.close()
    except:
        print("Pairplot failed")


    # # Detailed Amino Acid Frequency Breakdowns
    # This feature is most useful in peptide analysis, but could be useful for general repertoire comparisons

    # In[27]:


    # Need to get back just our sequences of interest.
    # Need to get back just our sequences of interest.
    seq1 = sub_seqF.loc[ref_sub[ref_sub[len(sub_matF)-1] == 0].index]
    seq2 = sub_seqF.loc[ref_sub[ref_sub[len(sub_matF)-1] == 1].index]
    # Calculate both position-insensitive amino acid frequency and digram frequencies
    # Can either normalize to the # of sequences or total number of AA (num_seq or num_AA)
    AA_freq_all1, digram_all1 = aims.full_AA_freq(seq1,norm='num_AA')
    AA_freq_all2, digram_all2 = aims.full_AA_freq(seq2,norm='num_AA')

    freqAll1 = np.transpose(pandas.DataFrame(AA_freq_all1))
    freqAll1.columns = AA_key
    freqAll2 = np.transpose(pandas.DataFrame(AA_freq_all2))
    freqAll2.columns = AA_key

    freqAll1_df = freqAll1[fin_key]
    freqAll2_df = freqAll2[fin_key]

    fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(16,8))
    ax[0,0].bar(np.arange(len(AA_freq_all1)),freqAll1_df.values[0], color=colors[0],alpha=0.5)
    ax[0,0].bar(np.arange(len(AA_freq_all2)),freqAll2_df.values[0],color=colors[1],alpha=0.5)
    xax=pl.setp(ax,xticks=np.arange(20),xticklabels=fin_key)
    ax[0,0].legend([label[0],label[1]])
    pl.ylabel('Frequency')
    pl.xlabel('Amino Acid')
    pl.savefig(outputDir+'/AIMS_AAnetProb.pdf',format='pdf')
    pl.close()

    #########################################################
    if DOstats:
        data1 = seqF[ref_sub[ref_sub[len(sub_matF)-1] == 0].index]
        data2 = seqF[ref_sub[ref_sub[len(sub_matF)-1] == 1].index]

        p_freq = aims.do_statistics(data1,data2,num_reps=boots,test='function',test_func=aims.full_AA_freq,func_val = 0)

        p_df = np.transpose(pandas.DataFrame(p_freq))
        p_df.columns = AA_key
        p_df_fin = p_df[fin_key]

        fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(16,8))
        pl.plot(p_df_fin.values[0],marker='o',linewidth=3,color='black')
        pl.plot(np.arange(20),np.ones(20)*0.05,linewidth=3,color='red')
        xax=pl.setp(ax,xticks=np.arange(20),xticklabels=fin_key)
        pl.savefig(outputDir+'/AIMS_AAnetProb_pval.pdf',format='pdf')
        pl.close()


    # # Plot the Digram Frequencies Per Sequence (or per AA)
    # There can occasionally be interesting information in these digram patterns, particularly if there are notabe differences between populations.
    # 
    # Note, while the matrix could possibly appear symmetric, it need not be so. The y-axis gives the first amino acid in the digram, the x-axis gives the second. So on the x,y coordinate map, E (x-axis) and D (y-axis) gives the frequency of the digram DE.

    # In[28]:

    fig, ax = pl.subplots(int(len(label)/2), 2,squeeze=False,figsize=(16,12))
    plot_max = np.max([np.max(digram_all1),np.max(digram_all2)])

    ax[0,0].set_title(str(label[0])+ ' Digram Frequency')
    ax[0,1].set_title(str(label[1])+ ' Digram Frequency')

    # Again for potentially altering the AA axis how you want:
    digAll1 = pandas.DataFrame(digram_all1)
    digAll1.columns = AA_key; digAll1.index = AA_key
    digAll2 = pandas.DataFrame(digram_all2)
    digAll2.columns = AA_key; digAll2.index = AA_key

    digF1 = digAll1[fin_key].loc[fin_key]
    digF2 = digAll2[fin_key].loc[fin_key]

    x1 = ax[0,0].imshow(digF1,  cmap = cm.Greys,vmin=0,vmax=plot_max)
    x2 = ax[0,1].imshow(digF2,  cmap = cm.Greys,vmin=0,vmax=plot_max)

    xax=pl.setp(ax[0,0],xticks=np.arange(20),xticklabels=fin_key)
    yax=pl.setp(ax[0,0],yticks=np.arange(20),yticklabels=fin_key)
    xax=pl.setp(ax[0,1],xticks=np.arange(20),xticklabels=fin_key)
    yax=pl.setp(ax[0,1],yticks=np.arange(20),yticklabels=fin_key)

    fig.colorbar(x1, ax=ax[0, 0], shrink=0.5)
    fig.colorbar(x2, ax=ax[0, 1], shrink=0.5)
    pl.savefig(outputDir+'/AIMS_digram_sep.pdf',format='pdf')
    pl.close()


    # In[29]:

    # Look at the difference in these digram frequencies
    # For now, hard to say if these are necessarily significant, 
    # more downstream analysis is needed to tease out conclusions
    fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(16,8))
    pl.title('Peptide Digram Difference')
    min_temp = np.min(digram_all2-digram_all1)
    max_temp = np.max(digram_all2-digram_all1)
    if abs(min_temp) > abs(max_temp):
        propMin = min_temp
        propMax = -min_temp
    elif abs(max_temp) > abs(min_temp):
        propMin = -max_temp
        propMax = max_temp
    else:
        propMin = min_temp
        propMax = max_temp
    zzz = pl.imshow(digF1-digF2, vmin = propMin, vmax = propMax, cmap = cm.PuOr)
    xax=pl.setp(ax,xticks=np.arange(20),xticklabels=fin_key)
    yax=pl.setp(ax,yticks=np.arange(20),yticklabels=fin_key)
    ax[0,0].set_ylabel('First Amino Acid')
    ax[0,0].set_xlabel('Second Amino Acid')
    pl.colorbar(zzz)
    pl.savefig(outputDir+'/AIMS_digramDiff.pdf',format='pdf')
    pl.close()

    ##################################################
    # Do statistics?
    if DOstats:
        data1 = seqF[ref_sub[ref_sub[len(sub_matF)-1] == 0].index]
        data2 = seqF[ref_sub[ref_sub[len(sub_matF)-1] == 1].index]

        p_digram = aims.do_statistics(data1,data2,num_reps=1000,test='function',test_func=aims.full_AA_freq,func_val = 1)

        dim1, dim2 = np.shape(p_digram)
        alpha = 0.05

        empty_mat_digram = np.zeros((dim1,dim2))
        for a in np.arange(dim1):
            for b in np.arange(dim2):
                if p_digram[a,b] < alpha:
                    empty_mat_digram[a,b] = 1

        fin_key = my_AA_key

        empty_dig = pandas.DataFrame(empty_mat_digram)
        empty_dig.columns = AA_key; empty_dig.index = AA_key
        empty_digF = empty_dig[fin_key].loc[fin_key]

        fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(10,8))
        pl.title('Peptide Digram Difference Stat. Sig.')
        x= pl.imshow(empty_digF,cmap=cm.Greys)
        xax=pl.setp(ax,xticks=np.arange(20),xticklabels=fin_key)
        xax=pl.setp(ax,yticks=np.arange(20),yticklabels=fin_key)
        ax[0,0].set_ylabel('First Amino Acid')
        ax[0,0].set_xlabel('Second Amino Acid')
        pl.savefig(outputDir+'/digram_diff_sig.pdf',format='pdf')
        pl.close()

if __name__ == '__main__':
    x=run()
    print(x)