##################################################################################################################
# The goal of this new file is to "hide" a bunch of the crap in the AIMS notebook
# I think that if people want to edit any of it, maybe they can pull it out in the tutorials
# or directly mess with this file.
# My hope is that the focus in the notebook can be more about plotting and making the data look pretty
# and/or adding to it in some meaningful way.

# As of now, these functions are ONLY called in the notebook. May eventually link the CLI
# and the GUI back to these, but not just yet [08/26/26]...
##################################################################################################################
import numpy as np
import pandas
import matplotlib.pyplot as pl
from matplotlib.lines import Line2D
from matplotlib import cm
import matplotlib.gridspec as gridspec
from aims_immune import aims_loader as aimsLoad
from aims_immune import aims_analysis as aims
from aims_immune import aims_classification as classy
from numba import njit, prange, set_num_threads, get_num_threads
from sklearn.utils import resample
import random

# Define colormap for plotting functions
import matplotlib as mpl
upper = mpl.cm.jet(np.arange(256))
lower = np.ones((int(256/4),4))
for i in range(3):
    lower[:,i] = np.linspace(1, upper[0,i], lower.shape[0])
cmap = np.vstack(( lower, upper ))
cmap = mpl.colors.ListedColormap(cmap, name='myColorMap', N=cmap.shape[0])

# This is where we actually do all the loading of the data:
def loadDat(fileName,datName,datDir,drop_duplicates,subStart=[],subEnd=[],subset=False):
    for i in np.arange(len(fileName)):
        seq_pre = aimsLoad.seq_loader(datDir+'/'+fileName[i],label=datName[i],subset=subset,
                                      drop_dups= drop_duplicates,subset_starts=subStart,subset_ends=subEnd)
        # Need to reset indices because otherwise they won't actually stack...
        seq_pre.index = np.arange(np.shape(seq_pre)[0])
        if i == 0:
            seqF = seq_pre
        else:
            seqF = pandas.concat([seqF,seq_pre],axis=1)

    # This is where we define the matrix dimensions
    mat_size = aims.get_sequence_dimension(seqF)

    # General changes that need to be done for every type of molecule
    AA_num_key = aims.get_props()[1]
    if np.shape(seqF)[0] != 1:
        for i in np.arange(len(mat_size)):
            if i == 0:
                xtick_loc = [mat_size[i]/2]
            else:
                pre_loc = sum(mat_size[:i])
                xtick_loc = xtick_loc + [mat_size[i]/2 + pre_loc]
    else:
        xtick_loc = np.array(mat_size)/2

    return(seqF,xtick_loc,AA_num_key,mat_size)

def get_keys(custom_key=''):
    if custom_key!='':
        if len(custom_key)!=20:
            print('Wrong number AAs! Switch to default')
            my_AA_key=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
        #rev_suggest = 'WYFMLIVAPGCSTNQDEHRK'
        rev_altered = 'WFMLIVPYHAGSTDECNQRK'
        my_AA_key = [a for a in rev_altered]
    else:
        # My AA key here is the "standard" AIMS key that has been used in previous papers
        # Note changing the key doesn't change anything BUT the ordering of Amino acids in some figures
        my_AA_key=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

    # Just in case you want to do MSA analysis, you need to use "my_AA_key_dash"
    my_AA_key_dash = my_AA_key + ['-']
    return(my_AA_key,my_AA_key_dash)

def grab_big(seqF,dsetF,seq_MIf,datName,normalize,renormalize,molecule,
             my_AA_key,my_AA_key_dash,mat_size,align,pad,nCores):

    if molecule.lower() == 'ig':
        special = ''
    elif molecule.lower() == 'peptide':
        special = ''
    elif molecule.lower() == 'msa':
        special = 'MSA'

    bigass = classy.get_bigass_matrix(dsetF,AA_key=my_AA_key,AA_key_dash=my_AA_key_dash, 
                                      giveSize = mat_size, alignment = align, norm = normalize,
                                      special=special,bulge_pad=pad,nCores=nCores)
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
            if str(j).find(datName[i]+'_') != -1:
                tokenized_dset.append(i)
    token_df = pandas.DataFrame(tokenized_dset,columns=['ID'])
    IDed_full_big = pandas.concat([full_big,token_df],axis=1)

    # Lastly, create a good-ole traditional averaged bphys property matrix
    # i.e. each sequence gets a single value for averaged charge, averaged flexibility, etc...
    posLen,seqLen = np.shape(seq_MIf)
    # The 61 is hardcoded here because it is our number of properties. Eventually we will let users define which properties to use
    seq_bigReshape = np.array(full_big).reshape(seqLen,61,posLen)
    return(full_big,parsed_mat,NonNorm_big,IDed_full_big,seq_bigReshape,token_df)

def do_cluster(chosen_dset,seq_MIf,reduce,clust,NClust=100,
               min_samples = 10,eps=0.5,nComp=3,n_neighbors=25,
               state=617):
    # As a note, eps is only used for DBSCAN, and HDBSCAN actually scans
    # across various eps values to find the best representation
    # (so actually, HDBSCAN should be the preferred way of doing things)
    if reduce == 'pca':
        from sklearn.decomposition import PCA
        pca = PCA(n_components=3, svd_solver='full')
        final=pca.fit_transform(chosen_dset)
        transform = pandas.DataFrame(np.transpose(final),columns = seq_MIf.columns)
        print("PCA Explained Variance Ratio:")
        print(pca.explained_variance_ratio_)
    elif reduce == 'umap':
        import umap
        reducer = umap.UMAP(n_components=nComp, n_neighbors = n_neighbors,random_state=state)
        final = reducer.fit_transform(chosen_dset)
        transform = pandas.DataFrame(np.transpose(final),columns = seq_MIf.columns)

    # Cluster the results:
    import sklearn.cluster as cluster
    clust_input = np.array(np.transpose(transform))
    if clust == 'kmean':
        clusts = cluster.KMeans(n_clusters=NClust).fit_predict(clust_input)
    elif clust == 'optics':
        clusts = cluster.OPTICS(min_samples=min_samples,eps=eps).fit_predict(clust_input)
    elif clust == 'dbscan':
        clusts = cluster.DBSCAN(min_samples = min_samples,eps=eps).fit_predict(clust_input)
    elif clust == 'hdbscan':
        # Note, this only exists in sklearn 1.3 and on
        clusts = cluster.HDBSCAN(min_cluster_size = min_samples).fit_predict(clust_input)

    cluster_dset = pandas.DataFrame(clusts,columns=['cluster'])
    return(clust_input,cluster_dset)

# Need to allow for an empty cluster dataset for when you skip clustering in the notebook (or elsewhere)
def get_metadata(token_df,datName,cluster_dset=[],compile_index=[],metaPath='',meta_form='category',gotMeta=False):
    # We can incorporate metadata either defining a categorical map or a quantitative map
    ###################################################################################
    # meta_form is either "category" or "quant". Check ReadTheDocs if more descriptions are needed
    # If using default metadata (i.e. loaded file), it should be "category"
    # GOT META SHOULD BE TRUE ONLY IF YOU HAVE A FILE LIKE "VDJDB META PARSE"
    ###################################################################################

    if meta_form == 'category':
        if gotMeta:
            ###################################################################################
            # This is an example of how to load in your own metadata
            # Your # metadat entries should match your # sequence entries
            metapre = pandas.read_csv(metaPath)
            metapre.index = compile_index
            meta = metapre.loc[compile_index]
            # A bunch of columns to choose from. Note though, you can
            # only visualize one at a time in the current version
            metadat = pandas.DataFrame(meta['antigen.species'])
            metadat.columns = [0]
            ###################################################################################
        else:
            ###################################################################################
            # This is the easiest default metadata definition.
            # Just based on the files that were loaded in (useless if only 1 file)
            metapre = token_df
            # Convert the metadata from numbers to strings
            meta_conv = []
            for i in token_df.values:
                meta_conv.append(datName[i[0]])
            metadat = pandas.DataFrame(meta_conv)
            ###################################################################################
    elif meta_form == 'quant':
        # QUANTITATIVE METADATA IS NOT READY YET! ON THE TO-DO!
        # This "quantitative" metadata will just count from 1 to the length of the dataset,
        # coloring the points on the plot in order.
        metaPRE = np.arange(len(token_df))
        metadat = pandas.DataFrame(metaPRE)
        # But, you could add things like MFI, binding affinity, or GEX data for a particular gene

    # From there, not much should change unless you want to give a unique name to your metadata
    meta_name = 'metadata'
    meta_map = aims.encode_meta(metadat)
    meta_map.columns = [meta_name]
    meta_leg = metadat.drop_duplicates().values

    # We also need to define our clusters more clearly
    # Need to pull some of these values even when you aren't doing clustering.
    clust_name='cluster'
    if len(cluster_dset) == 0:
        clust_map = []
        clust_leg = []
    else:
        clust_map = cluster_dset
        clust_leg = [a[0] for a in cluster_dset.drop_duplicates().sort_values(clust_name).values]
    return(clust_map,clust_leg,clust_name,meta_map,meta_leg,meta_name,metadat)

def back_decode(seq_MIf,outputDir,my_AA_key,AA_num_key,outName='AIMS_encodedSeqs.csv'):
    seqT = np.transpose(seq_MIf)
    fin_convert = []
    for i in np.arange(len(seqT)):
        tt = aims.decode_mat(seqT.iloc[i],num_key_AA=AA_num_key,key_AA=my_AA_key)
        fin_convert = fin_convert + [*tt]
    fff = np.array(fin_convert).reshape(len(seqT),len(tt))
    pandas.DataFrame(fff).to_csv(outputDir+'/'+outName)

def plot_clusters(clust_input,plot_metas,clust,
                  show_labels=False,proj_show='both',clust_show='both',
                  need_meta=True,labels=[]):

    # Grouped all these things into a single file since they come 
    # from a singular function
    clust_map=plot_metas[0];clust_leg=plot_metas[1]
    clust_name=plot_metas[2]; meta_map=plot_metas[3]
    meta_leg=plot_metas[4];meta_name=plot_metas[5]
    ##########################################################################
    ######### Optionally show labels to identify outlier sequences or sequences of interest##################
    # NOTE that these labels really only work well with very few datapoints. Future updates will try to
    # work around this issue a little bit better...
     # Could optionally plot other data or change legends if you would like
    chosen_map1 = clust_map; leg1 = clust_leg
    chosen_map2 = meta_map; leg2 = meta_leg

    # Define a colormap to color metadata. Can change this if you want, look at matplotlib colormap options
    # BUT NOTE, if your colormap is discrete, you can get repetitive colors in the plot if #meta_unique > #colors
    cmapF = pl.get_cmap('jet')
    # So there's no reason to have a function for this, other than to make this notebook look a bit prettier
    # Really just defining a bunch of stuff repeatedly with if statements
    fig3d,plotloc,plottype,plotem,legends,dattype = aims.get_plotdefs(clust_show,proj_show,chosen_map1,chosen_map2,leg1,leg2)

    # Now plot the actual stuff, and save the object handles in a list
    colorhandle= []; ax=[]
    for i in np.arange(len(plotloc)):
        if dattype[i] == 'clust':
            if clust == 'kmean':
                cmap_use = pl.get_cmap('Pastel1')
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
            ax[i].set_xlabel('AX1'); ax[i].set_ylabel('AX2')

    # We use the metadata mapped colors down the line, so if you aren't plotting it then you need to save some other way
    for i in np.arange(len(legends)):
        cmap_pre = pandas.DataFrame(colorhandle[i].get_facecolors())
        if dattype[i] == 'clust':
            mapDF = pandas.concat([clust_map,cmap_pre],axis=1)
            mapped_colors= mapDF.sort_values('cluster').drop_duplicates('cluster').values[:,1:]
        else:
            mapDF = pandas.concat([meta_map,cmap_pre],axis=1)
            mapped_colors= mapDF.sort_values('metadata').drop_duplicates('metadata').values[:,1:]
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
        cmap_discrete = cmapF(np.linspace(0, 1, len(meta_map.drop_duplicates())))
        cmap_pre = pandas.DataFrame(cmap_discrete)
        mapDF = pandas.concat([meta_map.drop_duplicates().reset_index(),cmap_pre],axis=1)
        mapped_colors= mapDF.sort_values('metadata').drop_duplicates('metadata').values[:,2:]
        keep_map = mapped_colors

    if show_labels:
        # ONLY show this for a 2D plot
        for num in np.arange(len(ax)):
            if plottype[num] == '2d':
                break
        a = 0; plot1 = clust_input[:,0]; plot2 = clust_input[:,1]
        plot_labels = labels
        for i,j in zip(plot1,plot2):
            ax[num].annotate(str(plot_labels[a]),xy=(i,j),fontsize=14)
            a+=1

    if need_meta:
        return(fig3d,ax,keep_map)
    else:
        return(fig3d,ax)

def plot_quantClust(cluster_dset,keep_map,plot_metas,norm=True,cmap3=pl.get_cmap('tab20b')):
    # Grouped all these things into a single file since they come 
    # from a singular function
    clust_map=plot_metas[0];clust_leg=plot_metas[1]
    clust_name=plot_metas[2]; meta_map=plot_metas[3]
    meta_leg=plot_metas[4];meta_name=plot_metas[5]

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

    meta_legF = np.array(meta_leg).reshape(len(meta_leg))
    if len(meta_leg) < 6:
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

    return(fig,final_breakdown,meta_legF)

def plotSubsets(seq_bigReshape,seq_MIf,plot_metas,mat_size,
                subset_sel='cluster',plot_props=False,showProp=1,
                show_lines = False):
    clust_map=plot_metas[0];clust_leg=plot_metas[1]; clust_name=plot_metas[2]; 
    meta_map=plot_metas[3]; meta_leg=plot_metas[4];meta_name=plot_metas[5]
    if subset_sel.lower() == 'cluster':
        chosen_map = clust_map; chosen_name = clust_name
    elif subset_sel.lower() == 'metadata':
        chosen_map = meta_map; chosen_name = meta_name
    else:
        print('ERROR: subset_sel = cluster or metadata')
        return()
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
    return(fig,chosen_map,chosen_name)

# Need to create a simplified version of the sequence distance calculation
# Previous versions used pandas dataframes too much... Very simple here
@njit(parallel=True)
def newParallel(full_big,num_threads=-1):
    # Need to allow an option such that users can opt for NO parallelization.
    # Will have parallel default now that we aren't using multiprocessing.
    # I trust numba to handle things moreso than my scrappy code.
    orig_threads = get_num_threads()
    # By default, numba uses all threads, so if user does not define num threads we will do the same
    if num_threads == -1:
        num_threads = orig_threads
    elif num_threads > orig_threads:
        print("Warning: You have requested more threads than are available. Setting to max threads.")
        num_threads = orig_threads  
    # If user defines number threads, set em.
    set_num_threads(num_threads)
    dist_calc = np.zeros((len(full_big),len(full_big)))
    for i in prange(len(full_big)):
        for j in prange(len(full_big)):
            dist_calc[i,j] = np.sqrt(sum((full_big[i] - full_big[j])**2))
            
    return(dist_calc)

def run_AIMSdist(full_big,seq_MIf,seqF,plot_metas,chosen_map,chosen_name,
                 parallel_dist=True,get_distClusts=True,maxD=5,nThreads=-1,reorder=False,noFig=False):
    clust_map=plot_metas[0];clust_leg=plot_metas[1]; clust_name=plot_metas[2]; 
    meta_map=plot_metas[3]; meta_leg=plot_metas[4];meta_name=plot_metas[5]
    metadat=plot_metas[6]
    for i in chosen_map.sort_values(chosen_name).drop_duplicates().values:
        if i == -1:
            continue
        sub_MI_temp = seq_MIf[seq_MIf.columns[chosen_map[chosen_map[chosen_name] == i[0]].index]]
        sub_seqs_temp = np.transpose(seqF[sub_MI_temp.columns])
        if i == 0:
            sorted_seqs = sub_seqs_temp
        else:
            sorted_seqs = pandas.concat([sorted_seqs,sub_seqs_temp])

    if reorder:
        dists=newParallel(full_big.loc[sorted_seqs.index].values,num_threads = nThreads)
    else:
        dists=newParallel(full_big.values,num_threads=nThreads)
    if noFig:
        fig = []
    else:
        fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(10,8))
        x = pl.imshow(np.transpose(dists), interpolation='nearest', aspect='auto')
        pl.colorbar(x)
    # Optionally can get back distance clusters:
    if get_distClusts:
        distance_clusters = aims.get_distClusts(dists,metadat,max_d=maxD)
        if noFig:
            return(dists,distance_clusters)
        else:
            return(fig,dists,distance_clusters)
    elif noFig:
        return(dists)
    else:
        return(fig,dists)

# For now, we're going to deprecate "seqLogo" but leave it here for 
def viz_subs(dset,seq_MIf,meta_legF,sub_sels,chosen_map,chosen_name,
             seqlogo=False,save_subSeqs=False,saveAll=False):
    if saveAll:
        # overwrite sub_sels if you want to visualize everything
        sub_sels = np.arange(len(chosen_map.drop_duplicates()))
    
    fig, ax = pl.subplots(len(sub_sels), 1,squeeze=False,figsize=(16,4*len(sub_sels)))
    label=[]
    a = 0
    for i in sub_sels:
        # Look at umap dset or pca dset
        sub_MI = seq_MIf[seq_MIf.columns[chosen_map[chosen_map[chosen_name] == i].index]]
        sub_seqs = np.transpose(dset[sub_MI.columns])
        if chosen_name.lower() == 'metadata' or chosen_name.lower()=='meta':
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
        return(fig,sub_matF,label,sub_seqF)
    else:
        return(fig,sub_matF,label)

@njit(parallel=True)
def netAvg_bootStrap(temp_bigReshape,boots,prop,num_threads=-1):
    # all this is standard across my numba scripts
    # search elsewhere in the code for more commments on code
    orig_threads = get_num_threads()
    if num_threads == -1:
        num_threads = orig_threads
    elif num_threads > orig_threads:
        print("Warning: You have requested more threads than are available. Setting to max threads.")
        num_threads = orig_threads  
    set_num_threads(num_threads)
    prop_avg = np.zeros(boots)
    re_big = np.zeros(np.shape(temp_bigReshape))
    for i in prange(boots):
        # Full dataset resample each time:
        for j in range(len(temp_bigReshape)):
            # Generate random index from 0 to n-1
            idx = np.random.randint(0, len(temp_bigReshape))
            re_big[j] = temp_bigReshape[idx]
        # Here is our resampling now
        prop_avg[i] = np.average(np.average(re_big[:,prop,:]))
    fin_avg = np.average(prop_avg)
    fin_std = np.std(prop_avg)
    return(fin_avg,fin_std)
# These need to be two separate functions because of how numba
# compiles the code. Isn't happy about ambiguous output sizes
@njit(parallel=True)
def posAvg_bootStrap(temp_bigReshape,boots,posLen,
                  prop,num_threads=-1):
    # all this is standard across my numba scripts
    # search elsewhere in the code for more commments on code
    orig_threads = get_num_threads()
    if num_threads == -1:
        num_threads = orig_threads
    elif num_threads > orig_threads:
        print("Warning: You have requested more threads than are available. Setting to max threads.")
        num_threads = orig_threads  
    set_num_threads(num_threads)
    # Start actual code
    prop_avg = np.zeros((boots,posLen))
    re_big = np.zeros(np.shape(temp_bigReshape))
    for i in prange(boots):
        # Full dataset resample each time:
        for j in range(len(temp_bigReshape)):
            # Generate random index from 0 to n-1
            idx = np.random.randint(0, len(temp_bigReshape))
            re_big[j] = temp_bigReshape[idx]
        # Here is our resampling now
        myBig = re_big[:,prop,:]
        result = np.zeros(myBig.shape[1])
        #result = np.average(myBig,axis=0)
        for k in np.arange(myBig.shape[1]):
            result[k] += np.transpose(myBig)[k].mean()
        prop_avg[i] = result
    # Unlike the net avg, just return the full matrix
    # we don't want to do more weird averages in numba
    return(prop_avg)


def netStats(ref_sub,sub_big,posLen,sub_matF,
             prop_names = ['Charge','Hydrophobicity','Bulkiness','Flexibility'],
             reps=1000):
    take_sub1 = ref_sub[ref_sub[len(sub_matF)-1] == 0].index
    take_big1 = sub_big.loc[take_sub1]
    temp_bigReshape1 = np.array(take_big1).reshape(len(take_big1),61,posLen)
    take_sub2 = ref_sub[ref_sub[len(sub_matF)-1] == 1].index
    take_big2 = sub_big.loc[take_sub2]
    temp_bigReshape2 = np.array(take_big2).reshape(len(take_big2),61,posLen)

    # Significance for Bar plots
    p_list = []
    for i in np.arange(4):
        propF = i+1
        data1 = np.average(temp_bigReshape1[:,propF,:],axis=1)
        data2 = np.average(temp_bigReshape2[:,propF,:],axis=1)
        p = aims.do_statistics(data1,data2,num_reps=reps,test='average')
        #print('p-value for '+prop_names[i]+' bar plot: '+str(p))
        p_list = p_list + [p]
    return(p_list)

def do_netAvg(sub_big,ref_sub,sub_matF,sub_sels,posLen,label,
              colors= ['Crimson','darkorchid'],
              bootstrap = False, boots = 1000,stats=False):
    fig, ax = pl.subplots(1, 1,squeeze=False,figsize=(16,8))

    # Since we aren't being as open with these functions, we need to have better handling for people.
    if len(colors)!=len(sub_sels):
        print('Define colors for each dataset if you do not want them auto-assigned')
        # Generate n evenly spaced values between 0 and 1
        # This should be consistent across all of our plots, since cmap
        # is defined at the top of this script
        colors = cmap(np.linspace(0.05, 1, len(sub_sels)))

    # For now, hard code the x-axis but definitely want to let
    # users pick which properties to plot at some point
    x_axis=[-0.2,0.9,2,3.1]
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
                fin_avg,fin_std = netAvg_bootStrap(temp_bigReshape,boots,propF)
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
    if stats:
        p_list = netStats(ref_sub,sub_big,posLen,sub_matF)
        return(fig,p_list)
    else:
        return(fig)

def do_posAvg(sub_sels,ref_sub,sub_matF,sub_big,posLen,mat_size,label,
              prop_sel = [1,2],colors=['Crimson','darkorchid'],
              bootstrap=False,boots=1000):
    # Now get the position sensitive avarege biophysical properties
    # DONT CHANGE THE FIGURE ASPECT RATIOS. AS-IS, IT SHOULD LINE UP VERY WELL WITH THE AIMS MATRIX VISUALIZATIONS
    fig, ax = pl.subplots(2, 1,squeeze=False,figsize=(16,10))
    full_avg = []; full_std = []

    if len(colors)!=len(sub_sels):
        print('Define colors for each dataset if you do not want them auto-assigned')
        colors = cmap(np.linspace(0.05, 1, len(sub_sels)))

    for dat in np.arange(len(sub_sels)):
        a=0
        for prop in prop_sel:
            take_sub = ref_sub[ref_sub[len(sub_matF)-1] == dat].index
            take_big = sub_big.loc[take_sub]
            temp_bigReshape = np.array(take_big).reshape(len(take_big),61,posLen)
            if bootstrap:
                pre_avg = posAvg_bootStrap(temp_bigReshape,boots,posLen,prop)
                fin_avg = np.average(pre_avg,axis=0)
                fin_std = np.std(pre_avg,axis=0)
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
    return(fig)

def posAvg_stats(ref_sub,sub_matF,sub_big,mat_size,posLen,prop_sel=[1,2]):
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

        p = aims.do_statistics(data1,data2,num_reps=1000,test='average')
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
    return(fig)

def pos_Shannon(seq_MIf,sub_sels,sub_matF,ref_sub,mat_size,label,
                colors=['Crimson','darkorchid'],
                bootstrap=False,boots=1000):
    fig = pl.figure(figsize=(16,8))
    gs = gridspec.GridSpec(2, 1,height_ratios=[1,4])
    if len(colors)!=len(sub_sels):
        print('Define colors for each dataset if you do not want them auto-assigned')
        colors = cmap(np.linspace(0.05, 1, len(sub_sels)))

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

    pl.legend(label); pl.xlabel('Sequence Position'); pl.ylabel('Shannon Entropy (Bits)')

    if type(mat_size) != int:
        for i in np.arange(len(mat_size)-1):
            ax2.plot( (mat_size[i] + sum(mat_size[:i]) - 0.5) * np.ones(100),np.linspace(0,4.2,100),'black',linewidth = 3)
    return(fig,entropy,frequencies)

def pos_MI(seq_MIf,sub_sels,sub_matF,ref_sub,label,mat_size):
    fig, ax = pl.subplots(1, len(sub_sels),squeeze=False,figsize=(18,5*len(sub_sels)))

    poses = len(seq_MIf)
    MI = []; ent_cond = []; count = []
    for dat in np.arange(len(sub_sels)):
        temp_MI = sub_matF[ref_sub[ref_sub[len(sub_matF)-1] == dat].index].iloc[0:-1]
        MI_pre,ent_cond_pre,count_pre = aims.calculate_MI(np.transpose(np.array(temp_MI)))
        ax[0,dat].imshow(MI_pre,vmin=0,vmax=2,cmap=cm.Greys)
        MI.append(MI_pre); ent_cond.append(ent_cond_pre); count.append(count_pre)
        ax[0,dat].set_title(label[dat])

    #pl.colorbar(x)
    # Help Guide the eyes a bit
    if type(mat_size) != int:
        for i in np.arange(len(mat_size)-1):
            for j in np.arange(len(sub_sels)):
                ax[0,j].plot( (mat_size[i] + sum(mat_size[:i]) - 0.5) * np.ones(100),np.linspace(0,poses,100),'black',linewidth = 3)
                ax[0,j].plot( np.linspace(0,poses,100), (mat_size[i] + sum(mat_size[:i]) - 0.5) * np.ones(100) ,'black',linewidth = 3)
    return(fig,MI,ent_cond)

def pos_freq(frequencies,my_AA_key,mat_size):
    # Calculate the probabilities of seeing each amino acid at each position
    fig, ax = pl.subplots(1, 2,squeeze=False,figsize=(18,10))
    #pl.title(str(label[0])+ ' AA Frequency - ' + str(label[1]) + ' AA Frequency')

    AA_key=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']

    freq1 = pandas.DataFrame(frequencies[0][:,1:])
    freq2 = pandas.DataFrame(frequencies[1][:,1:])
    # Remember that the "frequencies" are calculated with your custom
    # key in mind! So you need to carry that down here

    freq1.columns = my_AA_key
    freq2.columns = my_AA_key 
    fin_key = my_AA_key

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
    return(fig,freq1,freq2)