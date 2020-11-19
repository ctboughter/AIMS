import numpy as np
import pandas
import matplotlib.pyplot as pl
import math
import matplotlib as mpl
from matplotlib import cm

from sklearn.discriminant_analysis import LinearDiscriminantAnalysis
from sklearn.metrics import accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.decomposition import PCA
from sklearn.utils import resample
from sklearn.metrics import roc_curve as auc
from sklearn.model_selection import KFold
from sklearn.model_selection import StratifiedKFold
from sklearn.model_selection import LeaveOneOut
from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.decomposition import KernelPCA
from sklearn.feature_selection import SelectKBest
from sklearn.feature_selection import mutual_info_classif
from scipy.stats import pearsonr

# Custom Script
import aims_analysis as aims

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

def apply_matrix(mono_PCA,max_diffs,mat_size=100,props=properties[1:],ridZero=False,win_size = 3):
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
        mono_pca_NEW=np.transpose(mono_PCA)[~np.all(np.transpose(mono_PCA) == 0,axis=1)]
        mono_pca_NEW=np.transpose(mono_pca_NEW)
    else:
        mono_pca_NEW = mono_PCA

    # So this is where we should be able to do the averaging
    # Initialize the variable
    new_mat_mono=np.zeros((mat_size,len(mono_pca_NEW))) 
    for j in np.arange(len(mono_pca_NEW)): # for every clone
        for k in np.arange(len(max_diffs)): 
            for win_slide in np.arange(win_size):
                # This really shouldn't happen, but just in case...
                if max_diffs[k,2]+win_slide > len(mono_pca_NEW[0]):
                    print('Weird error, look closer at your data')
                    continue
                for m in AA_num_key:
                    if mono_pca_NEW[j,int(max_diffs[k,2]+win_slide)]==m:
                        # So I think that 0 and 1 should correspond to charge and hydrophobicity...
                        new_mat_mono[k,j]=new_mat_mono[k,j] + props[int(max_diffs[k,1]),m-1]

    return(new_mat_mono/win_size)

# CAN WE DO IT WITH ONE MATRIX???
def get_bigass_matrix(ALL_mono, OneChain = False, giveSize=[],manuscript_arrange=False):
    # THIS TERM IS ACTUALLY IRRELEVANT. NEED TO DEFINE A "GET MASK" function
    
    if OneChain:
        mono_PCA = aims.gen_1Chain_matrix(ALL_mono,key=AA_num_key_new/np.linalg.norm(AA_num_key_new),binary=False,giveSize=giveSize)
        mono_MI = aims.gen_1Chain_matrix(ALL_mono,key=AA_num_key,binary=False,giveSize=giveSize)
    else:
        mono_PCA = aims.gen_tcr_matrix(ALL_mono,key=AA_num_key_new/np.linalg.norm(AA_num_key_new),binary=False,giveSize=giveSize,manuscript_arrange = manuscript_arrange)
        mono_MI = aims.gen_tcr_matrix(ALL_mono,key=AA_num_key,binary=False,giveSize=giveSize,manuscript_arrange=manuscript_arrange)

    BIG_mono = aims.getBig(mono_MI)
    amono,bmono,cmono = np.shape(BIG_mono)

    #SO WE CANT JUST USE NP.RESHAPE
    #BECAUSE THAT COMMAND IS AGNOSTIC TO THE 
    #FACT THAT OUR DATA IS ARRANGED BY CLONE...
    BIG_mono_final = np.zeros((bmono,amono*cmono))
    for i in np.arange(bmono):
        BIG_mono_final[i] = BIG_mono[:,i,:].reshape(amono*cmono)

    mono_pca_stack = np.hstack([mono_PCA,BIG_mono_final])
    return(mono_pca_stack)

def do_classy_mda(ALL_mono, ALL_poly, matsize = 100, OneChain = False, 
                  xVal = 'kfold',ridCorr = False, feat_sel = 'none', classif = 'mda'):
        
    mono_dim=np.shape(ALL_mono)[1]
    poly_dim=np.shape(ALL_poly)[1]

    y_mono_all = np.ones(mono_dim)
    y_poly_all = 2*np.ones(poly_dim)

    y_all = np.hstack((y_mono_all,y_poly_all))
    seqs_all = np.hstack((ALL_mono,ALL_poly))
    
    # Stupid to recreate this matrix every single time... Should do it BEFORE, then splitting
    bigass_matrix = get_bigass_matrix(seqs_all, OneChain = OneChain)
        
    # Alright so here is where the data is actually split into the test/train/etc.
    if xVal == 'loo':
        crosser = LeaveOneOut()
        
    elif xVal == 'kfold':
        crosser = KFold(n_splits=10,shuffle=True) # Lets do 10x cross validation
        
    elif xVal == 'strat_kfold':
        # Random State can be defined if we want a reproducible shuffle
        crosser = StratifiedKFold(n_splits=10, random_state=None, shuffle=True)

    cycle = 0
    for train, test in crosser.split(bigass_matrix, y_all):

        X_train_pre, X_test_pre, y_train, y_test = bigass_matrix[train,:], bigass_matrix[test,:], y_all[train], y_all[test]

        # REMOVE HIGHLY CORRELATED FEATURES
        if ridCorr:
            pandaMatTrain = pandas.DataFrame(X_train_pre)
            pandaMatTest = pandas.DataFrame(X_test_pre)
            drop_zeros = [column for column in pandaMatTrain.columns if all(pandaMatTrain[column] == 0 )]
            NOzeroTrain = pandaMatTrain.drop(pandaMatTrain[drop_zeros], axis=1)
            NOzeroTest = pandaMatTest.drop(pandaMatTest[drop_zeros], axis=1)
            corrScoreTrain = NOzeroTrain.corr().abs()
            # Select upper triangle of correlation matrix
            upperTrain = corrScoreTrain.where(np.triu(np.ones(corrScoreTrain.shape), k=1).astype(np.bool))
            to_drop = [column for column in upperTrain.columns if ( any(upperTrain[column] > 0.75) ) ]
            finalTrain = NOzeroTrain.drop(NOzeroTrain[to_drop], axis=1)
            finalTest = NOzeroTest.drop(NOzeroTest[to_drop], axis=1)
            X_train = np.array(finalTrain)#; cols = final.columns
            X_test = np.array(finalTest)
        else:
            X_train = X_train_pre
            X_test = X_test_pre
        # !!!!!!!!!!! FEATURE SELECTION STEP !!!!!!!!!!!!!
        if feat_sel == 'PCA':
            pca = PCA(n_components=matsize, svd_solver='full')
            train_mat=pca.fit_transform(X_train)
        
        elif feat_sel == 'kPCA':   # What about Kernel PCA?
            pca = KernelPCA(n_components=matsize)
            train_mat=pca.fit_transform(X_train)
        
        elif feat_sel == 'kbest':   # Sklearn's "Kbest" module
            pca = SelectKBest(mutual_info_classif,k=matsize)
            train_mat = pca.fit_transform(X_train,y_train)
        
        elif feat_sel == 'max_diff':
            # Go back to my way of picking out the best features?
            # Need to split the poly and mono into separate matrices for this one...
            max_diffs = aims.parse_props(X_train,y_train,mat_size=matsize)
            indices = [int(a) for a in max_diffs[:,1]]
            train_mat = X_train[:,indices]
        
        elif feat_sel == 'none':
            test_mat = X_test
            train_mat = X_train
        
        # Apply the feature selection to the training dataset...
        if feat_sel == 'max_diff':
            test_mat = X_test[:,indices]
        elif feat_sel == 'none':
            test_mat = X_test
        else:
            test_mat = pca.transform(X_test)
        
        # Apply the actual classifier to the data
        if classif == 'mda':  # MDA
            clf_all = LinearDiscriminantAnalysis(n_components=1,solver='svd')
            clf_all.fit(train_mat,y_train)
        
        elif classif == 'svm': # SVM
            clf_all = SVC(gamma='auto',kernel='linear')
            clf_all.fit(train_mat, y_train)
        
        elif classif == 'logReg':  #Logistic Regression
            clf_all = LogisticRegression(solver = 'lbfgs')
            clf_all.fit(train_mat, y_train)
            
        elif classif == 'forest':   # Random Forest
            num_TREES = 500
            clf_all=RandomForestClassifier(n_estimators=num_TREES)
            clf_all.fit(train_mat,y_train)

        p_all = clf_all.predict(test_mat)
        acc_all = accuracy_score(y_test.flatten(),p_all)
        #fpr, tpr, thresholds = auc(y_test.flatten(), pAll_proba[:,0],pos_label=1)
        if cycle == 0:
            acc_fin = acc_all
            cycle = cycle + 1
        else:
            acc_fin = np.vstack((acc_fin,acc_all))
            cycle = cycle + 1
        #print(acc_all)
    return(acc_fin)

def apply_pretrained_LDA(bigass_mono,did_drop,indices,weights):
    # Take already bigass matrices and drop entries to look indentical to 
    # Need to have all these pre-defined variables in there
    prop_list_old = ['Phobic1','Charge','Phobic2','Bulk','Flex','Kid1','Kid2','Kid3','Kid4',
    'Kid5','Kid6','Kid7','Kid8','Kid9','Kid10']
    prop_list_new = ['Hot'+str(b+1) for b in range(46)]

    prop_names = prop_list_old + prop_list_new
    num_locs = int(np.shape(bigass_mono)[1]/61)
    Bigass_names = []
    for i in prop_names:
        for j in np.arange(num_locs):
            Bigass_names = Bigass_names + [ i + '-' + str(j) ]

    x = pandas.DataFrame(bigass_mono,columns = Bigass_names)

    dropped = x.drop(x[did_drop], axis=1)
    uncorr_mat = np.array(dropped)

    # The original indices were 
    pre_final = uncorr_mat[:,indices]
    final_apply=np.matmul(pre_final,np.transpose(weights))

    return(final_apply)

# DISTINCT FROM CLASSIFICATION. HERE IS HOW WE FIND THE PROPERTIES WHICH
# BEST DISCRIMINATE THE DATASET
def do_linear_split(test_mono,test_poly,ridCorr = True,giveSize=[],matSize=75,manuscript_arrange=False,pca_split=False):
    num_mono = np.shape(test_mono)[1]
    num_poly = np.shape(test_poly)[1]

    mat = np.hstack((test_mono,test_poly))
    if manuscript_arrange:
        total_mat = get_bigass_matrix(mat,giveSize = giveSize,manuscript_arrange=True)
    else:
        total_mat = get_bigass_matrix(mat,giveSize = giveSize,manuscript_arrange=False)
    prop_list_old = ['Phobic1','Charge','Phobic2','Bulk','Flex','Kid1','Kid2','Kid3','Kid4','Kid5','Kid6','Kid7','Kid8','Kid9','Kid10']
    prop_list_new = ['Hot'+str(b+1) for b in range(46)]

    prop_names = prop_list_old + prop_list_new
    num_locs = int(np.shape(total_mat)[1]/61)
    Bigass_names = []
    for i in prop_names:
        for j in np.arange(num_locs):
            Bigass_names = Bigass_names + [ i + '-' + str(j) ]
            
    if ridCorr:
        x = pandas.DataFrame(total_mat,columns = Bigass_names)
        drop_zeros = [column for column in x.columns if all(x[column] == 0 )]
        y = x.drop(x[drop_zeros], axis=1)
        z = y.corr().abs()
        # Select upper triangle of correlation matrix
        upper = z.where(np.triu(np.ones(z.shape), k=1).astype(np.bool))
        to_drop = [column for column in upper.columns if ( any(upper[column] > 0.75) ) ]
        final = y.drop(y[to_drop], axis=1)
        X_train = np.array(final); cols = final.columns

        did_drop = to_drop + drop_zeros
    else:
        X_train = total_mat; cols = np.array(Bigass_names)

    Y_train = np.hstack((np.ones(num_mono),2*np.ones(num_poly)))
    
    if pca_split:
        pca = PCA(n_components=matSize, svd_solver='full')
        train_mat=pca.fit_transform(X_train)
    else:
        Class_mat = aims.parse_props(np.transpose(X_train),Y_train,matSize)
        indices = [int(a) for a in Class_mat[:,1]]

        train_mat = X_train[:,indices]

    clf_all = LinearDiscriminantAnalysis(n_components=4,solver='svd')    
    mda_all=clf_all.fit_transform(train_mat,Y_train)

    # NOTE, THIS IS LIKE A "BEST CASE-SCENARIO" Accuracy
    # BUT, this does tell you how to best discriminate two classes.
    p_all = clf_all.predict(train_mat)
    acc_all = accuracy_score(Y_train.flatten(),p_all)
    # Give me the coefficients
    weights=clf_all.coef_
    
    if ridCorr and pca_split == False:
        return(acc_all,weights,cols,indices,mda_all,did_drop)
    elif pca_split:
        return(acc_all,mda_all)
    else:
        return(acc_all,weights,cols,indices,mda_all)

######### SAME AS DO_CLASSY_MDA, BUT NOW SPECIFY TEST-TRAIN ##########################
# Inputs should be complete matrices, not split into poly vs non-poly
def classy_apply(test_mat,y_test,train_mat,y_train, matsize = 100, OneChain = False, 
                  ridCorr = False, feat_sel = 'none', classif = 'mda'):
        
    max_len1 = aims.get_sequence_dimension(test_mat)[0]
    max_len2 = aims.get_sequence_dimension(train_mat)[0]
    max_len = np.zeros(6)
    for i in np.arange(6):
        max_len[i]=max(max_len1[i],max_len2[i])
    
    # Get the massive property matrices. These names are set to match
    # the names in the do_classy_mda code
    X_test_pre = get_bigass_matrix(test_mat, OneChain = OneChain,giveSize=max_len)
    X_train_pre = get_bigass_matrix(train_mat, OneChain = OneChain,giveSize=max_len)

    # REMOVE HIGHLY CORRELATED FEATURES
    if ridCorr:
        pandaMatTrain = pandas.DataFrame(X_train_pre)
        pandaMatTest = pandas.DataFrame(X_test_pre)
        drop_zeros = [column for column in pandaMatTrain.columns if all(pandaMatTrain[column] == 0 )]
        NOzeroTrain = pandaMatTrain.drop(pandaMatTrain[drop_zeros], axis=1)
        NOzeroTest = pandaMatTest.drop(pandaMatTest[drop_zeros], axis=1)
        corrScoreTrain = NOzeroTrain.corr().abs()
        # Select upper triangle of correlation matrix
        upperTrain = corrScoreTrain.where(np.triu(np.ones(corrScoreTrain.shape), k=1).astype(np.bool))
        to_drop = [column for column in upperTrain.columns if ( any(upperTrain[column] > 0.75) ) ]
        finalTrain = NOzeroTrain.drop(NOzeroTrain[to_drop], axis=1)
        finalTest = NOzeroTest.drop(NOzeroTest[to_drop], axis=1)
        X_train = np.array(finalTrain)#; cols = final.columns
        X_test = np.array(finalTest)
    else:
        X_train = X_train_pre
        X_test = X_test_pre
    # !!!!!!!!!!! FEATURE SELECTION STEP !!!!!!!!!!!!!
    if feat_sel == 'PCA':
        pca = PCA(n_components=matsize, svd_solver='full')
        train_mat=pca.fit_transform(X_train)
    
    elif feat_sel == 'kPCA':   # What about Kernel PCA?
        pca = KernelPCA(n_components=matsize)
        train_mat=pca.fit_transform(X_train)
    
    elif feat_sel == 'kbest':   # Sklearn's "Kbest" module
        pca = SelectKBest(mutual_info_classif,k=matsize)
        train_mat = pca.fit_transform(X_train,y_train)
    
    elif feat_sel == 'max_diff':
        # Go back to my way of picking out the best features?
        # Need to split the poly and mono into separate matrices for this one...
        max_diffs = aims.parse_props(X_train,y_train,mat_size=matsize)
        indices = [int(a) for a in max_diffs[:,1]]
        train_mat = X_train[:,indices]
    
    elif feat_sel == 'none':
        test_mat = X_test
        train_mat = X_train
    
    # Apply the feature selection to the training dataset...
    if feat_sel == 'max_diff':
        test_mat = X_test[:,indices]
    elif feat_sel == 'none':
        test_mat = X_test
    else:
        test_mat = pca.transform(X_test)
    
    # Apply the actual classifier to the data
    if classif == 'mda':  # MDA
        clf_all = LinearDiscriminantAnalysis(n_components=1,solver='svd')
        clf_all.fit(train_mat,y_train)
    
    elif classif == 'svm': # SVM
        clf_all = SVC(gamma='auto',kernel='linear')
        clf_all.fit(train_mat, y_train)
    
    elif classif == 'logReg':  #Logistic Regression
        clf_all = LogisticRegression(solver = 'lbfgs')
        clf_all.fit(train_mat, y_train)
        
    elif classif == 'forest':   # Random Forest
        num_TREES = 500
        clf_all=RandomForestClassifier(n_estimators=num_TREES)
        clf_all.fit(train_mat,y_train)

    p_all = clf_all.predict(test_mat)
    acc_all = accuracy_score(y_test.flatten(),p_all)
    #fpr, tpr, thresholds = auc(y_test.flatten(), pAll_proba[:,0],pos_label=1)
    return(acc_all)
