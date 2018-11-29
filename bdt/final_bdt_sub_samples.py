import os, sys
sys.path.insert(0, os.environ['HOME'] + '/.local/lib/python2.7/site-packages')

import numpy as np
import matplotlib as mpl
mpl.use('Agg')
# mpl.use('TkAgg')
from scipy.stats import ks_2samp
import matplotlib.pyplot as plt
import pandas, root_numpy
from copy import deepcopy
import pickle
import sklearn
from   sklearn.externals import joblib
from   sklearn.ensemble  import GradientBoostingClassifier
from   sklearn.metrics   import roc_curve
from   sklearn.model_selection import train_test_split
from   array       import array

from   xgboost import XGBClassifier, plot_importance
from   ROOT    import TFile, TTree, TH1F, gROOT, TChain


mc_sigma = 0.0400 
mc_mass  = 5.27783 

gROOT.SetBatch(True)

## weight input mass
min_m  = 4.8
max_m  = 5.8
n_bins = 500


###############################################
#####   SIGNAL AND BACKGROUND SELECTION   #####
###############################################
sig_selection_cutbased = '((tagB0==1 && (bMass    > {M} - 2.5*{S} && bMass    < {M} + 2.5*{S})) || \
                           (tagB0==0 && (bBarMass > {M} - 2.5*{S} && bBarMass < {M} + 2.5*{S}))) &&\
                            truthMatchMum==1 && truthMatchMup ==1 && truthMatchTrkm==1 && truthMatchTrkp ==1 && \
                           (tagB0==1 && genSignal==1) | (tagB0==0 && genSignal==2) && \
                            trig==0 && \
                           (mumuMass < 2.702 || mumuMass > 4.) '.format( M=mc_mass,S=mc_sigma)


bkg_selection_cutbased = '((tagB0==1 && ((bMass    > {M}-6.*{S} && bMass    < {M}-3.*{S}) || (bMass    > {M}+3.*{S} && bMass    < {M}+6*{S}))) || \
                           (tagB0==0 && ((bBarMass > {M}-6.*{S} && bBarMass < {M}-3.*{S}) || (bBarMass > {M}+3.*{S} && bBarMass < {M}+6*{S}))) && \
                           (mumuMass < 2.702 || mumuMass > 4. )) '.format(M=mc_mass,S=mc_sigma)


#############################
#####   INPUT SAMPLES   #####
#############################

sig_list = [
          'sub_samples/sample_mc_LMNR_events_0.root',
          'sub_samples/sample_mc_LMNR_events_1.root',
          'sub_samples/sample_mc_LMNR_events_2.root',
          'sub_samples/sample_mc_LMNR_events_3.root',
          'sub_samples/sample_mc_LMNR_events_4.root',
          'sub_samples/sample_mc_LMNR_events_5.root',
          'sub_samples/sample_mc_LMNR_events_6.root',
          'sub_samples/sample_mc_LMNR_events_7.root',
          'sub_samples/sample_mc_LMNR_events_8.root',
          'sub_samples/sample_mc_LMNR_events_9.root',
          'sub_samples/sample_mc_LMNR_events_10.root',
]

bkg_list = [
          'sub_samples/sample_data_LMNR_events_0.root',
          'sub_samples/sample_data_LMNR_events_1.root',
          'sub_samples/sample_data_LMNR_events_2.root',
          'sub_samples/sample_data_LMNR_events_3.root',
          'sub_samples/sample_data_LMNR_events_4.root',
          'sub_samples/sample_data_LMNR_events_5.root',
          'sub_samples/sample_data_LMNR_events_6.root',
          'sub_samples/sample_data_LMNR_events_7.root',
          'sub_samples/sample_data_LMNR_events_8.root',
          'sub_samples/sample_data_LMNR_events_9.root',
          'sub_samples/sample_data_LMNR_events_10.root',
]


for isample in range(11):

    tag = '_forPlots_%s'%isample

#####################################
#####   FEATURES AND BRANCHES   #####
#####################################
    
    features = [
        'bCosAlphaBS',
        'bLBS/bLBSE',
        'kstTrkmDCABS/kstTrkmDCABSE',
        'kstTrkpDCABS/kstTrkpDCABSE',
        'bVtxCL',
        'bDCABS/bDCABSE',
        'kstTrkmMinIP2D', ## min impact parameter of the track from any PV (in 2D)
        'kstTrkpMinIP2D', ## min impact parameter of the track from any PV (in 2D)
    ]
    
    branches = features + [
        'bMass',
        'bBarMass',
        'tagB0',
        'mumNTrkLayers',
        'mupNTrkLayers',
        'mumNPixLayers',
        'mupNPixLayers',
        'mupdxyVtx',
        'mumdxyVtx',
        'mumdzVtx',
        'mupdzVtx',
        'mupHighPurity',
        'mumHighPurity',
        'mumTMOneStationTight',
        'mupTMOneStationTight',
        'mupPt',
        'mumPt',
        'kstTrkpPt',
        'kstTrkmPt',
        'kstTrkpHighPurity',
        'kstTrkmHighPurity',
        'kstTrkmDCABSE',
        'kstTrkmDCABS',
        'kstTrkpDCABS',
        'kstTrkpDCABSE',
        'kstMass',
        'kstBarMass',
        'kkMass',
        'bLBS',
        'bLBSE',    
        'kstTrkmGlobalMuon',
        'kstTrkmNTrkLayers',
        'kstTrkmNPixHits',
        'kstTrkpGlobalMuon',
        'kstTrkpNTrkLayers',
        'kstTrkpNPixHits',
        'mumuMass',
        'mumIsoPt_dr04',
        'mupIsoPt_dr04',
        'kstTrkmIsoPt_dr04',
        'kstTrkpIsoPt_dr04',
    ]

    branches = list(set(branches))
    
    presel = '( mumNTrkLayers >= 6)         && ( mupNTrkLayers >= 6 ) && \
              ( mumNPixLayers >= 1)         && ( mupNPixLayers >= 1 ) && \
              ( mumdxyVtx < 0.3)            && ( mupdxyVtx < 0.3    ) && \
              ( mumdzVtx < 20 )             && ( mupdzVtx  < 20     ) && \
              ( mumHighPurity == 1 )        && ( mupHighPurity == 1 ) && \
              ( mumTMOneStationTight == 1 ) &&  mupTMOneStationTight == 1  && \
              ( kstTrkmHighPurity == 1 )    &&  kstTrkpHighPurity == 1     && \
              ( kkMass > 1.035 ) && \
              (!((kstTrkmGlobalMuon == 1) && ( kstTrkmNTrkLayers > 5 ) && (kstTrkmNPixHits > 0))) && \
              (!((kstTrkpGlobalMuon == 1) && ( kstTrkpNTrkLayers > 5 ) && (kstTrkpNPixHits > 0)))'
              

    #################################
    #####   WEIGHT INPUT MASS   #####
    #################################

    hmass_bkg = TH1F('hmass_bkg', 'hmass_bkg', n_bins, min_m, max_m)
    t = TChain('ntuple')
    for i, ifile in enumerate(bkg_list):
        if i == isample:  
            continue 
        t.Add(ifile)
    t.Draw('bMass*tagB0 + bBarMass*(1-tagB0) >> hmass_bkg', '(mumuMass < 2.702 || mumuMass > 4.) && %s '%presel)
    clone_hm_bkg = deepcopy(hmass_bkg)
    
    weights_bkg = {}
    for i in range(hmass_bkg.GetNbinsX()):
        weights_bkg[hmass_bkg.GetBinCenter(i+1)] = hmass_bkg.Integral(hmass_bkg.FindBin(min_m),hmass_bkg.FindBin(max_m))/hmass_bkg.GetBinContent(i+1)/n_bins

    @np.vectorize
    def add_m_weights_bkg(theMass):
        m_weight  = weights_bkg[clone_hm_bkg.GetBinCenter(clone_hm_bkg.FindBin(theMass))]
        return m_weight
    
    t = TChain('ntuple')
    for i, ifile in enumerate(sig_list):
        if i == isample:  
            continue 
        t.Add(ifile)
    
    hmass_sig = TH1F('hmass_sig', 'hmass_sig', n_bins, min_m, max_m)
    t.Draw('bMass*tagB0 + bBarMass*(1-tagB0) >> hmass_sig', 'weight*((mumuMass < 2.702 || mumuMass > 4.) && trig == 0 && ((tagB0==1 && genSignal==1) || (tagB0==0 && genSignal==2)) && %s )'%presel)
    clone_hm_sig = deepcopy(hmass_sig)
    
    weights_sig = {}
    for i in range(hmass_sig.GetNbinsX()):
        weights_sig[hmass_sig.GetBinCenter(i+1)] = hmass_sig.Integral(hmass_sig.FindBin(min_m),hmass_sig.FindBin(max_m))/hmass_sig.GetBinContent(i+1)/n_bins
    
    @np.vectorize
    def add_m_weights_sig(theMass):
        m_weight  = weights_sig[clone_hm_sig.GetBinCenter(clone_hm_sig.FindBin(theMass))]
        return m_weight


    #################################
    #####   LOAD DATASETS       #####
    #################################
    sig = pandas.DataFrame(
        root_numpy.root2array(
            sig_list[:isample] + sig_list[(isample + 1):], 
            'ntuple',
            branches  = branches + ['weight'],
            selection = sig_selection_cutbased,
        )
    )
    
    bkg = pandas.DataFrame(
        root_numpy.root2array(
            bkg_list[:isample] + bkg_list[(isample + 1):], 
            'ntuple',
            branches  = branches,
            selection = bkg_selection_cutbased,
        )
    )

    ##################################
    #####   DEFINE THE TARGETS   #####
    ##################################
    
    sig['target'] = np.ones (sig.shape[0]).astype(np.int)
    bkg['target'] = np.zeros(bkg.shape[0]).astype(np.int)
    
    
    ## define isolation: # tracks with pt in a cone
    sig['isopt_mum_04']  = sig.mumIsoPt_dr04/sig.mumPt
    sig['isopt_mup_04']  = sig.mupIsoPt_dr04/sig.mupPt
    sig['isopt_trkm_04'] = sig.kstTrkmIsoPt_dr04/sig.kstTrkmPt
    sig['isopt_trkp_04'] = sig.kstTrkpIsoPt_dr04/sig.kstTrkpPt
    sig['sum_isopt_04']  = sig.isopt_mum_04 + sig.isopt_mup_04 + sig.isopt_trkm_04 + sig.isopt_trkp_04
    
    bkg['isopt_mum_04']  = bkg.mumIsoPt_dr04/bkg.mumPt
    bkg['isopt_mup_04']  = bkg.mupIsoPt_dr04/bkg.mupPt
    bkg['isopt_trkm_04'] = bkg.kstTrkmIsoPt_dr04/bkg.kstTrkmPt
    bkg['isopt_trkp_04'] = bkg.kstTrkpIsoPt_dr04/bkg.kstTrkpPt
    bkg['sum_isopt_04']  = bkg.isopt_mum_04 + bkg.isopt_mup_04 + bkg.isopt_trkm_04 + bkg.isopt_trkp_04
    
    features.append('sum_isopt_04')
    
    
    ## use only one combination for K* mass
    sig['kstarmass']  = sig.tagB0*sig.kstMass +(1- sig.tagB0)*sig.kstBarMass
    bkg['kstarmass']  = bkg.tagB0*bkg.kstMass +(1- bkg.tagB0)*bkg.kstBarMass
    features.append('kstarmass')
    
    ## add b mass calculation for plotting only
    bkg['themass']  = bkg.bMass*bkg.tagB0 + bkg.bBarMass*(1-bkg.tagB0)
    sig['themass']  = sig.bMass*sig.tagB0 + sig.bBarMass*(1-sig.tagB0)
    branches.append('themass')
    
    bkg.hist(column='themass',bins=400)
    plt.savefig('mass_bkg.pdf')
    sig.hist(column='themass',bins=400)   
    plt.savefig('mass_sig.pdf')
    
    bkg['m_weight']  = add_m_weights_bkg(bkg.themass)
    sig['m_weight']  = add_m_weights_sig(sig.themass)

    ## add normfactor here
    sig['normfactor'] = sig.weight*sig.m_weight #/len(sig) *10000
    bkg['normfactor'] = bkg.m_weight# /len(bkg) #* 10000
    
    
    ###########################################################
    #####   PRESELCTION & SPLIT TRAIN AND TEST SAMPLES    #####
    ###########################################################
    data_all = pandas.concat([sig, bkg])
    
    ## add column for pass-preselection
    data_all['pass_preselection'] = ( data_all.mumNTrkLayers >= 6)  & ( data_all.mupNTrkLayers >= 6 ) & \
                                    ( data_all.mumNPixLayers >= 1)  & ( data_all.mupNPixLayers >= 1 ) & \
                                    ( data_all.mumdxyVtx < 0.3)     & ( data_all.mupdxyVtx < 0.3    ) & \
                                    ( data_all.mumdzVtx < 20 )      & ( data_all.mupdzVtx  < 20     ) & \
                                    ( data_all.mumHighPurity == 1 ) & ( data_all.mupHighPurity == 1 ) & \
                                    ( data_all.mumTMOneStationTight == 1 ) & ( data_all.mupTMOneStationTight == 1 ) & \
                                    ( data_all.kstTrkmHighPurity == 1 )    & ( data_all.kstTrkpHighPurity == 1    ) & \
                                    ( data_all.kkMass > 1.035 ) & \
                                    (~((data_all.kstTrkmGlobalMuon == 1) & ( data_all.kstTrkmNTrkLayers > 5 ) & ( data_all.kstTrkmNPixHits > 0))) & \
                                    (~((data_all.kstTrkpGlobalMuon == 1) & ( data_all.kstTrkpNTrkLayers > 5 ) & ( data_all.kstTrkpNPixHits > 0))) 
    
    '''
     split the sample to have:
     test_unbias = X% of original data (no preselection), not used for training nor validating
     data = 100-X% of original data, also passing preselection
    '''
    
    data = data_all[data_all.pass_preselection==1]
    train, test = train_test_split(data, test_size=0.30, random_state = 17)
    

    ###########################################################
    #####   CLASSIFIER DECLARATION & TRAINING             #####
    ###########################################################
    clf = XGBClassifier(
        max_depth        = 4,
        learning_rate    = 0.05,
        n_estimators     = 700,
        subsample        = 0.5,
        colsample_bytree = 0.6,
        min_child_weight = 4.,
        gamma            = 4,     ## optimized
        reg_alpha        = 0,
        reg_lambda       = 1,
        seed             = 1986,
        silent           = False,
        )
        
    clf.fit(
        train[features], 
        train.target,
        eval_set              = [(train[features], train.target), (test[features], test.target)],
        early_stopping_rounds = 50,
        eval_metric           = 'auc',
        verbose               = True,
        sample_weight         = train['normfactor'],
    )
    
    joblib.dump(clf, 'results/classifier_%s.pkl' %(tag), compress=True)
    
    ####################################################
    # #####   PREDICT ON THE TEST UNBIAS SAMPLE    #####
    # ##################################################
    pred = clf.predict_proba(test[features])[:, 1] #### !!! changed from test_unbias
    
    train_sig = clf.predict_proba(train[features][train.target>0.5])[:,1]
    train_bkg = clf.predict_proba(train[features][train.target<0.5])[:,1]
    
    test_sig = clf.predict_proba(test[features][test.target>0.5])[:,1]
    test_bkg = clf.predict_proba(test[features][test.target<0.5])[:,1]
    
    
    ###########################
    #####   ROC CURVE    ######
    ###########################
    plt.clf()
    
    import itertools
    xy = [i*j for i,j in itertools.product([10.**i for i in range(-8, 0)], [1,2,4,8])]+[1]
    plt.plot(xy, xy, color='grey', linestyle='--')
    plt.xlim([10**-3, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    
    # ## draw ROC         
    fpr, tpr, threshold = roc_curve(test.target, (pred * test.pass_preselection) ) ### !! changed from test unbias
    plt.plot(fpr, tpr, label='BDT', color='r', label='ROC test')
    
    pred2 = clf.predict_proba(train[features])[:, 1]
    fpr2, tpr2, sara2 = roc_curve(train.target, pred2)
    plt.plot(fpr2, tpr2, label='BDT', color='b', label='ROC train')

    plt.xscale('log')
    plt.grid()
    
    roc_file = open('roc_%s.pck' %(tag), 'w+')
    pickle.dump((tpr, fpr), roc_file)
    roc_file.close()
    
    plt.legend(loc='best')
    plt.grid()
    plt.title('ROC')
    plt.tight_layout()
    plt.savefig('results/roc_train_test_%s.pdf' %(tag))
    plt.clf()

    ### add plot of BDT vs tpr     
    plt.plot(threshold, tpr, color='b')
    
    plt.legend(loc='best')
    plt.grid()
    plt.title('turn on')
    plt.tight_layout()
    plt.xlabel('BDT output')
    plt.ylabel('True Positive Rate')
    plt.savefig('results/turn_on_%s.pdf' %(tag))
    plt.clf()
    
    roc_file = open('results/turnon_%s.pck' %(tag), 'w+')
    pickle.dump((threshold, tpr), roc_file)
    roc_file.close()
    
    ### add plot of BDT vs fpr  
    plt.plot(threshold, fpr, color='b')
    plt.legend(loc='best')
    plt.grid()
    plt.title('rate')
    plt.xlabel('BDT output')
    plt.ylabel('False Positive Rate')
    plt.yscale('log')
    
    plt.tight_layout()
    plt.savefig('results/rate_%s.pdf' %(tag))
    plt.clf()
    
    roc_file = open('results/rate_%s.pck' %(tag), 'w+')
    pickle.dump((threshold, fpr), roc_file)
    roc_file.close()
    
    ###################################
    #####   OVERTRAINING TEST    ######
    ###################################
    
    low  = 0
    high = 1
    low_high = (low,high)
    bins = 50
    
    plt.hist(
        train_sig,
        color='r', 
        alpha=0.5, 
        range=low_high, 
        bins=bins,
        histtype='stepfilled', 
        normed=True,
        log=True,
        label='S (train)'
    )
    
    #################################################
    plt.hist(
        train_bkg,
        color='b', 
        alpha=0.5, 
        range=low_high, 
        bins=bins,
        histtype='stepfilled', 
        normed=True,
        log=True,
        label='B (train)'
    )
    
    #################################################
    hist, bins = np.histogram(
        test_sig,
        bins=bins, 
        range=low_high, 
        normed=True,
    )
    
    width  = (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    scale  = len(test_sig) / sum(hist)
    err    = np.sqrt(hist * scale) / scale
    
    plt.errorbar(
        center, 
        hist, 
        yerr=err, 
        fmt='o', 
        c='r', 
        label='S (test)'
    )
    
    #################################################
    hist, bins = np.histogram(
        test_bkg,
        bins=bins, 
        range=low_high, 
        normed=True
    )
    
    width  = (bins[1] - bins[0])
    center = (bins[:-1] + bins[1:]) / 2
    scale  = len(test_bkg) / sum(hist)
    err    = np.sqrt(hist * scale) / scale
    
    plt.errorbar(
        center, 
        hist, 
        yerr=err, 
        fmt='o', 
        c='b', 
        label='B (test)'
    )
    
    #################################################
    plt.xlabel('BDT output')
    plt.ylabel('Arbitrary units')
    plt.legend(loc='best')
    ks_sig = ks_2samp(train_sig, test_sig)
    ks_bkg = ks_2samp(train_bkg, test_bkg)
    plt.suptitle('KS p-value: sig = %.3f%s - bkg = %.2f%s' %(ks_sig.pvalue * 100., '%', ks_bkg.pvalue * 100., '%'))
    plt.savefig('results/overtrain_%s.pdf' %(tag))
    
    ###################################
    #####   FEATURE IMPORTANCE    #####
    ###################################
    plot_importance(clf)
    plt.tight_layout()
    plt.savefig('results/feat_importance_%s.pdf' %(tag))
    
    
    ###################################
    #####   OVERTRAINING SCORE    #####
    ###################################
    plt.clf()
    
    auc_train = clf.evals_result()['validation_0']['auc']
    auc_test  = clf.evals_result()['validation_1']['auc']
    
    n_estimators = np.arange(len(auc_train))
    
    plt.plot(n_estimators, auc_train, color='r', label='AUC train')
    plt.plot(n_estimators, auc_test , color='b', label='AUC test' )
    
    plt.xlabel('# tree')
    plt.ylabel('Area Under ROC')
    
    plt.xscale('log')
    plt.grid()
    
    # plt.xlim([1, 1000])
    # plt.ylim([0.985, 1.0])
    
    plt.legend(loc='best')
    plt.tight_layout()
    plt.savefig('results/auc_score_%s.pdf' %(tag))
    
    
    
    ######################################
    #### Correlation Matrix Plot    ######
    ######################################
    sig  ['bdt'] = clf.predict_proba(sig[features])[:, 1] 
    bkg  ['bdt'] = clf.predict_proba(bkg[features])[:, 1] 

    sig_corr = sig[features+ ['bdt', 'themass']]
    bkg_corr = bkg[features+ ['bdt', 'themass']]
    
    sig_correlations = sig_corr.corr()
    bkg_correlations = bkg_corr.corr()
    
    fig = plt.figure()
    ticks = np.arange(0,len(features)+2,1)
    
    ax  = fig.add_subplot(121)
    cax = ax.matshow(sig_correlations, vmax=1., vmin=-1, cmap=plt.cm.get_cmap('Blues', 9))

    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.set_xticklabels(features+ ['bdt', 'themass'],fontsize=5)
    ax.set_yticklabels(features+ ['bdt', 'themass'],fontsize=5)
    plt.xticks(rotation=90)
    
    # plot correlation matrix
    ax2  = fig.add_subplot(122)
    cax2 = ax2.matshow(bkg_correlations, vmax=1., vmin=-1, cmap=plt.cm.get_cmap('Blues', 11))
    ax2.set_xticks(ticks)
    ax2.set_yticks(ticks)
    ax2.set_xticklabels(features+ ['bdt', 'themass'],fontsize=5)
    ax2.set_yticklabels(features+ ['bdt', 'themass'],fontsize=0)
    plt.xticks(rotation=90)
    
    cbar = fig.colorbar(cax, orientation="horizontal", pad=0.045)
    cbar.ax.tick_params(labelsize=5)
    
    plt.show()
    plt.savefig('results/bkg_correlation_%s.pdf' %(tag))
    
    plt.close()
    
