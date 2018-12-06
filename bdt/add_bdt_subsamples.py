import os, sys
sys.path.insert(0, os.environ['HOME'] + '/.local/lib/python2.7/site-packages')

import ROOT
import numpy as np
import xgboost as xgb
import pandas, root_numpy
from sklearn.externals import joblib
from sklearn import preprocessing
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt 


for i in range(0,11):  
  
    tag = 'Nov21_'+ str(i)
     
    ifile = 'sub_samples/sample_mc_PSIPRIME_events_' + str(i) + '.root'  
  
    print 'adding bdt score from classifier%s.pkl' %tag
    classifier = joblib.load('classifiers/classifier_%s.pkl' %tag)
  
    feat_names = [
        'bCosAlphaBS',
        'bLBS/bLBSE',
        'kstTrkmDCABS/kstTrkmDCABSE',
        'kstTrkpDCABS/kstTrkpDCABSE',
        'bVtxCL',
        'bDCABS/bDCABSE',
        'kstTrkmMinIP2D', ## min impact parameter of the track from any PV (in 2D)
        'kstTrkpMinIP2D', ## min impact parameter of the track from any PV (in 2D)
    ]

    additional = [
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
        'bPt',
        'kstPt',
        'kstTrkpTrackerMuonArbitrated',
        'kstTrkmTrackerMuonArbitrated',
        'kstTrkpHighPurity',
        'kstTrkmHighPurity',
        'kstTrkpPt',
        'kstTrkmPt',
        'kstTrkmDCABSE',
        'kstTrkmDCABS',
        'kstTrkpDCABS',
        'kstTrkpDCABSE',
        'kstBarMass',
        'kstMass',
        'kkMass',
        'bLBS',
        'bLBSE',    
        'kstTrkmIsoN',
        'kstTrkpIsoN',
        'kstTrkmGlobalMuon',
        'kstTrkmNTrkLayers',
        'kstTrkmNPixHits',
        'kstTrkpGlobalMuon',
        'kstTrkpNTrkLayers',
        'kstTrkpNPixHits',
        'eventN',
        'mumuMass',
        'mumIsoPt_dr04',
        'mupIsoPt_dr04',
        'kstTrkmIsoPt_dr04',
        'kstTrkpIsoPt_dr04',
    ]
    
    print 'loading support dataset...'
    dataset_support = pandas.DataFrame(
        root_numpy.root2array(
            ifile,
            'ntuple',
            branches=feat_names + additional,
        )
    )
    print '\t...done'
    
    
    print 'loading dataset...'
    dataset = pandas.DataFrame(
        root_numpy.root2array(
            ifile,
            'ntuple',
    #         branches = feat_names + additional + ['kstTrkmIsoN']
        )
    )
    print '\t...done'
    
    ## define isolation: # tracks with iso < 0.6 GeV
    ## define isolation: # tracks with pt in a cone
    dataset_support['isopt_mum_04']  = dataset_support.mumIsoPt_dr04    /dataset_support.mumPt
    dataset_support['isopt_mup_04']  = dataset_support.mupIsoPt_dr04    /dataset_support.mupPt
    dataset_support['isopt_trkm_04'] = dataset_support.kstTrkmIsoPt_dr04/dataset_support.kstTrkmPt
    dataset_support['isopt_trkp_04'] = dataset_support.kstTrkpIsoPt_dr04/dataset_support.kstTrkpPt
    dataset_support['sum_isopt_04']  = dataset_support.isopt_mum_04 + dataset_support.isopt_mup_04 + dataset_support.isopt_trkm_04 + dataset_support.isopt_trkp_04
    
    dataset['isopt_mum_04']  = dataset.mumIsoPt_dr04    /dataset.mumPt
    dataset['isopt_mup_04']  = dataset.mupIsoPt_dr04    /dataset.mupPt
    dataset['isopt_trkm_04'] = dataset.kstTrkmIsoPt_dr04/dataset.kstTrkmPt
    dataset['isopt_trkp_04'] = dataset.kstTrkpIsoPt_dr04/dataset.kstTrkpPt
    dataset['sum_isopt_04']  = dataset.isopt_mum_04 + dataset.isopt_mup_04 + dataset.isopt_trkm_04 + dataset.isopt_trkp_04
    
    feat_names.append('sum_isopt_04')
    
    additional.append('isopt_mum_04')
    additional.append('isopt_mup_04')
    additional.append('isopt_trkm_04')
    additional.append('isopt_trkp_04')
    
    dataset_support['kstarmass']  = dataset_support.tagB0*dataset_support.kstMass +(1- dataset_support.tagB0)*dataset_support.kstBarMass
    dataset['kstarmass']          = dataset.tagB0*dataset.kstMass +(1- dataset.tagB0)*dataset.kstBarMass
    feat_names.append('kstarmass')
    
    dataset['pass_preselection'] =  ( dataset.mumNTrkLayers >= 6)  & ( dataset.mupNTrkLayers >= 6 ) & \
                                    ( dataset.mumNPixLayers >= 1)  & ( dataset.mupNPixLayers >= 1 ) & \
                                    ( dataset.mumdxyVtx < 0.3)     & ( dataset.mupdxyVtx < 0.3    ) & \
                                    ( dataset.mumdzVtx < 20 )      & ( dataset.mupdzVtx  < 20     ) & \
                                    ( dataset.mumHighPurity == 1 ) & ( dataset.mupHighPurity == 1 ) & \
                                    ( dataset.mumTMOneStationTight == 1 ) & ( dataset.mupTMOneStationTight == 1 ) & \
                                    ( dataset.kstTrkmHighPurity == 1 )    & ( dataset.kstTrkpHighPurity == 1    ) & \
                                    ( dataset.kkMass > 1.035 ) & \
                                    (~((dataset.kstTrkmGlobalMuon == 1) & ( dataset.kstTrkmNTrkLayers > 5 ) & ( dataset.kstTrkmNPixHits > 0))) & \
                                    (~((dataset.kstTrkpGlobalMuon == 1) & ( dataset.kstTrkpNTrkLayers > 5 ) & ( dataset.kstTrkpNPixHits > 0))) 
    
    
    
    ## add branch for BDT 
    print 'computing probabilities...'
    bdt_prob_array = classifier.predict_proba(dataset_support[feat_names])[:,1]
    print '\t...done'
    
    print 'adding new column to the dataset...'
    dataset['bdt_prob'] = bdt_prob_array
    print '\t...done'	
    
    
    # https://github.com/scikit-hep/root_pandas
    import root_pandas
    dataset.to_root(ifile.replace('.root', '_addBDT.root'), key='ntuple', store_index=False)
    # ifile.GetName().replace('.root', '_%s.root'%tag)
