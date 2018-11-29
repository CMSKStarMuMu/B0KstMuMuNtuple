import ROOT, math, sys
import numpy as np
import pandas, root_numpy
from copy import deepcopy as dc

from ROOT import gSystem
gSystem.Load('libRooFit')
from ROOT import RooFit, RooRealVar, RooDataSet, RooArgList, RooTreeData, RooArgSet, RooAddPdf, RooFormulaVar, RooExtendPdf
from ROOT import RooGaussian, RooExponential, RooChebychev

ROOT.gROOT.SetBatch(True)
# ROOT.gErrorIgnoreLevel = ROOT.kFatal
# ROOT.RooMsgService.instance().setSilentMode(True)

B0Mass_   = 5.27958
JPsiMass_ = 3.096916
PsiPMass_ = 3.686109
KStMass_  = 0.896

B0Mass     = RooRealVar("B0Mass"    , "B0Mass"  , B0Mass_  )
JPsiMass   = RooRealVar("JPsiMass"  , "JPsiMass", JPsiMass_)
PsiPMass   = RooRealVar("PsiPMass"  , "PsiPMass", PsiPMass_)
KStMass    = RooRealVar("KStMass"   , "KStMass" , KStMass_ )

nSigma_psiRej = 3.

outfile   = ROOT.TFile('outcome_wp_finding_massFlatten_onlyRT_NoB0PsiCut.root', 'recreate')

n_bdt_points = 15
sample_range = 11

def fitNSig(ibdt, fulldata):
    mean        = RooRealVar ("mass"          , "mean"           ,  B0Mass_,   3,    7, "GeV")
    sigma       = RooRealVar ("#sigma_{1}"    , "sigma"          ,  0.028,     0,   10, "GeV")
    signalGauss = RooGaussian("signalGauss"   , "signal gauss"   ,  theBMass,  mean, sigma)
    
    sigma2       = RooRealVar ("#sigma_{2}"   , "sigma2"         ,  0.048,     0,   0.09, "GeV")
    signalGauss2 = RooGaussian("signalGauss2" , "signal gauss2"  ,  theBMass,  mean,sigma2)
    f1           = RooRealVar ("f1"           , "f1"             ,   0.8 ,     0.,   1.)
    gaus         = RooAddPdf  ("gaus"         , "gaus1+gaus2"    , RooArgList(signalGauss,signalGauss2), RooArgList(f1))
    
    pol_c1      = RooRealVar ("p1"            , "coeff x^0 term" ,    0.5,   -10, 10)
    pol_c2      = RooRealVar ("p2"            , "coeff x^1 term" ,    0.5,   -10, 10)
    # pol_c3      = RooRealVar ("p3"           , "coeff x^2 term",    0.5,   -10, 10)
    # slope       = RooRealVar ("slope"        , "slope"         ,    0.5,   -10, 10)
    # bkg_exp     = RooExponential("bkg_exp"   , "exponential"   ,  slope,   theBMass  )
    bkg_pol     = RooChebychev("bkg_pol"      , "2nd order pol"  ,  theBMass, RooArgList(pol_c1))
    
    nsig        = RooRealVar("Yield"          , "signal frac"    ,  40000,     0,   1000000)
    nbkg        = RooRealVar("nbkg"           , "bkg fraction"   ,   1000,     0,   550000)
    
    cut = cut_base + '&& bdt_prob > %s'%(ibdt)
    
    data        = fulldata.reduce(RooArgSet(theBMass,mumuMass,mumuMassE), cut)
    fitFunction = RooAddPdf ("fitfunction" , "fit function"  ,  RooArgList(gaus, bkg_pol), RooArgList(nsig, nbkg))
    r = fitFunction.fitTo(data, RooFit.Extended(True), RooFit.Save(), RooFit.Range(4.9,5.6))
    
    frame = theBMass.frame()
    data.plotOn(frame, RooFit.Binning(70), RooFit.MarkerSize(.7))
    fitFunction.plotOn(frame, )
    fitFunction.plotOn(frame, RooFit.Components("bkg_pol")     , RooFit.LineStyle(ROOT.kDashed))
    fitFunction.plotOn(frame, RooFit.Components("signalGauss") , RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kGreen+1))
    fitFunction.plotOn(frame, RooFit.Components("signalGauss2"), RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kOrange+1))
    
    parList = RooArgSet (nsig,sigma,sigma2, mean)
    # fitFunction.plotOn(frame, RooFit.Components("signalGauss2"), RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kGreen+2))
    
    fitFunction.paramOn(frame, RooFit.Parameters(parList), RooFit.Layout(0.62,0.86,0.88))
    frame.Draw()
    
    dict_s_v1[ibdt]  = [nsig.getVal(), nsig.getError()]
    dict_sigma[ibdt] = math.sqrt(f1.getVal()*(sigma.getVal()**2) + (1-f1.getVal())*(sigma2.getVal()**2))



def fitNBkg(ibdt, fullbkg):

    slope       = RooRealVar    ("slope"     , "slope"         ,    0.5,   -10, 10 )
    bkg_exp     = RooExponential("bkg_exp"   , "exponential"   ,  slope,   theBMass)
    
    cut = cut_base + '&& bdt_prob > %s'%(ibdt)
    
    theBMass.setRange('sigRangeMC', B0Mass_ - 3*dict_sigma[ibdt], B0Mass_ + 3*dict_sigma[ibdt])
    
    databkg = fullbkg.reduce(RooArgSet(theBMass,mumuMass,mumuMassE), cut)
    r       = bkg_exp.fitTo(databkg, RooFit.Save(), ROOT.RooFit.Range('left,right'))
    
    frame = theBMass.frame()
    databkg.plotOn(frame, RooFit.Binning(70), RooFit.MarkerSize(.7))
    bkg_exp.plotOn(frame, )
    
    # bkg_exp.fixCoefRange('left,right')
    
    nbkg = RooRealVar  ('nbkg', 'bkg n'  ,  1000,       0,   550000)
    ebkg = RooExtendPdf('ebkg','ebkg'    ,  bkg_exp, nbkg, 'sigRangeMC') 
    ebkg.fitTo(databkg, ROOT.RooFit.Range('left,right'))
    ebkg.plotOn(frame, RooFit.LineStyle(ROOT.kDashed), RooFit.LineColor(ROOT.kGreen+1), RooFit.Range(4.9,5.6))
    frame.Draw()
    
    dict_b_v1[ibdt] =  [nbkg.getVal(), nbkg.getError()]



for isample in range(sample_range):
    ifileMC   = '../sub_samples/sample_mc_LMNR_events_%s_addBDT_yesMassFlatten_yesSubSample_yesKstarMass_onlyRTevents.root'%isample
    ifileData = '../sub_samples/sample_data_LMNR_events_%s_addBDT_yesMassFlatten_yesSubSample_yesKstarMass_onlyRTevents.root'%isample

    dict_s_v1  = {}
    dict_b_v1  = {}
    dict_sigma = {}
        
    ### retrieve S from fitting the MC sample
    bMass     = RooRealVar("bMass"             , "#mu^{+}#mu^{-}K* mass", 2, 20, "GeV")
    bBarMass  = RooRealVar("bBarMass"          , "#mu^{+}#mu^{-}K* mass", 2, 20, "GeV")
    mumuMass  = RooRealVar("mumuMass"          , "mumuMass"             , 0, 1000)
    mumuMassE = RooRealVar("mumuMassE"         , "mumuMassE"            , 0, 10000)
    tagB0     = RooRealVar("tagB0"             , "tagB0"                , -1, 2)
    bdt_prob  = RooRealVar("bdt_prob"          , "bdt_prob"             ,  0.5, 2) ## already cut here on BDT !!!
    pass_pre  = RooRealVar("pass_preselection" , "pass_preselection"    ,  1  , 2) ## cut on preselection 
    
    thevars = RooArgSet()
    thevars.add(bMass)
    thevars.add(bBarMass)
    thevars.add(mumuMass)
    thevars.add(mumuMassE)
    thevars.add(tagB0)
    thevars.add(bdt_prob)

    tree = ROOT.TChain('ntuple')
    tree.AddFile(ifileMC)
    fulldata   = RooDataSet('fulldata', 'fulldataset', tree,  RooArgSet(thevars))
    
    ## add to the input tree the combination of the variables for the B0 arb. mass
    theBMassfunc = RooFormulaVar("theBMass", "#mu^{+}#mu^{-}K^{#pm}#pi^{#mp} mass [GeV]", "@0*@1 + (1-@0)*@2", RooArgList(tagB0,bMass,bBarMass) )
    ## add to the input tree the combination of the variables, to be used for the cuts on the dimuon mass
    theBMass     = fulldata.addColumn(theBMassfunc) 
    deltaB0Mfunc = RooFormulaVar("deltaB0M", "deltaB0M", "@0 - @1", RooArgList(theBMass,B0Mass) )
    deltaJMfunc  = RooFormulaVar("deltaJpsiM" , "deltaJpsiM" , "@0 - @1", RooArgList(mumuMass,JPsiMass) )
    deltaPMfunc  = RooFormulaVar("deltaPsiPM" , "deltaPsiPM" , "@0 - @1", RooArgList(mumuMass,PsiPMass) )
    
    theBMass.setRange(4.9,5.6)
    deltaB0M     = fulldata.addColumn(deltaB0Mfunc) 
    deltaJpsiM   = fulldata.addColumn(deltaJMfunc) 
    deltaPsiPM   = fulldata.addColumn(deltaPMfunc) 
    
#     cut_base = '( abs(mumuMass - {JPSIM}) > {CUT}*mumuMassE && abs(mumuMass - {PSIM}) > {CUT}*mumuMassE &&  \
#                 (( mumuMass < {JPSIM} && !( abs(deltaB0M - deltaJpsiM) < 0.16 || abs(deltaB0M - deltaPsiPM) < 0.06) ) || \
#                  ( mumuMass > {PSIM}  && !( abs(deltaB0M - deltaJpsiM) < 0.06 || abs(deltaB0M - deltaPsiPM) < 0.03) ) || \
#                  ( mumuMass > {JPSIM} && mumuMass < {PSIM} && !( abs(deltaB0M - deltaJpsiM) < 0.06 || abs(deltaB0M - deltaPsiPM) < 0.06 ))))'.format(JPSIM=JPsiMass_, PSIM=PsiPMass_,  CUT=nSigma_psiRej)  
    cut_base = '( mumuMass < 2.702 || mumuMass > 4. )'
    
    
    for i in range(n_bdt_points):
#         ibdt = i*0.002+0.96
        ibdt = i*0.01+0.85
        fitNSig(ibdt, fulldata)
    
    
    ### retrieve B from fitting the data sample
    treeBkg = ROOT.TChain('ntuple')
    treeBkg.AddFile(ifileData)
    fullbkg     = RooDataSet('fullbkg', 'fullbkg', treeBkg,  RooArgSet(thevars))
    theBMass    = fullbkg.addColumn(theBMassfunc) 
    deltaB0M    = fullbkg.addColumn(deltaB0Mfunc) 
    deltaJpsiM  = fullbkg.addColumn(deltaJMfunc)  
    deltaPsiPM  = fullbkg.addColumn(deltaPMfunc)  
    
    theBMass.setRange(4.9,5.6)
    theBMass.setRange('left'    , 4.9, 5.13)
    theBMass.setRange('right'   , 5.4, 5.6 )
    
    
    for i in range(n_bdt_points):
#         ibdt = i*0.002+0.96
        ibdt = i*0.01+0.85
        fitNBkg(ibdt, fullbkg)
    
    
    
    sign_vec_v1 = []
    bkg_vec_v1  = []
    bdt_vec     = []  
    for k,v in dict_s_v1.iteritems():
        bdt_vec.append(k)
        sign_vec_v1.append(v[0])
    for k,v in dict_b_v1.iteritems():
        bkg_vec_v1.append(v[0])
    
    sign_array_v1 = np.asarray(sign_vec_v1)
    bkg_array_v1  = np.asarray(bkg_vec_v1)
    bdt_array     = np.asarray(bdt_vec )
    

    graph_v1 = ROOT.TGraph(len(sign_array_v1),bdt_array, sign_array_v1 ) 
    graph_v1.GetXaxis().SetTitle('bdt cut')
    graph_v1.GetYaxis().SetTitle('S')
    graph_v1.SetTitle('')
    graph_v1.SetMarkerStyle(8)
    graph_v1.SetName('graph_signal_sample%s'%isample)

    graph_v2 = ROOT.TGraph(len(bkg_array_v1),bdt_array, bkg_array_v1 ) 
    graph_v2.GetXaxis().SetTitle('bdt cut')
    graph_v2.GetYaxis().SetTitle('B')
    graph_v2.SetTitle('')
    graph_v2.SetMarkerStyle(8)
    graph_v2.SetName('graph_background_sample%s'%isample)

    outfile.cd()
    graph_v1.Write()
    graph_v2.Write()


outfile.Close()