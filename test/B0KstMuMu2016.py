#################
### Variables ###
#################
runDataMC          = 1 # 1 = Data; 2 = MC (Reco + Gen); 3 = MC (Gen)
useJSON            = False
printMsg           = False
HLT_LMNR           = True
HLT_Psi            = False
triggerProcessName = 'HLT' # 'GEN' or 'HLT' or 'RECO' or 'TEST' or ...

print "\n@@@ CMSSW run configuration flags @@@"
print "runDataMC          : ", runDataMC
print "useJSON            : ", useJSON
print "printMsg           : ", printMsg
print "dataset LMNR       : ", HLT_LMNR
print "dataset Charmonium : ", HLT_Psi
print "triggerProcessName : ", triggerProcessName


#####################
### CMSSW configs ###
#####################
import FWCore.ParameterSet.Config as cms
process = cms.Process('B0KSTMUMUNTUPLIZER')

process.load('FWCore.MessageLogger.MessageLogger_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.MessageLogger.suppressWarning = cms.untracked.vstring('B0KstMuMu')


#################
### GlobalTag ###
#################
if (runDataMC == 1):
    process.GlobalTag.globaltag = cms.string('80X_dataRun2_2016SeptRepro_v7') 
#     process.GlobalTag.globaltag = cms.string('80X_dataRun2_Prompt_v16')
else:
    process.GlobalTag.globaltag = cms.string('80X_mcRun2_asymptotic_2016_TrancheIV_v6')


process.source = cms.Source('PoolSource',
                            duplicateCheckMode = cms.untracked.string('noDuplicateCheck'),
                            skipEvents = cms.untracked.uint32(0),
                            fileNames = cms.untracked.vstring(
                              'root://cms-xrd-global.cern.ch//store/data/Run2016G/Charmonium/AOD/23Sep2016-v1/110000/0E66894A-E39B-E611-A30C-848F69FD45CB.root'
                            ))


##################################
### Use JSON file interatively ###
##################################
# if (runDataMC == 1 and useJSON == True):
#     import FWCore.PythonUtilities.LumiList as LumiList
#     import FWCore.ParameterSet.Types as CfgTypes
#     if (run2016 == True):
#         myLumis = LumiList.LumiList(filename = 'Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_MuonPhys.txt').getCMSSWString().split(',')
#     else:
#         myLumis = LumiList.LumiList(filename = 'Cert_271036-284044_13TeV_PromptReco_Collisions16_JSON_MuonPhys.txt').getCMSSWString().split(',')
#     process.source.lumisToProcess = CfgTypes.untracked(CfgTypes.VLuminosityBlockRange())
#     process.source.lumisToProcess.extend(myLumis)


###################
### Output file ###
###################
process.maxEvents = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.maxLuminosityBlocks = cms.untracked.PSet(input = cms.untracked.int32(-1))
process.TFileService = cms.Service('TFileService', fileName = cms.string(
    'B0ToKstMuMu.root'
    ), closeFileFast = cms.untracked.bool(True))


###################################
### Import PAT skeleton process ###
###################################
process.load("PhysicsTools.PatAlgos.patSequences_cff")


############################
### Add track candidates ###
############################
process.load("PhysicsTools.PatAlgos.cleaningLayer1.genericTrackCleaner_cfi")
# Require NOT to check overlap with muons and electrons
process.cleanPatTracks.checkOverlaps.muons.requireNoOverlaps     = cms.bool(False)
process.cleanPatTracks.checkOverlaps.electrons.requireNoOverlaps = cms.bool(False)


#########################
### Set up PAT tracks ###
#########################
from PhysicsTools.PatAlgos.tools.trackTools import *
from PhysicsTools.PatAlgos.tools.coreTools  import *

# import pdb ; pdb.set_trace()
makeTrackCandidates(process,
                    label        = 'TrackCands',                  # output collection
                    tracks       = cms.InputTag('generalTracks'), # input track collection
                    particleType = 'pi+',                         # particle type (for assigning a mass)
                    preselection = 'pt > 0.5',                    # preselection cut on candidates
                    selection    = 'pt > 0.5',                    # selection cut on candidates
                    isolation    = {},                            # isolations to use (set to {} for None)
                    isoDeposits  = [],
                    mcAs         = None)                          # replicate MC match as the one used for Muons


removeMCMatching(process, ['All'], outputModules = [])


#####################################
### Do trigger matching for muons ###
#####################################

if (HLT_LMNR):
	process.cleanMuonTriggerMatchHLT = cms.EDProducer(
		'PATTriggerMatcherDRLessByR',
		src                   = cms.InputTag('cleanPatMuons'),
		matched               = cms.InputTag('patTrigger'),
		matchedCuts           = cms.string('path("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v*")'),
		maxDeltaR             = cms.double(0.1),
		resolveAmbiguities    = cms.bool(True),
		resolveByMatchQuality = cms.bool(True))
	process.cleanTrackTriggerMatchHLT = cms.EDProducer(
		'PATTriggerMatcherDRLessByR',
		src                   = cms.InputTag('cleanPatTrackCands'),
		matched               = cms.InputTag('patTrigger'),
		matchedCuts           = cms.string('path("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v*")'),
		maxDeltaR             = cms.double(0.1),
		resolveAmbiguities    = cms.bool(True),
		resolveByMatchQuality = cms.bool(True))
elif (HLT_Psi):
	process.cleanMuonTriggerMatchHLT = cms.EDProducer(
		'PATTriggerMatcherDRLessByR',
		src                   = cms.InputTag('cleanPatMuons'),
		matched               = cms.InputTag('patTrigger'),
		matchedCuts           = cms.string('path("HLT_DoubleMu4_JpsiTrk_Displaced_v*") || path("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v*")'),
		maxDeltaR             = cms.double(0.1),
		resolveAmbiguities    = cms.bool(True),
		resolveByMatchQuality = cms.bool(True))
	process.cleanTrackTriggerMatchHLT = cms.EDProducer(
		'PATTriggerMatcherDRLessByR',
		src                   = cms.InputTag('cleanPatTrackCands'),
		matched               = cms.InputTag('patTrigger'),
		matchedCuts           = cms.string('path("HLT_DoubleMu4_JpsiTrk_Displaced_v*") || path("HLT_DoubleMu4_PsiPrimeTrk_Displaced_v*")'),
		maxDeltaR             = cms.double(0.1),
		resolveAmbiguities    = cms.bool(True),
		resolveByMatchQuality = cms.bool(True))


##################################
### Switch on PAT trigger info ###
##################################
from PhysicsTools.PatAlgos.tools.trigTools import *
switchOnTriggerMatchEmbedding(process, triggerMatchers = ['cleanMuonTriggerMatchHLT', 'cleanTrackTriggerMatchHLT' ], hltProcess = triggerProcessName, outputModule = '')


#########################
### Remove unused PAT ###
#########################
process.patDefaultSequence.remove(process.muonMatch)
process.patDefaultSequence.remove(process.electronMatch)
process.patDefaultSequence.remove(process.photonMatch)
process.patDefaultSequence.remove(process.tauMatch)
process.patDefaultSequence.remove(process.tauGenJetMatch)

process.patDefaultSequence.remove(process.tauGenJets)
process.patDefaultSequence.remove(process.tauGenJetsSelectorAllHadrons)

process.patDefaultSequence.remove(process.patJetPartonMatch)
process.patDefaultSequence.remove(process.patJetGenJetMatch)
process.patDefaultSequence.remove(process.patJetCorrFactors)
process.patDefaultSequence.remove(process.patJetCharge)
process.patDefaultSequence.remove(process.patJetPartons)
process.patDefaultSequence.remove(process.patJetPartonsLegacy)            
process.patDefaultSequence.remove(process.patJetPartonAssociationLegacy)  
process.patDefaultSequence.remove(process.patJetFlavourAssociationLegacy) 
process.patDefaultSequence.remove(process.patJetFlavourAssociation)
process.patDefaultSequence.remove(process.patJets)

process.patDefaultSequence.remove(process.selectedPatJets)
process.patDefaultSequence.remove(process.cleanPatJets)
process.patDefaultSequence.remove(process.countPatJets)

process.patDefaultSequence.remove(process.selectedPatTaus)
process.patDefaultSequence.remove(process.cleanPatTaus)
process.patDefaultSequence.remove(process.countPatTaus)

process.patDefaultSequence.remove(process.selectedPatElectrons)
process.patDefaultSequence.remove(process.cleanPatElectrons)
process.patDefaultSequence.remove(process.countPatElectrons)

process.patDefaultSequence.remove(process.selectedPatPhotons)
process.patDefaultSequence.remove(process.cleanPatPhotons)
process.patDefaultSequence.remove(process.countPatPhotons)

process.patDefaultSequence.remove(process.patMETs)

## filter for good PVs
process.primaryVertexFilter = cms.EDFilter(
	"VertexSelector",
	src = cms.InputTag("offlinePrimaryVertices"),
	cut = cms.string("!isFake & ndof > 4 & abs(z) <= 24 & position.Rho <= 2"),
	filter = cms.bool(True)
	)


# ##########################
# ### B0 --> K*0 mu+ mu- ###
# ##########################
process.B0KstMuMu = cms.EDAnalyzer('B0KstMuMu',
                                   HLTriggerResults = cms.untracked.string(triggerProcessName),
                                   VtxSample        = cms.InputTag("offlinePrimaryVertices"),
                                   BeamSpot         = cms.InputTag("offlineBeamSpot"),
                                   GenParticles     = cms.untracked.InputTag('genParticles'),
                                   MuonType         = cms.untracked.InputTag('cleanPatMuonsTriggerMatch'),
                                   TrackType        = cms.untracked.InputTag('cleanPatTrackCandsTriggerMatch'),
                                   TriggerResult    = cms.untracked.InputTag("TriggerResults::HLT"),
                                   l1results        = cms.InputTag("gtStage2Digis"),
                                   PuInfoTag        = cms.untracked.InputTag("addPileupInfo"),
                                   GenFilterTag     = cms.untracked.InputTag("genFilterEfficiencyProducer"),
                                   doGenReco        = cms.untracked.uint32(runDataMC),
                                   TriggerNames     = cms.vstring("HLT_DoubleMu4_LowMassNonResonantTrk_Displaced_v", 
                                                                  "HLT_DoubleMu4_JpsiTrk_Displaced_v",
                                                                  "HLT_DoubleMu4_PsiPrimeTrk_Displaced_v"
                                                                  ),
                                   L1Names  = cms.vstring("L1_DoubleMu_11_4", 
                                                          "L1_DoubleMu_12_5",
                                                          "L1_DoubleMu_10_0_dEta_Max1p8",
                                                          "L1_DoubleMu0er1p6_dEta_Max1p8", 
                                                          "L1_DoubleMu0er1p6_dEta_Max1p8_OS",
                                                          "L1_DoubleMu0er1p4_dEta_Max1p8_OS", ### -> should be added for the next round if any
                                                   ),
                                   ## HLT selections
                                   MuMuVtxCL        = cms.untracked.double(0.1),    # mu-mu Vtx CL [0.1]             
                                   MuMuLsBS         = cms.untracked.double(3.0),    # mu-mu L/sigma w/respect to BS [3.0]    
                                   DCAMuMu          = cms.untracked.double(0.5),    # mu-mu DCA w/respect to each other [0.5 cm]
                                   DCAMuBS          = cms.untracked.double(2.0),    # mu DCA w/respect to BS [2.0 cm] 
                                   cosAlphaMuMuBS   = cms.untracked.double(0.9),    # mu-mu cos(alpha) w/respect to BS [0.9]   
                                   MinMupT          = cms.untracked.double(4.0),    # mu min pT [4.0 GeV/c]
                                   MuEta            = cms.untracked.double(2.4),    # mu max eta [2.4]    
                                   MuMupT           = cms.untracked.double(6.9),    # mu-mu min pT [6.9 GeV/c]    
                                   MinMuMuMass      = cms.untracked.double(1.0),    # mu-mu min inv. mass [1.0 GeV/c2]    
                                   MaxMuMuMass      = cms.untracked.double(4.8),    # mu-mu max inv. mass [4.8 GeV/c2]    
                                   ## Cand pre-selections
                                   MinB0Mass        = cms.untracked.double(4.5),    # B0 mass lower limit [4.5 GeV/c2]   
                                   MaxB0Mass        = cms.untracked.double(6.5),    # B0 mass upper limit [6.5 GeV/c2]
                                   B0VtxCL          = cms.untracked.double(0.01),   # B0 Vtx CL [0.01]
                                   KstMass          = cms.untracked.double(3.0),    # K*0 (OR K*0bar) mass window sigma [3.0]
                                   HadDCASBS        = cms.untracked.double(0.8),    # hadron DCA/sigma w/respect to BS [0.8] (also in HLT, now is 2) 
                                   HadpT            = cms.untracked.double(0.8),    ####   # hadron min pT [0.8 GeV/c] (also in HLT)
                                   MaxB0RoughMass   = cms.untracked.double(25.),    # B0 mass upper limit  before performing the fit
                                   
                                   printMsg         = cms.untracked.bool(printMsg))

# ###########
# ### RUN ###
# ###########
process.options = cms.untracked.PSet(allowUnscheduled = cms.untracked.bool(True))

process.patPath  = cms.Path(process.patDefaultSequence)
process.ntupPath = cms.Path(process.primaryVertexFilter*process.B0KstMuMu)

if (runDataMC == 3):
    process.schedule = cms.Schedule(process.ntupPath)
else:
    process.schedule = cms.Schedule(process.patPath, process.ntupPath)
