// ##################################################
// # Description:                                   #
// # Make rootTuple for b --> s mu+ mu- analysis    #
// # Original Author:  Mauro Dinardo                #
// #         Created:  Mon Apr 27 09:53:19 MDT 2011 #
// ##################################################

#ifndef B0KSTMUMU_H
#define B0KSTMUMU_H

// System include files
#include <memory>

// User include files
#include "FWCore/Framework/interface/ConsumesCollector.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CompositeCandidate.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/GenericParticle.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/MuonReco/interface/MuonChamberMatch.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/HitPattern.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"

#include "SimDataFormats/GeneratorProducts/interface/GenFilterInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicParticleFitter.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticleFactoryFromTransientTrack.h"
#include "TrackingTools/PatternTools/interface/ClosestApproachInRPhi.h"

#include "TrackingTools/IPTools/interface/IPTools.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "RecoVertex/VertexPrimitives/interface/TransientVertex.h"


#include <vector>
#include <string>
#include <TVector3.h>
#include "TLorentzVector.h"
#include "TTree.h"

#include "Utils.h"
#include "B0KstMuMuTreeContent.h"


class B0KstMuMu : public edm::EDAnalyzer
{
 public:
  
  explicit B0KstMuMu (const edm::ParameterSet&);
  ~B0KstMuMu ();
  
  std::string getMuCat (reco::Muon const& muon);

  const reco::Candidate* skipOscillations (const reco::Candidate* Mother);

  const reco::Candidate* findDaughter (edm::Handle<reco::GenParticleCollection> genParticles,
				       std::vector<std::vector<unsigned int>* >* Dau,
				       unsigned int it);
  
  void searchForStableChargedDaughters (const reco::Candidate* Mother,
					unsigned int motherIndex,
					std::vector<std::vector<unsigned int>* >* posDau,
					std::vector<std::vector<unsigned int>* >* negDau);

  float calculateIsolation(ClosestApproachInRPhi , 
                           reco::TransientTrack, 
                           reco::TransientTrack, 
                           const ParticleMass,
                           float, float, float );

  std::pair<double,double> pionImpactParameter(reco::TransientTrack , reco::Vertex    );
  int findPV(int , const reco::VertexCollection* );

 private:

  virtual void beginJob ();
  virtual void analyze (const edm::Event&, const edm::EventSetup&);
  virtual void endJob ();
  virtual void endLuminosityBlock (const edm::LuminosityBlock& lumiBlock, const edm::EventSetup& iSetup);

  // ######################################
  // # Ntuplizer configuration parameters #
  // ######################################
  std::string hltTriggerResults_;

  edm::InputTag vtxSampleTag_;
  edm::EDGetTokenT<reco::VertexCollection> vtxSampleToken_;
  edm::InputTag beamSpotTag_;
  edm::EDGetTokenT<reco::BeamSpot> beamSpotToken_;
  edm::InputTag genParticlesTag_;
  edm::EDGetTokenT<reco::GenParticleCollection> genParticlesToken_;

  edm::InputTag muonTypeTag_;   
  edm::EDGetTokenT< std::vector<pat::Muon> > muonTypeToken_;
  edm::InputTag trackTypeTag_;   
  edm::EDGetTokenT< std::vector<pat::GenericParticle> > trackTypeToken_;
  edm::InputTag triggerResultTag_;
  edm::EDGetTokenT<edm::TriggerResults>   triggerResultToken_;
  edm::InputTag puTag_;
  edm::EDGetTokenT<std::vector< PileupSummaryInfo>> puToken_;
  edm::InputTag genFilterTag_;
  edm::EDGetTokenT<GenFilterInfo>   genFilterToken_;

  unsigned int doGenReco_;

  std::vector<std::string> TrigTable_;
  std::vector<std::string> l1Table_;

  edm::EDGetTokenT<GlobalAlgBlkBxCollection> l1results_;
  edm::EDGetTokenT<GlobalExtBlkBxCollection> l1ext_;

  // ####################
  // # HLT-trigger cuts #
  // ####################
  double CLMUMUVTX;
  double LSMUMUBS;
  double DCAMUMU;
  double DCAMUBS;
  double COSALPHAMUMUBS;
  double MUMINPT;
  double MUMAXETA;
  double MINMUMUPT;
  double MINMUMUINVMASS;
  double MAXMUMUINVMASS;

  // ######################
  // # Pre-selection cuts #
  // ######################
  double B0MASSLOWLIMIT;
  double B0MASSUPLIMIT;
  double CLB0VTX;
  double KSTMASSWINDOW;
  double HADDCASBS;
  double MINHADPT;
  double MAXB0PREMASS;

  bool printMsg;
  
  TTree* theTree;
  B0KstMuMuTreeContent* NTuple;
  Utils* Utility;
  l1t::L1TGlobalUtil *fGtUtil;
};

#endif
