// ##################################################
// # Description:                                   #
// # Make rootTuple for b --> s mu+ mu- analysis    #
// # Original Author:  Mauro Dinardo                #
// #         Created:  Mon Apr 27 09:53:19 MDT 2011 #
// ##################################################

// System include files
#include <memory>

// User include files
#include "../interface/B0KstMuMu.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoVertex/AdaptiveVertexFit/interface/AdaptiveVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/KinematicConstrainedVertexFitter.h"
#include "RecoVertex/KinematicFit/interface/TwoTrackMassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/MultiTrackKinematicConstraint.h"
#include "RecoVertex/KinematicFit/interface/MassKinematicConstraint.h"
#include "RecoVertex/KinematicFitPrimitives/interface/KinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/RefCountedKinematicParticle.h"
#include "RecoVertex/KinematicFitPrimitives/interface/TransientTrackKinematicParticle.h"

#include "TrackingTools/Records/interface/TransientTrackRecord.h"
#include "TrackingTools/IPTools/interface/IPTools.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Math/interface/LorentzVector.h"

#include "TMath.h"

#include <sstream>
#include <utility>


// ####################
// # Global constants #
// ####################
#define TRKMAXR 110.0 // [cm]
#define TRKMAXZ 280.0 // [cm]

#define PRIVTXNDOF  4.0
#define PRIVTXMAXZ 50.0 // [cm]
#define PRIVTXMAXR  2.0 // [cm]

#define MUVARTOLE 0.01  // [From 0 to 1]
#define HADVARTOLE 0.10 // [From 0 to 1]


// #######################
// # Truth matching cuts #
// #######################
#define RCUTMU 0.004 // [eta-phi]
#define RCUTTRK 0.3  // [eta-phi]


B0KstMuMu::B0KstMuMu (const edm::ParameterSet& iConfig) :
  hltTriggerResults_ ( iConfig.getUntrackedParameter<std::string>("HLTriggerResults", std::string("HLT")) ),
  vtxSampleTag_      ( iConfig.getParameter<edm::InputTag>("VtxSample")), 
    vtxSampleToken_     (consumes<reco::VertexCollection>(vtxSampleTag_)), 
  beamSpotTag_       ( iConfig.getParameter<edm::InputTag>("BeamSpot")), 
    beamSpotToken_      (consumes<reco::BeamSpot>(beamSpotTag_)), 
  genParticlesTag_   ( iConfig.getUntrackedParameter<edm::InputTag>("GenParticles")),
    genParticlesToken_  (consumes<reco::GenParticleCollection>(genParticlesTag_)), 
  muonTypeTag_       ( iConfig.getUntrackedParameter<edm::InputTag>("MuonType")),
    muonTypeToken_      (consumes< std::vector<pat::Muon> >(muonTypeTag_)), 
  trackTypeTag_      ( iConfig.getUntrackedParameter<edm::InputTag>("TrackType")),
    trackTypeToken_     (consumes< std::vector<pat::GenericParticle> >(trackTypeTag_)), 
  triggerResultTag_  (iConfig.getUntrackedParameter<edm::InputTag>("TriggerResult")), 
    triggerResultToken_ (consumes<edm::TriggerResults>(triggerResultTag_)),
  puTag_             (iConfig.getUntrackedParameter<edm::InputTag>("PuInfoTag")),
    puToken_            (consumes<std::vector< PileupSummaryInfo>>(puTag_)), 
  genFilterTag_      (iConfig.getUntrackedParameter<edm::InputTag>("GenFilterTag")),
    genFilterToken_     (consumes<GenFilterInfo,edm::InLumi>(genFilterTag_)),

  doGenReco_         ( iConfig.getUntrackedParameter<unsigned int>("doGenReco",       1) ),
  TrigTable_         ( iConfig.getParameter<std::vector<std::string> >("TriggerNames")),   

  // # Load HLT-trigger cuts #
  CLMUMUVTX          ( iConfig.getUntrackedParameter<double>("MuMuVtxCL")      ),
  LSMUMUBS           ( iConfig.getUntrackedParameter<double>("MuMuLsBS")       ),
  DCAMUMU            ( iConfig.getUntrackedParameter<double>("DCAMuMu")        ),
  DCAMUBS            ( iConfig.getUntrackedParameter<double>("DCAMuBS")        ),
  COSALPHAMUMUBS     ( iConfig.getUntrackedParameter<double>("cosAlphaMuMuBS") ),
  MUMINPT            ( iConfig.getUntrackedParameter<double>("MinMupT")        ),
  MUMAXETA           ( iConfig.getUntrackedParameter<double>("MuEta")          ),
  MINMUMUPT          ( iConfig.getUntrackedParameter<double>("MuMupT")         ),
  MINMUMUINVMASS     ( iConfig.getUntrackedParameter<double>("MinMuMuMass")    ),
  MAXMUMUINVMASS     ( iConfig.getUntrackedParameter<double>("MaxMuMuMass")    ),
 
  // # Load pre-selection cuts #
  B0MASSLOWLIMIT     ( iConfig.getUntrackedParameter<double>("MinB0Mass")      ),
  B0MASSUPLIMIT      ( iConfig.getUntrackedParameter<double>("MaxB0Mass")      ),
  CLB0VTX            ( iConfig.getUntrackedParameter<double>("B0VtxCL")        ),
  KSTMASSWINDOW      ( iConfig.getUntrackedParameter<double>("KstMass")        ),
  HADDCASBS          ( iConfig.getUntrackedParameter<double>("HadDCASBS")      ),
  MINHADPT           ( iConfig.getUntrackedParameter<double>("HadpT")          ),
  MAXB0PREMASS       ( iConfig.getUntrackedParameter<double>("MaxB0RoughMass") ),
  

  printMsg           ( iConfig.getUntrackedParameter<bool>("printMsg",                false) ),
  
  theTree(0)
{
  std::cout << "\n@@@ Ntuplizer configuration parameters @@@" << std::endl;
  std::cout << __LINE__ << " : hltTriggerResults    = " << hltTriggerResults_ << std::endl;
  std::cout << __LINE__ << " : doGenReco     = " << doGenReco_     << std::endl;
  std::cout << __LINE__ << " : printMsg      = " << printMsg       << std::endl;
  for (auto itrig : TrigTable_ ) std::cout << __LINE__ << " : HLT paths     = " << itrig << std::endl;
  
  NTuple = new B0KstMuMuTreeContent();
  NTuple->Init();
  NTuple->ClearNTuple();

  Utility = new Utils();

}


B0KstMuMu::~B0KstMuMu ()
{
  delete NTuple;
  delete Utility;
}


void B0KstMuMu::analyze (const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  // ######################
  // # Internal Variables #
  // ######################
  std::stringstream myString;
  std::string tmpString1, tmpString2, tmpString3, tmpString4;

  std::string MuMCat;
  std::string MuPCat;

  const ParticleMass muonMass = Utility->muonMass;
  const ParticleMass pionMass = Utility->pionMass;
  const ParticleMass kaonMass = Utility->kaonMass;
  float muonMassErr           = Utility->muonMassErr;
  float pionMassErr           = Utility->pionMassErr;
  float kaonMassErr           = Utility->kaonMassErr;
  
  double chi;
  double ndf;

  double LSVtx;
  double LSVtxErr;
  double LSBS;
  double LSBSErr;

  double cosAlphaVtx;
  double cosAlphaVtxErr;
  double cosAlphaBS;
  double cosAlphaBSErr;

  double deltaEtaPhi;
  double pT;
  double eta;

  KinematicParticleFactoryFromTransientTrack partFactory;

  AdaptiveVertexFitter theVtxFitter;                              // Vertex fitter in nominal reconstruction
  KinematicParticleVertexFitter PartVtxFitter;                    // Vertex fit with vtx constraint

  ClosestApproachInRPhi ClosestApp;
  GlobalPoint XingPoint;

  reco::TrackRef muTrackTmp;
  reco::TrackRef muTrackm;
  reco::TrackRef muTrackp;
  reco::TrackRef Trackm;
  reco::TrackRef Trackp;

  std::pair <bool,Measurement1D> theDCAXVtx;
  TrajectoryStateClosestToPoint theDCAXBS;
  
  TrajectoryStateClosestToPoint traj;
  double mumMind0  = 100;
  double mupMind0  = 100;
  double TrkmMind0 = 100;
  double TrkpMind0 = 100;
  double mumMind0E, mupMind0E, TrkmMind0E, TrkpMind0E;
  double mumMinIP  = 100;
  double mupMinIP  = 100;
  double TrkmMinIP = 100;
  double TrkpMinIP = 100;
  double mumMinIPE, mupMinIPE, TrkmMinIPE, TrkpMinIPE;
  GlobalPoint vert;

  
  TLorentzVector  mum_lv;
  TLorentzVector  mup_lv;
  TLorentzVector jpsi_lv;

  TLorentzVector tkm_lv;
  TLorentzVector tkp_lv;
  TLorentzVector kst_lv;
  TLorentzVector kstbar_lv;

  std::vector<float> mum_isovec, mup_isovec, trkm_isovec, trkp_isovec; 


  if (printMsg) std::cout << "\n\n" << __LINE__ << " : @@@@@@ Start Analyzer @@@@@@" << std::endl;

  // Get Gen-Particles
  edm::Handle<reco::GenParticleCollection> genParticles;
  if ((doGenReco_ == 2) || (doGenReco_ == 3)) iEvent.getByToken(genParticlesToken_, genParticles);

  // Get magnetic field
  edm::ESHandle<MagneticField> bFieldHandle;
  iSetup.get<IdealMagneticFieldRecord>().get(bFieldHandle);      
  
  // Get HLT results
  edm::Handle<edm::TriggerResults> hltTriggerResults;
  try
  {
    iEvent.getByToken(triggerResultToken_, hltTriggerResults);
  }
  catch ( ... )
  {
    if (printMsg) std::cout << __LINE__ << " : couldn't get handle on HLT Trigger" << std::endl;
  }
  

  // ############################
  // # Save trigger information #
  // ############################
  HLTConfigProvider hltConfig_;
  bool changed = true;
  if (((hltTriggerResults.isValid() == false) || (hltConfig_.init(iEvent.getRun(), iSetup, hltTriggerResults_, changed) == false)) &&
       (printMsg)) std::cout << __LINE__ << " : no trigger results" << std::endl;
  else
  {
    // Get hold of trigger names - based on TriggerResults object
    const edm::TriggerNames& triggerNames_ = iEvent.triggerNames(*hltTriggerResults);
      
    for (unsigned int itrig = 0; itrig < hltTriggerResults->size(); itrig++)
    {
      if ((*hltTriggerResults)[itrig].accept() == 1)
      {
        std::string trigName = triggerNames_.triggerName(itrig);
        int trigPrescale = hltConfig_.prescaleValue(itrig, trigName);
          
        if (printMsg) std::cout << __LINE__ << " : Trigger name in the event: "  << trigName << "\twith prescale: " << trigPrescale << std::endl;
          
        // ################################
        // # Save HLT trigger information #
        // ################################
        for (unsigned int it = 0; it < TrigTable_.size(); it++){
          if (trigName.find(TrigTable_[it]) != std::string::npos)
          {
            NTuple->TrigTable->push_back(trigName);
            NTuple->TrigPrescales->push_back(trigPrescale);
               break;
          }
        }
      }  
    }
    if (NTuple->TrigTable->size() == 0)
    {
      NTuple->TrigTable->push_back("NotInTable");
      NTuple->TrigPrescales->push_back(-1);
    }
  }


  if ((doGenReco_ == 1) || (doGenReco_ == 2))
  {
    // Get BeamSpot
    edm::Handle<reco::BeamSpot> beamSpotH;
    iEvent.getByToken(beamSpotToken_, beamSpotH);
    reco::BeamSpot beamSpot = *beamSpotH;
      
    // Get primary vertex with BeamSpot constraint: ordered by the scalar sum of the |pT|^2 of the tracks that compose the vertex
    edm::Handle<reco::VertexCollection> recVtx;
    iEvent.getByToken(vtxSampleToken_, recVtx);
    reco::Vertex bestVtx = *(recVtx->begin());
    reco::Vertex bestVtxReFit;
      
    // Get PAT Muons
    edm::Handle< std::vector<pat::Muon> > thePATMuonHandle;
    iEvent.getByToken(muonTypeToken_, thePATMuonHandle);
      
    // Get PAT Tracks
    edm::Handle< std::vector<pat::GenericParticle> > thePATTrackHandle;
    iEvent.getByToken(trackTypeToken_, thePATTrackHandle);

    if (printMsg) std::cout << __LINE__ << " : the event has: " << thePATMuonHandle->size() << " muons AND " << thePATTrackHandle->size() << " tracks" << std::endl;


    // #################################
    // # Search for first valid vertex #
    // #################################
    for (std::vector<reco::Vertex>::const_iterator iVertex = recVtx->begin(); iVertex != recVtx->end(); iVertex++) { 
      bestVtx = *(iVertex); if (bestVtx.isValid() == true) break; 
    }


    if (bestVtx.isValid() == true)
    {
      // ###########
      // # Get mu- #
      // ###########
      for (std::vector<pat::Muon>::const_iterator iMuonM = thePATMuonHandle->begin(); iMuonM != thePATMuonHandle->end(); iMuonM++)
      {
        bool skip = false;

        // ########################
        // # Check mu- kinematics #
        // ########################
        muTrackm = iMuonM->innerTrack();
        if ((muTrackm.isNull() == true) || (muTrackm->charge() != -1)) continue;
 
        // #########################
        // # Muon- pT and eta cuts #
        // #########################
        pT  = muTrackm -> pt();
        eta = muTrackm -> eta();
        if ((pT < (MUMINPT*(1.0-MUVARTOLE))) || (fabs(eta) > (MUMAXETA*(1.0+MUVARTOLE))))
        {
          if (printMsg) std::cout << __LINE__ << " : break --> too low pT of mu- : " << pT << " or too high eta : " << eta << std::endl;
          break;
        }

        const reco::TransientTrack muTrackmTT(muTrackm, &(*bFieldHandle));
         
        // ###############################
        // # Compute mu- DCA to BeamSpot #
        // ###############################
        theDCAXBS = muTrackmTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
        if (theDCAXBS.isValid() == false)
        {
          if (printMsg) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for mu-" << std::endl;
          continue;
        }
        double DCAmumBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
        double DCAmumBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
        if (fabs(DCAmumBS) > DCAMUBS)
        {
          if (printMsg) std::cout << __LINE__ << " : continue --> bad absolute impact parameter 2D for mu- : " << DCAmumBS << std::endl;
          continue;
        }


        // ###########
        // # Get mu+ #
        // ###########
        for (std::vector<pat::Muon>::const_iterator iMuonP = thePATMuonHandle->begin(); iMuonP != thePATMuonHandle->end(); iMuonP++)
        {
          // ########################
          // # Check mu+ kinematics #
          // ########################
          muTrackp = iMuonP->innerTrack();
          if ((muTrackp.isNull() == true) || (muTrackp->charge() != 1)) continue;

          // #########################
          // # Muon+ pT and eta cuts #
          // #########################
          pT  = muTrackp->pt();
          eta = muTrackp->eta();
          if ((pT < (MUMINPT*(1.0-MUVARTOLE))) || (fabs(eta) > (MUMAXETA*(1.0+MUVARTOLE))))
          {
            if (printMsg) std::cout << __LINE__ << " : break --> too low pT of mu+ : " << pT << " or too high eta : " << eta << std::endl;
            break;
          }

          const reco::TransientTrack muTrackpTT(muTrackp, &(*bFieldHandle));

          // ###############################
          // # Compute mu+ DCA to BeamSpot #
          // ###############################
          theDCAXBS = muTrackpTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
          if (theDCAXBS.isValid() == false)
          {
            if (printMsg) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for mu+" << std::endl;
            continue;
          }
          double DCAmupBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
          double DCAmupBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
          if (fabs(DCAmupBS) > DCAMUBS)
          {
            if (printMsg) std::cout << __LINE__ << " : continue --> bad absolute impact parameter 2D for mu+: " << DCAmupBS << std::endl;
            continue;
          }


          // ############################################
          // # Check goodness of muons closest approach #
          // ############################################
          ClosestApp.calculate(muTrackpTT.initialFreeState(),muTrackmTT.initialFreeState());
          if (ClosestApp.status() == false)
          {
            if (printMsg) std::cout << __LINE__ << " : continue --> bad status of closest approach" << std::endl;
            continue;
          }
          XingPoint = ClosestApp.crossingPoint();
          if ((sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y()) > TRKMAXR) || (fabs(XingPoint.z()) > TRKMAXZ))
          {
            if (printMsg) std::cout << __LINE__ << " : continue --> closest approach crossing point outside the tracker volume" << std::endl;
            continue;
          }


          // #####################################################
          // # Cut on the mumu 3D-DCA with respect to each other #
          // #####################################################
          double mumuDCA = ClosestApp.distance();
          if (mumuDCA > DCAMUMU)
          {
            if (printMsg) std::cout << __LINE__ << " : continue --> bad 3D-DCA of mu+(-) with respect to mu-(+): " << mumuDCA << std::endl;
            continue;
          }


          // ############################################
          // # Cut on the dimuon inviariant mass and pT #
          // ############################################
 
          mum_lv.SetPxPyPzE( muTrackmTT.track().px(), 
                             muTrackmTT.track().py(),
                             muTrackmTT.track().pz(),
                             sqrt( pow(muTrackmTT.track().p(),2) + pow(Utility->muonMass,2) ));
          mup_lv.SetPxPyPzE( muTrackpTT.track().px(), 
                             muTrackpTT.track().py(),
                             muTrackpTT.track().pz(),
                             sqrt( pow(muTrackpTT.track().p(),2) + pow(Utility->muonMass,2) ));
          jpsi_lv = mup_lv + mum_lv;
 
          if ((jpsi_lv.Pt() < (MINMUMUPT*(1.0-MUVARTOLE))) || (jpsi_lv.M() < (MINMUMUINVMASS*(1.0-MUVARTOLE))) || (jpsi_lv.M() > (MAXMUMUINVMASS*(1.0+MUVARTOLE))))
          {
            if (printMsg) std::cout << __LINE__ << " : continue --> no good mumu pair pT: " << jpsi_lv.Pt() << "\tinv. mass: " << jpsi_lv.M() << std::endl;
            continue;
          }


          // #######################################################
          // # @@@ Make mu-mu and implement pre-selection cuts @@@ #
          // #######################################################
          if (printMsg) std::cout << "\n" << __LINE__ << " : @@@ I have 2 good oppositely-charged muons. I'm trying to vertex them @@@" << std::endl;


          chi = 0.;
          ndf = 0.;
          // ####################################################
          // # Try to vertex the two muons to get dimuon vertex #
          // ####################################################
          std::vector<RefCountedKinematicParticle> muonParticles;
          muonParticles.push_back(partFactory.particle(muTrackmTT,muonMass,chi,ndf,muonMassErr));
          muonParticles.push_back(partFactory.particle(muTrackpTT,muonMass,chi,ndf,muonMassErr));
      
          RefCountedKinematicTree mumuVertexFitTree = PartVtxFitter.fit(muonParticles);
          if (mumuVertexFitTree->isValid() == false)
          {
            if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the mu+ mu- vertex fit" << std::endl;
            continue; 
          }
      
          mumuVertexFitTree->movePointerToTheTop();
          RefCountedKinematicVertex mumu_KV   = mumuVertexFitTree->currentDecayVertex();
          if (mumu_KV->vertexIsValid() == false)
          {
            if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the mu+ mu- vertex fit" << std::endl;
            continue;
          }
          if (TMath::Prob(static_cast<double>(mumu_KV->chiSquared()), static_cast<int>(rint(mumu_KV->degreesOfFreedom()))) < CLMUMUVTX)
            {
            if (printMsg) std::cout << __LINE__ << " : continue --> bad vtx CL from mu+ mu- fit: " << TMath::Prob(static_cast<double>(mumu_KV->chiSquared()), static_cast<int>(rint(mumu_KV->degreesOfFreedom()))) << std::endl;
            continue;
          }


          // #########################################################
          // # Extract the re-fitted tracks after the dimuon vtx fit #
          // #########################################################
          RefCountedKinematicParticle mumu_KP = mumuVertexFitTree->currentParticle();
          mumuVertexFitTree->movePointerToTheTop();

          mumuVertexFitTree->movePointerToTheFirstChild();
          RefCountedKinematicParticle refitMum  = mumuVertexFitTree->currentParticle();
          const reco::TransientTrack refitMumTT = refitMum->refittedTransientTrack();

          mumuVertexFitTree->movePointerToTheNextChild();
          RefCountedKinematicParticle refitMup  = mumuVertexFitTree->currentParticle();
          const reco::TransientTrack refitMupTT = refitMup->refittedTransientTrack();


          // ########################
          // # Muon pT and eta cuts #
          // ########################
          pT  = refitMupTT.track().pt();
          eta = refitMupTT.track().eta();
          if ((pT < MUMINPT) || (fabs(eta) > MUMAXETA))
          {
            if (printMsg) std::cout << __LINE__ << " : continue --> too low pT of mu+ : " << pT << " or too high eta : " << eta << std::endl;
            continue;
          }

          pT  = refitMumTT.track().pt();
          eta = refitMumTT.track().eta();
          if ((pT < MUMINPT) || (fabs(eta) > MUMAXETA))
          {
            if (printMsg) std::cout << __LINE__ << " : continue --> too low pT of mu- : " << pT << " or too high eta : " << eta << std::endl;
            skip = true;
            continue;
          }


          // ############################################
          // # Cut on the dimuon invariant mass and pT #
          // ############################################
          pT = sqrt((refitMumTT.track().momentum().x() + refitMupTT.track().momentum().x()) * (refitMumTT.track().momentum().x() + refitMupTT.track().momentum().x()) +
                    (refitMumTT.track().momentum().y() + refitMupTT.track().momentum().y()) * (refitMumTT.track().momentum().y() + refitMupTT.track().momentum().y()));
          double MuMuInvMass = mumu_KP->currentState().mass();
          if ((pT < MINMUMUPT) || (MuMuInvMass < MINMUMUINVMASS) || (MuMuInvMass > MAXMUMUINVMASS))
          {
            if (printMsg) std::cout << __LINE__ << " : continue --> no good mumu pair pT: " << pT << "\tinv. mass: " << MuMuInvMass << std::endl;
            continue;
          }


          // ######################################################
          // # Compute the distance between mumu vtx and BeamSpot #
          // ######################################################
          double MuMuLSBS;
          double MuMuLSBSErr;
          Utility->computeLS (mumu_KV->position().x(),mumu_KV->position().y(),0.0,
                      beamSpot.position().x(),beamSpot.position().y(),0.0,
                      mumu_KV->error().cxx(),mumu_KV->error().cyy(),0.0,
                      mumu_KV->error().matrix()(0,1),0.0,0.0,
                      beamSpot.covariance()(0,0),beamSpot.covariance()(1,1),0.0,
                      beamSpot.covariance()(0,1),0.0,0.0,
                      &MuMuLSBS,&MuMuLSBSErr);
          if (MuMuLSBS/MuMuLSBSErr < LSMUMUBS)     
          {
            if (printMsg) std::cout << __LINE__ << " : continue --> bad mumu L/sigma with respect to BeamSpot: " << MuMuLSBS << "+/-" << MuMuLSBSErr << std::endl;
            continue;
          }
      

          // ###################################################################
          // # Compute cos(alpha) between mumu momentum and mumuVtx - BeamSpot #
          // ###################################################################
          double MuMuCosAlphaBS;
          double MuMuCosAlphaBSErr;
          Utility->computeCosAlpha (mumu_KP->currentState().globalMomentum().x(),mumu_KP->currentState().globalMomentum().y(),0.0,
                                    mumu_KV->position().x() - beamSpot.position().x(),mumu_KV->position().y() - beamSpot.position().y(),0.0,
                                    mumu_KP->currentState().kinematicParametersError().matrix()(3,3),mumu_KP->currentState().kinematicParametersError().matrix()(4,4),0.0,
                                    mumu_KP->currentState().kinematicParametersError().matrix()(3,4),0.0,0.0,
                                    mumu_KV->error().cxx() + beamSpot.covariance()(0,0),mumu_KV->error().cyy() + beamSpot.covariance()(1,1),0.0,
                                    mumu_KV->error().matrix()(0,1) + beamSpot.covariance()(0,1),0.0,0.0,
                                    &MuMuCosAlphaBS,&MuMuCosAlphaBSErr);
          if (MuMuCosAlphaBS < COSALPHAMUMUBS)
          {
            if (printMsg) std::cout << __LINE__ << " : continue --> bad mumu cos(alpha) with respect to BeamSpot: " << MuMuCosAlphaBS << "+/-" << MuMuCosAlphaBSErr << std::endl;
            continue;
          }


          // ########################### convert KinFit vertex to reco Vertex
          reco::Vertex::Point mumu_GP  = reco::Vertex::Point(mumu_KV->position().x(), mumu_KV->position().y(), mumu_KV->position().z());
          const reco::Vertex::Error mumu_error = mumu_KV->vertexState().error().matrix();
          float mumu_chi2      = mumu_KV -> chiSquared();
          float mumu_ndof      = mumu_KV -> degreesOfFreedom();
          reco::Vertex mumu_rv =  reco::Vertex( mumu_GP, mumu_error, mumu_chi2, mumu_ndof, 2 );

          // ##############
          // # Get Track- #
          // ##############
          for (std::vector<pat::GenericParticle>::const_iterator iTrackM = thePATTrackHandle->begin(); iTrackM != thePATTrackHandle->end(); iTrackM++)
            {
              bool skip = false;

              // ###########################
              // # Check Track- kinematics #
              // ###########################
              Trackm = iTrackM->track();
              if ((Trackm.isNull() == true) || (Trackm->charge() != -1)) continue;

              const reco::TransientTrack TrackmTT(Trackm, &(*bFieldHandle));

              // ##########################
              // # Track- pT and eta cuts #
              // ##########################
              pT = TrackmTT.track().pt();
              if (pT < (MINHADPT*(1.0-HADVARTOLE)))
              {
                if (printMsg) std::cout << __LINE__ << " : break --> too low pT of track- : " << pT << std::endl;
                break;
              }
              
              // ######################################
              // # Compute K*0 track- DCA to BeamSpot #
              // ######################################
              theDCAXBS = TrackmTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
              if (theDCAXBS.isValid() == false)
              {
                if (printMsg) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for track-" << std::endl;
                continue;
              }
              double DCAKstTrkmBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
              double DCAKstTrkmBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
              if (fabs(DCAKstTrkmBS/DCAKstTrkmBSErr) < HADDCASBS)
              {
                if (printMsg) std::cout << __LINE__ << " : continue --> track- DCA/sigma with respect to BeamSpot is too small: " << DCAKstTrkmBS << "+/-" << DCAKstTrkmBSErr << std::endl;
                continue;
              }
              
              
              // ##############
              // # Get Track+ #
              // ##############
              for (std::vector<pat::GenericParticle>::const_iterator iTrackP = thePATTrackHandle->begin(); iTrackP != thePATTrackHandle->end(); iTrackP++)
              {
                // ###########################
                // # Check Track+ kinematics #
                // ###########################
                Trackp = iTrackP->track();
                if ((Trackp.isNull() == true) || (Trackp->charge() != 1)) continue;
  
                const reco::TransientTrack TrackpTT(Trackp, &(*bFieldHandle));

                // ##########################
                // # Track+ pT and eta cuts #
                // ##########################
                pT = TrackpTT.track().pt();
                if (pT < (MINHADPT*(1.0-HADVARTOLE)))
                {
                  if (printMsg) std::cout << __LINE__ << " : break --> too low pT of track+ : " << pT << std::endl;
                  break;
                }


                // ######################################
                // # Compute K*0 track+ DCA to BeamSpot #
                // ######################################
                theDCAXBS = TrackpTT.trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
                if (theDCAXBS.isValid() == false)
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for track+" << std::endl;
                  continue;
                }
                double DCAKstTrkpBS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
                double DCAKstTrkpBSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
                if (fabs(DCAKstTrkpBS/DCAKstTrkpBSErr) < HADDCASBS)
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> track+ DCA/sigma with respect to BeamSpot is too small: " << DCAKstTrkpBS << "+/-" << DCAKstTrkpBSErr << std::endl;
                  continue;
                }


                // ##############################################
                // # Check goodness of hadrons closest approach #
                // ##############################################
                ClosestApp.calculate(TrackpTT.initialFreeState(),TrackmTT.initialFreeState());
                if (ClosestApp.status() == false)
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> bad status of closest approach" << std::endl;
                  continue;
                }
                XingPoint = ClosestApp.crossingPoint();
                if ((sqrt(XingPoint.x()*XingPoint.x() + XingPoint.y()*XingPoint.y()) > TRKMAXR) || (fabs(XingPoint.z()) > TRKMAXZ))
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> closest approach crossing point outside the tracker volume" << std::endl;
                  continue;
                }


                // ######################################################
                // # Check if K*0 (OR K*0bar) mass is within acceptance #
                // ######################################################
                tkm_lv.SetPxPyPzE(TrackmTT.track().momentum().x(),
                                  TrackmTT.track().momentum().y(),
                                  TrackmTT.track().momentum().z(),
                                  sqrt( pow(TrackmTT.track().p(),2) + pow(Utility->pionMass,2) ));
                tkp_lv.SetPxPyPzE(TrackpTT.track().momentum().x(),
                                  TrackpTT.track().momentum().y(),
                                  TrackpTT.track().momentum().z(),
                                  sqrt( pow(TrackpTT.track().p(),2) + pow(Utility->kaonMass,2) ) );
                double kstInvMass = (tkm_lv + tkp_lv).M();
                
                tkm_lv.SetE(sqrt( pow(TrackmTT.track().p(),2) + pow(Utility->kaonMass,2) ));
                tkp_lv.SetE(sqrt( pow(TrackpTT.track().p(),2) + pow(Utility->pionMass,2) ));
                double kstBarInvMass = (tkm_lv + tkp_lv).M();
                
                if ((fabs(kstInvMass - Utility->kstMass)    > (KSTMASSWINDOW*Utility->kstSigma*(1.0+HADVARTOLE))) &&
                    (fabs(kstBarInvMass - Utility->kstMass) > (KSTMASSWINDOW*Utility->kstSigma*(1.0+HADVARTOLE))))
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> bad K*0 mass: " << kstInvMass << " AND K*0bar mass: " << kstBarInvMass << std::endl;
                  continue;
                }



                // ####################################################
                // # @@@ Make K* and implement pre-selection cuts @@@ #
                // ####################################################
                if (printMsg) std::cout << "\n" << __LINE__ << " : @@@ I have 2 good oppositely-charged tracks. I'm trying to vertex them @@@" << std::endl;

                chi = 0.;
                ndf = 0.;
                // ##############################################################################
                // # Try to vertex the two Tracks to get K*0 vertex: pion = track- | k = track+ #
                // ##############################################################################
                std::vector<RefCountedKinematicParticle> kstParticles;
                kstParticles.push_back(partFactory.particle(TrackmTT,pionMass,chi,ndf,pionMassErr));
                kstParticles.push_back(partFactory.particle(TrackpTT,kaonMass,chi,ndf,kaonMassErr));
                
                RefCountedKinematicTree kstVertexFitTree = PartVtxFitter.fit(kstParticles);
                if (kstVertexFitTree->isValid() == false)
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the K*0 vertex fit" << std::endl;
                  continue;
                }
                
                kstVertexFitTree->movePointerToTheTop();
                RefCountedKinematicParticle kst_KP = kstVertexFitTree->currentParticle();
                RefCountedKinematicVertex kst_KV   = kstVertexFitTree->currentDecayVertex();
                if (kst_KV->vertexIsValid() == false)
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the K*0 vertex fit" << std::endl;
                  continue;
                }


                chi = 0.;
                ndf = 0.;
                // #################################################################################
                // # Try to vertex the two Tracks to get K*0bar vertex: pion = track+ | k = track- #
                // #################################################################################
                std::vector<RefCountedKinematicParticle> kstBarParticles;
                kstBarParticles.push_back(partFactory.particle(TrackmTT,kaonMass,chi,ndf,kaonMassErr));
                kstBarParticles.push_back(partFactory.particle(TrackpTT,pionMass,chi,ndf,pionMassErr));
                
                RefCountedKinematicTree kstBarVertexFitTree = PartVtxFitter.fit(kstBarParticles);
                if (kstBarVertexFitTree->isValid() == false)
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the K*0bar vertex fit" << std::endl;
                  continue;
                }
                
                kstBarVertexFitTree->movePointerToTheTop();
                RefCountedKinematicParticle kstBar_KP = kstBarVertexFitTree->currentParticle();
                RefCountedKinematicVertex kstBar_KV   = kstBarVertexFitTree->currentDecayVertex();
                if (kstBar_KV->vertexIsValid() == false)
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the K*0bar vertex fit" << std::endl;
                  continue;
                }
  
                // ######################################################
                // # Check if K*0 (OR K*0bar) mass is within acceptance #
                // ######################################################
                kstInvMass    = kst_KP->currentState().mass();
                kstBarInvMass = kstBar_KP->currentState().mass();
                if ((fabs(kstInvMass - Utility->kstMass) > KSTMASSWINDOW*Utility->kstSigma) && (fabs(kstBarInvMass - Utility->kstMass) > KSTMASSWINDOW*Utility->kstSigma))
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> bad K*0 mass: " << kstInvMass << " AND K*0bar mass: " << kstBarInvMass << std::endl;
                  continue;
                }
  
  
                // ###########################################################
                // # Extract the re-fitted tracks after the dihadron vtx fit #
                // ###########################################################
                kstVertexFitTree->movePointerToTheTop();
  
                kstVertexFitTree->movePointerToTheFirstChild();
                RefCountedKinematicParticle refitTrkm  = kstVertexFitTree->currentParticle();
                const reco::TransientTrack refitTrkmTT = refitTrkm->refittedTransientTrack();
  
                kstVertexFitTree->movePointerToTheNextChild();
                RefCountedKinematicParticle refitTrkp  = kstVertexFitTree->currentParticle();
                const reco::TransientTrack refitTrkpTT = refitTrkp->refittedTransientTrack();


                // ##########################
                // # Hadron pT and eta cuts #
                // ##########################
                pT = refitTrkpTT.track().pt();
                if (pT < MINHADPT)
                {
                  if (printMsg) std::cout << __LINE__ << " : break --> too low pT of track+ : " << pT << std::endl;
                  continue;
                }
  
                pT = refitTrkmTT.track().pt();
                if (pT < MINHADPT)
                {
                  if (printMsg) std::cout << __LINE__ << " : break --> too low pT of track- : " << pT << std::endl;
                  skip = true;
                  continue;
                }
  
  
                // ####################################################
                // # @@@ Make B0 and implement pre-selection cuts @@@ #
                // ####################################################
                if (printMsg) std::cout << "\n" << __LINE__ << " : @@@ I have 4 good charged tracks. I'm trying to vertex them @@@" << std::endl;
  
                TLorentzVector a_lv, b_lv, c_lv, d_lv, tot_lv;
                a_lv.SetPxPyPzE(muTrackmTT.track().momentum().x(),
                                muTrackmTT.track().momentum().y(),
                                muTrackmTT.track().momentum().z(),
                                sqrt( pow(muTrackmTT.track().p(),2) + pow(Utility->muonMass,2) ));
                b_lv.SetPxPyPzE(muTrackpTT.track().momentum().x(),
                                muTrackpTT.track().momentum().y(),
                                muTrackpTT.track().momentum().z(),
                                sqrt( pow(muTrackpTT.track().p(),2) + pow(Utility->muonMass,2) ));
                c_lv.SetPxPyPzE(TrackmTT.track().momentum().x(),
                                TrackmTT.track().momentum().y(),
                                TrackmTT.track().momentum().z(),
                                sqrt( pow(TrackmTT.track().p(),2) + pow(Utility->pionMass,2) ));
                d_lv.SetPxPyPzE(TrackpTT.track().momentum().x(),
                                TrackpTT.track().momentum().y(),
                                TrackpTT.track().momentum().z(),
                                sqrt( pow(TrackpTT.track().p(),2) + pow(Utility->kaonMass,2) ));
                            
                tot_lv = a_lv + b_lv + c_lv + d_lv;
                if  (tot_lv.M() > MAXB0PREMASS) {
                  if (printMsg) std::cout << __LINE__ << " : continue --> b0 mass before fit is > max value" << std::endl;
                    continue;    
                }         
  
  
  
                // #################################################
                // # Check if the hadron Track- is actually a muon #
                // #################################################
                MuMCat.clear();
                MuMCat = "NotMatched";
                for (std::vector<pat::Muon>::const_iterator iMuon = thePATMuonHandle->begin(); iMuon != thePATMuonHandle->end(); iMuon++)
                {
                  muTrackTmp = iMuon->innerTrack();
                  if ((muTrackTmp.isNull() == true) || (muTrackTmp->charge() != -1)) continue;
                  if (Trackm == muTrackTmp)
                  {
                    MuMCat.clear();
                    MuMCat.append(getMuCat(*iMuon));
                    if (printMsg) std::cout << __LINE__ << " : negative charged hadron is actually a muon (momentum: " << Trackm->p() << ") whose category is: " << MuMCat.c_str() << std::endl;
                    break;
                  }
                }
  
  
                // #################################################
                // # Check if the hadron Track+ is actually a muon #
                // #################################################
                MuPCat.clear();
                MuPCat = "NotMatched";
                for (std::vector<pat::Muon>::const_iterator iMuon = thePATMuonHandle->begin(); iMuon != thePATMuonHandle->end(); iMuon++)
                {
                  muTrackTmp = iMuon->innerTrack();
                  if ((muTrackTmp.isNull() == true) || (muTrackTmp->charge() != 1)) continue;
                  if (Trackp == muTrackTmp)
                  {
                    MuPCat.clear();
                    MuPCat.append(getMuCat(*iMuon));
                    if (printMsg) std::cout << __LINE__ << " : positive charged hadron is actually a muon (momentum: " << Trackp->p() << ") whose category is: " << MuPCat.c_str() << std::endl;
                    break;
                  }
                }


                chi = 0.;
                ndf = 0.;
                // #####################################
                // # B0 vertex fit with vtx constraint #
                // #####################################
                std::vector<RefCountedKinematicParticle> bParticles;
                bParticles.push_back(partFactory.particle(muTrackmTT,muonMass,chi,ndf,muonMassErr));
                bParticles.push_back(partFactory.particle(muTrackpTT,muonMass,chi,ndf,muonMassErr));
                bParticles.push_back(partFactory.particle(TrackmTT,pionMass,chi,ndf,pionMassErr));
                bParticles.push_back(partFactory.particle(TrackpTT,kaonMass,chi,ndf,kaonMassErr));
  
                RefCountedKinematicTree bVertexFitTree = PartVtxFitter.fit(bParticles);
                if (bVertexFitTree->isValid() == false)
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the B0 vertex fit" << std::endl;
                  continue;
                }
  
                bVertexFitTree->movePointerToTheTop();
                RefCountedKinematicParticle b_KP = bVertexFitTree->currentParticle();
                RefCountedKinematicVertex b_KV   = bVertexFitTree->currentDecayVertex();
                if (b_KV->vertexIsValid() == false)
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the B0 vertex fit" << std::endl;
                  continue;
                }
  
  
                chi = 0.;
                ndf = 0.;
                // ########################################
                // # B0bar vertex fit with vtx constraint #
                // ########################################
                std::vector<RefCountedKinematicParticle> bBarParticles;
                bBarParticles.push_back(partFactory.particle(muTrackmTT,muonMass,chi,ndf,muonMassErr));
                bBarParticles.push_back(partFactory.particle(muTrackpTT,muonMass,chi,ndf,muonMassErr));
                bBarParticles.push_back(partFactory.particle(TrackmTT,kaonMass,chi,ndf,kaonMassErr));
                bBarParticles.push_back(partFactory.particle(TrackpTT,pionMass,chi,ndf,pionMassErr));
                    
                RefCountedKinematicTree bBarVertexFitTree = PartVtxFitter.fit(bBarParticles);
                if (bBarVertexFitTree->isValid() == false)
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the B0bar vertex fit" << std::endl;
                  continue;
                }
  
                bBarVertexFitTree->movePointerToTheTop();
                RefCountedKinematicParticle bBar_KP = bBarVertexFitTree->currentParticle();
                RefCountedKinematicVertex bBar_KV   = bBarVertexFitTree->currentDecayVertex();
                if (bBar_KV->vertexIsValid() == false)
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> invalid vertex from the B0bar vertex fit" << std::endl;
                  continue;
                }
  
                    // ########################
                // # Cuts on B0 AND B0bar #
                // ########################
                if (((b_KP->currentState().mass() < B0MASSLOWLIMIT) || (b_KP->currentState().mass() > B0MASSUPLIMIT)) &&
                  ((bBar_KP->currentState().mass() < B0MASSLOWLIMIT) || (bBar_KP->currentState().mass() > B0MASSUPLIMIT)))
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> bad B0 mass: " << b_KP->currentState().mass() << " AND B0bar mass: " << bBar_KP->currentState().mass() << std::endl;
                  continue;
                }
                if ((TMath::Prob(static_cast<double>(b_KV->chiSquared()), static_cast<int>(rint(b_KV->degreesOfFreedom()))) < CLB0VTX) &&
                    (TMath::Prob(static_cast<double>(bBar_KV->chiSquared()), static_cast<int>(rint(bBar_KV->degreesOfFreedom()))) < CLB0VTX))
                {
                  if (printMsg)
                  {
                    std::cout << __LINE__ << " : continue --> bad vtx CL from B0 fit: " << TMath::Prob(static_cast<double>(b_KV->chiSquared()), static_cast<int>(rint(b_KV->degreesOfFreedom())));
                    std::cout << " AND bad vtx CL from B0bar fit: " << TMath::Prob(static_cast<double>(bBar_KV->chiSquared()), static_cast<int>(rint(bBar_KV->degreesOfFreedom()))) << std::endl;
                  }
                  continue;
                }


                // ###############################
                // # Cuts on B0 L/sigma BeamSpot #
                // ###############################
                Utility->computeLS (b_KV->position().x(),b_KV->position().y(),0.0,
                                    beamSpot.position().x(),beamSpot.position().y(),0.0,
                                    b_KV->error().cxx(),b_KV->error().cyy(),0.0,
                                    b_KV->error().matrix()(0,1),0.0,0.0,
                                    beamSpot.covariance()(0,0),beamSpot.covariance()(1,1),0.0,
                                    beamSpot.covariance()(0,1),0.0,0.0,
                                    &LSBS,&LSBSErr);
  
  
                // ##############################
                // # Compute B0 DCA to BeamSpot #
                // ##############################
                theDCAXBS = b_KP->refittedTransientTrack().trajectoryStateClosestToPoint(GlobalPoint(beamSpot.position().x(),beamSpot.position().y(),beamSpot.position().z()));
                if (theDCAXBS.isValid() == false)
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 2D for B0" << std::endl;
                  continue;
                }      
                double DCAB0BS    = theDCAXBS.perigeeParameters().transverseImpactParameter();
                double DCAB0BSErr = theDCAXBS.perigeeError().transverseImpactParameterError();
  
          
                // #####################################
                // # Compute B0 cos(alpha) to BeamSpot #
                // #####################################
                Utility->computeCosAlpha (b_KP->currentState().globalMomentum().x(),b_KP->currentState().globalMomentum().y(),0.0,
                                          b_KV->position().x() - beamSpot.position().x(),b_KV->position().y() - beamSpot.position().y(),0.0,
                                          b_KP->currentState().kinematicParametersError().matrix()(3,3),b_KP->currentState().kinematicParametersError().matrix()(4,4),0.0,
                                          b_KP->currentState().kinematicParametersError().matrix()(3,4),0.0,0.0,
                                          b_KV->error().cxx() + beamSpot.covariance()(0,0),b_KV->error().cyy() + beamSpot.covariance()(1,1),0.0,
                                          b_KV->error().matrix()(0,1) + beamSpot.covariance()(0,1),0.0,0.0,
                                          &cosAlphaBS,&cosAlphaBSErr);
  
  
  
                // #############################################################################################################
                // # If the tracks in the B reco candidate belongs to the primary vertex, then remove them and redo the Vertex #
                // #############################################################################################################
                std::vector<reco::TransientTrack> vertexTracks;
                for (std::vector<reco::TrackBaseRef>::const_iterator iTrack = bestVtx.tracks_begin(); iTrack != bestVtx.tracks_end(); iTrack++)
                {
                  // Compare primary tracks to check for matches with B candidate
                  reco::TrackRef trackRef = iTrack->castTo<reco::TrackRef>();
                        
                  if ((Trackm != trackRef) && (Trackp != trackRef) && (muTrackm != trackRef) && (muTrackp != trackRef))
                  {
                    reco::TransientTrack TT(trackRef, &(*bFieldHandle));
                    vertexTracks.push_back(TT);
                  }
                  else if (printMsg) std::cout << __LINE__ << " : some of the B0 reco candidate tracks belong to the Pri.Vtx" << std::endl;
                }
  
                if (vertexTracks.size() < 2)
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> number of tracks of the new Pri.Vtx is too small: " << vertexTracks.size() << std::endl;
                  continue;
                }
                        
                TransientVertex bestTransVtx = theVtxFitter.vertex(vertexTracks);
                if (bestTransVtx.isValid() == false)
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> invalid new Pri.Vtx" << std::endl;
                  continue;
                }
                bestVtxReFit = (reco::Vertex)bestTransVtx;
  
  
  
  
                // ##########################
                // # Compute mu- DCA to Vtx #
                // ##########################
                theDCAXVtx = IPTools::absoluteImpactParameter3D(muTrackmTT, bestVtxReFit);
                if (theDCAXVtx.first == false)
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 3D for mu-" << std::endl;
                  continue;
                }
                double DCAmumVtx    = theDCAXVtx.second.value();
                double DCAmumVtxErr = theDCAXVtx.second.error();
              
                                    
                // ##########################
                // # Compute mu+ DCA to Vtx #
                // ##########################
                theDCAXVtx = IPTools::absoluteImpactParameter3D(muTrackpTT, bestVtxReFit);
                if (theDCAXVtx.first == false)
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 3D for mu+" << std::endl;
                  continue;
                }
                double DCAmupVtx    = theDCAXVtx.second.value();
                double DCAmupVtxErr = theDCAXVtx.second.error();
  
              
                // #################################
                // # Compute K*0 track- DCA to Vtx #
                // #################################
                theDCAXVtx = IPTools::absoluteImpactParameter3D(TrackmTT, bestVtxReFit);
                if (theDCAXVtx.first == false)
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 3D for track-" << std::endl;
                  continue;
                }
                double DCAKstTrkmVtx    = theDCAXVtx.second.value();
                double DCAKstTrkmVtxErr = theDCAXVtx.second.error();
  
  
                // #################################
                // # Compute K*0 track+ DCA to Vtx #
                // #################################
                theDCAXVtx = IPTools::absoluteImpactParameter3D(TrackpTT, bestVtxReFit);
                if (theDCAXVtx.first == false)
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 3D for track+" << std::endl;
                  continue;
                }
                double DCAKstTrkpVtx    = theDCAXVtx.second.value();
                double DCAKstTrkpVtxErr = theDCAXVtx.second.error();
  
  
                // #########################
                // # Compute B0 DCA to Vtx #
                // #########################
                theDCAXVtx = IPTools::absoluteImpactParameter3D(b_KP->refittedTransientTrack(), bestVtxReFit);
                if (theDCAXVtx.first == false)
                {
                  if (printMsg) std::cout << __LINE__ << " : continue --> invalid absolute impact parameter 3D for B0" << std::endl;
                  continue;
                }
                double DCAB0Vtx    = theDCAXVtx.second.value();
                double DCAB0VtxErr = theDCAXVtx.second.error();
  
  
                // ################################
                // # Compute B0 cos(alpha) to Vtx #
                // ################################
                Utility->computeCosAlpha (b_KP->currentState().globalMomentum().x(),b_KP->currentState().globalMomentum().y(),b_KP->currentState().globalMomentum().z(),
                                          b_KV->position().x() - bestVtxReFit.x(),b_KV->position().y() - bestVtxReFit.y(), b_KV->position().z() - bestVtxReFit.z(),
                                          b_KP->currentState().kinematicParametersError().matrix()(3,3),b_KP->currentState().kinematicParametersError().matrix()(4,4),b_KP->currentState().kinematicParametersError().matrix()(5,5),
                                          b_KP->currentState().kinematicParametersError().matrix()(3,4),b_KP->currentState().kinematicParametersError().matrix()(3,5),b_KP->currentState().kinematicParametersError().matrix()(4,5),
                                          b_KV->error().cxx() + bestVtxReFit.error()(0,0),b_KV->error().cyy() + bestVtxReFit.error()(1,1),b_KV->error().czz() + bestVtxReFit.error()(2,2),
                                          b_KV->error().matrix()(0,1) + bestVtxReFit.error()(0,1),b_KV->error().matrix()(0,2) + bestVtxReFit.error()(0,2),b_KV->error().matrix()(1,2) + bestVtxReFit.error()(1,2),
                                          &cosAlphaVtx,&cosAlphaVtxErr);
  
  
                // #############################
                // # Compute B0 L/sigma to Vtx #
                // #############################
                Utility->computeLS (b_KV->position().x(),b_KV->position().y(),b_KV->position().z(),
                                    bestVtxReFit.x(),bestVtxReFit.y(),bestVtxReFit.z(),
                                    b_KV->error().cxx(),b_KV->error().cyy(),b_KV->error().czz(),
                                    b_KV->error().matrix()(0,1),b_KV->error().matrix()(0,2),b_KV->error().matrix()(1,2),
                                    bestVtxReFit.error()(0,0),bestVtxReFit.error()(1,1),bestVtxReFit.error()(2,2),
                                    bestVtxReFit.error()(0,1),bestVtxReFit.error()(0,2),bestVtxReFit.error()(1,2),
                                    &LSVtx,&LSVtxErr);
  
  
                // #####################################################
                // # Compute ctau = mB * (BVtx - PVtx) dot pB / |pB|^2 #
                // #####################################################
                GlobalVector Bmomentum = b_KP->currentState().globalMomentum();
                double Bmomentum2 = Bmomentum.dot(Bmomentum);
                AlgebraicSymMatrix33 BmomentumErr2;
                BmomentumErr2(0,0) = b_KP->currentState().kinematicParametersError().matrix()(3,3);
                BmomentumErr2(0,1) = b_KP->currentState().kinematicParametersError().matrix()(3,4);
                BmomentumErr2(0,2) = b_KP->currentState().kinematicParametersError().matrix()(3,5);
                BmomentumErr2(1,0) = b_KP->currentState().kinematicParametersError().matrix()(4,3);
                BmomentumErr2(1,1) = b_KP->currentState().kinematicParametersError().matrix()(4,4);
                BmomentumErr2(1,2) = b_KP->currentState().kinematicParametersError().matrix()(4,5);
                BmomentumErr2(2,0) = b_KP->currentState().kinematicParametersError().matrix()(5,3);
                BmomentumErr2(2,1) = b_KP->currentState().kinematicParametersError().matrix()(5,4);
                BmomentumErr2(2,2) = b_KP->currentState().kinematicParametersError().matrix()(5,5);
  
                GlobalVector BVtxPVtxSep = GlobalPoint(b_KV->position()) - GlobalPoint(bestVtxReFit.position().x(), bestVtxReFit.position().y(), bestVtxReFit.position().z());
                double BVtxPVtxSepDOTBmomentum = BVtxPVtxSep.dot(Bmomentum);
  
                double bctauPVBS    = Utility->B0Mass * BVtxPVtxSepDOTBmomentum / Bmomentum2;
                double bctauPVBSErr = sqrt(Utility->B0Mass * Utility->B0Mass / (Bmomentum2 * Bmomentum2) *
                         
                                (Bmomentum.x() * Bmomentum.x() * bestVtxReFit.error()(0,0) +
                                 Bmomentum.y() * Bmomentum.y() * bestVtxReFit.error()(1,1) +
                                 Bmomentum.z() * Bmomentum.z() * bestVtxReFit.error()(2,2) +
                          
                                 ((BVtxPVtxSep.x()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.x()) *
                                  (BVtxPVtxSep.x()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.x()) *
                                  BmomentumErr2(0,0) +
                                  (BVtxPVtxSep.y()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.y()) *
                                  (BVtxPVtxSep.y()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.y()) *
                                  BmomentumErr2(1,1) +
                                  (BVtxPVtxSep.z()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.z()) *
                                  (BVtxPVtxSep.z()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.z()) *    
                                  BmomentumErr2(2,2) +
                          
                                  (BVtxPVtxSep.x()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.x()) * (BVtxPVtxSep.y()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.y()) * 2.*BmomentumErr2(0,1) +
                                  (BVtxPVtxSep.x()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.x()) * (BVtxPVtxSep.z()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.z()) * 2.*BmomentumErr2(0,2) +
                                  (BVtxPVtxSep.y()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.y()) * (BVtxPVtxSep.z()*Bmomentum2 - 2.*BVtxPVtxSepDOTBmomentum*Bmomentum.z()) * 2.*BmomentumErr2(1,2)) /
   
                                 (Bmomentum2 * Bmomentum2)));
  
  
  
  
                // #######################################
                // # @@@ Fill B0-candidate variables @@@ #
                // #######################################
                if (printMsg) std::cout << __LINE__ << " : @@@ Filling B0 candidate variables @@@\n\n" << std::endl;
  
  
                // ############
                // # Save: B0 #
                // ############
                NTuple->nB++;
  
                NTuple->bMass->push_back(b_KP->currentState().mass());
                NTuple->bMassE->push_back(sqrt(b_KP->currentState().kinematicParametersError().matrix()(6,6)));
                NTuple->bBarMass->push_back(bBar_KP->currentState().mass());
                NTuple->bBarMassE->push_back(sqrt(bBar_KP->currentState().kinematicParametersError().matrix()(6,6)));
  
                NTuple->bPx->push_back(b_KP->currentState().globalMomentum().x());
                NTuple->bPy->push_back(b_KP->currentState().globalMomentum().y());
                NTuple->bPz->push_back(b_KP->currentState().globalMomentum().z());          
  
                NTuple->bVtxCL->push_back(TMath::Prob(static_cast<double>(b_KV->chiSquared()), static_cast<int>(rint(b_KV->degreesOfFreedom()))));
                NTuple->bVtxX->push_back(b_KV->position().x());
                NTuple->bVtxY->push_back(b_KV->position().y());
                NTuple->bVtxZ->push_back(b_KV->position().z());
      
                NTuple->bCosAlphaVtx->push_back(cosAlphaVtx);
                NTuple->bCosAlphaVtxE->push_back(cosAlphaVtxErr);
                NTuple->bCosAlphaBS->push_back(cosAlphaBS);
                NTuple->bCosAlphaBSE->push_back(cosAlphaBSErr);
  
                NTuple->bLVtx->push_back(LSVtx);
                NTuple->bLVtxE->push_back(LSVtxErr);
                NTuple->bLBS->push_back(LSBS);
                NTuple->bLBSE->push_back(LSBSErr);
  
                NTuple->bDCAVtx->push_back(DCAB0Vtx);
                NTuple->bDCAVtxE->push_back(DCAB0VtxErr);
                NTuple->bDCABS->push_back(DCAB0BS);
                NTuple->bDCABSE->push_back(DCAB0BSErr);
   
                NTuple->bctauPVBS->push_back(bctauPVBS);
                NTuple->bctauPVBSE->push_back(bctauPVBSErr);
  
  
                // #############
                // # Save: K*0 #
                // #############
                NTuple->kstMass->push_back(kstInvMass);
                NTuple->kstMassE->push_back(sqrt(kst_KP->currentState().kinematicParametersError().matrix()(6,6)));
                NTuple->kstBarMass->push_back(kstBarInvMass);
                NTuple->kstBarMassE->push_back(sqrt(kstBar_KP->currentState().kinematicParametersError().matrix()(6,6)));
  
                NTuple->kstPx->push_back(kst_KP->currentState().globalMomentum().x());
                NTuple->kstPy->push_back(kst_KP->currentState().globalMomentum().y());
                NTuple->kstPz->push_back(kst_KP->currentState().globalMomentum().z());
  
                NTuple->kstVtxCL->push_back(TMath::Prob(static_cast<double>(kst_KV->chiSquared()), static_cast<int>(rint(kst_KV->degreesOfFreedom()))));
                NTuple->kstVtxX->push_back(kst_KV->position().x());
                NTuple->kstVtxY->push_back(kst_KV->position().y());
                NTuple->kstVtxZ->push_back(kst_KV->position().z());
  
  
                // #################
                // # Save: mu+ mu- #
                // #################
                NTuple->mumuMass->push_back(MuMuInvMass);
                NTuple->mumuMassE->push_back(sqrt(mumu_KP->currentState().kinematicParametersError().matrix()(6,6)));
  
                NTuple->mumuPx->push_back(mumu_KP->currentState().globalMomentum().x());
                NTuple->mumuPy->push_back(mumu_KP->currentState().globalMomentum().y());
                NTuple->mumuPz->push_back(mumu_KP->currentState().globalMomentum().z());
  
                NTuple->mumuVtxCL->push_back(TMath::Prob(static_cast<double>(mumu_KV->chiSquared()), static_cast<int>(rint(mumu_KV->degreesOfFreedom()))));
                NTuple->mumuVtxX->push_back(mumu_KV->position().x());
                NTuple->mumuVtxY->push_back(mumu_KV->position().y());
                NTuple->mumuVtxZ->push_back(mumu_KV->position().z());
  
                NTuple->mumuCosAlphaBS->push_back(MuMuCosAlphaBS);
                NTuple->mumuCosAlphaBSE->push_back(MuMuCosAlphaBSErr);
                NTuple->mumuLBS->push_back(MuMuLSBS);
                NTuple->mumuLBSE->push_back(MuMuLSBSErr);
                NTuple->mumuDCA->push_back(mumuDCA);
  
  
                // #############
                // # Save: mu- #
                // #############
                NTuple->mumHighPurity->push_back( (int)muTrackm->quality(reco::Track::highPurity));
                NTuple->mumCL->push_back(TMath::Prob(muTrackmTT.chi2(), static_cast<int>(rint(muTrackmTT.ndof()))));
                NTuple->mumNormChi2->push_back(muTrackm->normalizedChi2());
                NTuple->mumPx->push_back(refitMumTT.track().momentum().x());
                NTuple->mumPy->push_back(refitMumTT.track().momentum().y());
                NTuple->mumPz->push_back(refitMumTT.track().momentum().z());
  
                NTuple->mumDCAVtx->push_back(DCAmumVtx);
                NTuple->mumDCAVtxE->push_back(DCAmumVtxErr);
                NTuple->mumDCABS->push_back(DCAmumBS);
                NTuple->mumDCABSE->push_back(DCAmumBSErr);
  
                NTuple->mumKinkChi2->push_back(iMuonM->combinedQuality().trkKink);
                NTuple->mumFracHits->push_back(static_cast<double>(muTrackm->hitPattern().numberOfValidHits()) / static_cast<double>(muTrackm->hitPattern().numberOfValidHits() +
                                                                   muTrackm->hitPattern().numberOfLostHits(reco::HitPattern::TRACK_HITS) +
                                                                   muTrackm->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) +
                                                                   muTrackm->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS)));
                theDCAXVtx = IPTools::absoluteTransverseImpactParameter(muTrackmTT, bestVtxReFit);
                NTuple->mumdxyVtx->push_back(theDCAXVtx.second.value());
                NTuple->mumdzVtx->push_back(muTrackmTT.track().dz(bestVtxReFit.position()));
  
                NTuple->mumCat->push_back(getMuCat(*iMuonM));
  
                NTuple->mumNPixHits->push_back(muTrackm->hitPattern().numberOfValidPixelHits());
                NTuple->mumNPixLayers->push_back(muTrackm->hitPattern().pixelLayersWithMeasurement());  
                NTuple->mumNTrkHits->push_back(muTrackm->hitPattern().numberOfValidTrackerHits());
                NTuple->mumNTrkLayers->push_back(muTrackm->hitPattern().trackerLayersWithMeasurement());
                if (iMuonM->isGlobalMuon() == true) NTuple->mumNMuonHits->push_back(iMuonM->globalTrack()->hitPattern().numberOfValidMuonHits());
                else NTuple->mumNMuonHits->push_back(0);
                NTuple->mumNMatchStation->push_back(iMuonM->numberOfMatchedStations());
  
  
  
                // #############
                // # Save: mu+ #
                // #############
                NTuple->mupHighPurity->push_back( (int) muTrackp->quality(reco::Track::highPurity));
                NTuple->mupCL->push_back(TMath::Prob(muTrackpTT.chi2(), static_cast<int>(rint(muTrackpTT.ndof()))));
                NTuple->mupNormChi2->push_back(muTrackp->normalizedChi2());
                NTuple->mupPx->push_back(refitMupTT.track().momentum().x());
                NTuple->mupPy->push_back(refitMupTT.track().momentum().y());
                NTuple->mupPz->push_back(refitMupTT.track().momentum().z());
  
                NTuple->mupDCAVtx->push_back(DCAmupVtx);
                NTuple->mupDCAVtxE->push_back(DCAmupVtxErr);
                NTuple->mupDCABS->push_back(DCAmupBS);
                NTuple->mupDCABSE->push_back(DCAmupBSErr);
                
                NTuple->mupKinkChi2->push_back(iMuonP->combinedQuality().trkKink);
                NTuple->mupFracHits->push_back(static_cast<double>(muTrackp->hitPattern().numberOfValidHits()) / static_cast<double>(muTrackp->hitPattern().numberOfValidHits() +
                                                                   muTrackp->hitPattern().numberOfLostHits(reco::HitPattern::TRACK_HITS) +
                                                                   muTrackp->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS) +
                                                                   muTrackp->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_OUTER_HITS)));
                theDCAXVtx = IPTools::absoluteTransverseImpactParameter(muTrackpTT, bestVtxReFit);
                NTuple->mupdxyVtx->push_back(theDCAXVtx.second.value());
                NTuple->mupdzVtx->push_back(muTrackpTT.track().dz(bestVtxReFit.position()));
  
                NTuple->mupCat->push_back(getMuCat(*iMuonP));
  
                NTuple->mupNPixHits->push_back(muTrackp->hitPattern().numberOfValidPixelHits());
                NTuple->mupNPixLayers->push_back(muTrackp->hitPattern().pixelLayersWithMeasurement());  
                NTuple->mupNTrkHits->push_back(muTrackp->hitPattern().numberOfValidTrackerHits());
                NTuple->mupNTrkLayers->push_back(muTrackp->hitPattern().trackerLayersWithMeasurement());
                if (iMuonP->isGlobalMuon() == true) NTuple->mupNMuonHits->push_back(iMuonP->globalTrack()->hitPattern().numberOfValidMuonHits());
                else NTuple->mupNMuonHits->push_back(0);
                NTuple->mupNMatchStation->push_back(iMuonP->numberOfMatchedStations());
  
  
                // ################
                // # Save: Track- #
                // ################
                NTuple->kstTrkmHighPurity->push_back( (int)Trackm->quality(reco::Track::highPurity));
                NTuple->kstTrkmCL->push_back(TMath::Prob(TrackmTT.chi2(), static_cast<int>(rint(TrackmTT.ndof()))));
                NTuple->kstTrkmNormChi2->push_back(Trackm->normalizedChi2());
                NTuple->kstTrkmPx->push_back(refitTrkmTT.track().momentum().x());
                NTuple->kstTrkmPy->push_back(refitTrkmTT.track().momentum().y());
                NTuple->kstTrkmPz->push_back(refitTrkmTT.track().momentum().z());
  
                NTuple->kstTrkmDCAVtx->push_back(DCAKstTrkmVtx);
                NTuple->kstTrkmDCAVtxE->push_back(DCAKstTrkmVtxErr);
                NTuple->kstTrkmDCABS->push_back(DCAKstTrkmBS);
                NTuple->kstTrkmDCABSE->push_back(DCAKstTrkmBSErr);
  
                NTuple->kstTrkmFracHits->push_back(static_cast<double>(Trackm->hitPattern().numberOfValidHits()) / static_cast<double>(Trackm->hitPattern().numberOfValidHits() +
                                                                       Trackm->hitPattern().numberOfLostHits(reco::HitPattern::TRACK_HITS) +
                                                                       Trackm->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)));
                theDCAXVtx = IPTools::absoluteTransverseImpactParameter(TrackmTT, bestVtxReFit);
                NTuple->kstTrkmdxyVtx->push_back(theDCAXVtx.second.value());
                NTuple->kstTrkmdzVtx->push_back(TrackmTT.track().dz(bestVtxReFit.position()));
  
                NTuple->kstTrkmMuMatch->push_back(MuMCat);
  
                // I do NOT include the number of missing outer hits because the hadron might interact
                NTuple->kstTrkmNPixHits->push_back(Trackm->hitPattern().numberOfValidPixelHits());
                NTuple->kstTrkmNPixLayers->push_back(Trackm->hitPattern().pixelLayersWithMeasurement());  
                NTuple->kstTrkmNTrkHits->push_back(Trackm->hitPattern().numberOfValidTrackerHits());
                NTuple->kstTrkmNTrkLayers->push_back(Trackm->hitPattern().trackerLayersWithMeasurement());
  
  
                // ################
                // # Save: Track+ #
                // ################
                NTuple->kstTrkpHighPurity->push_back((int)Trackp->quality(reco::Track::highPurity));
                NTuple->kstTrkpCL->push_back(TMath::Prob(TrackpTT.chi2(), static_cast<int>(rint(TrackpTT.ndof()))));
                NTuple->kstTrkpNormChi2->push_back(Trackp->normalizedChi2());
                NTuple->kstTrkpPx->push_back(refitTrkpTT.track().momentum().x());
                NTuple->kstTrkpPy->push_back(refitTrkpTT.track().momentum().y());
                NTuple->kstTrkpPz->push_back(refitTrkpTT.track().momentum().z());
  
                NTuple->kstTrkpDCAVtx->push_back(DCAKstTrkpVtx);
                NTuple->kstTrkpDCAVtxE->push_back(DCAKstTrkpVtxErr);
                NTuple->kstTrkpDCABS->push_back(DCAKstTrkpBS);
                NTuple->kstTrkpDCABSE->push_back(DCAKstTrkpBSErr);
  
                NTuple->kstTrkpFracHits->push_back(static_cast<double>(Trackp->hitPattern().numberOfValidHits()) / static_cast<double>(Trackp->hitPattern().numberOfValidHits() +
                                                                       Trackp->hitPattern().numberOfLostHits(reco::HitPattern::TRACK_HITS) +
                                                                       Trackp->hitPattern().numberOfLostHits(reco::HitPattern::MISSING_INNER_HITS)));
                theDCAXVtx = IPTools::absoluteTransverseImpactParameter(TrackpTT, bestVtxReFit);
                NTuple->kstTrkpdxyVtx->push_back(theDCAXVtx.second.value());
                NTuple->kstTrkpdzVtx->push_back(TrackpTT.track().dz(bestVtxReFit.position()));
  
                NTuple->kstTrkpMuMatch->push_back(MuPCat);
  
                // I do NOT include the number of missing outer hits because the hadron might interact
                NTuple->kstTrkpNPixHits->push_back(Trackp->hitPattern().numberOfValidPixelHits());
                NTuple->kstTrkpNPixLayers->push_back(Trackp->hitPattern().pixelLayersWithMeasurement());  
                NTuple->kstTrkpNTrkHits->push_back(Trackp->hitPattern().numberOfValidTrackerHits());
                NTuple->kstTrkpNTrkLayers->push_back(Trackp->hitPattern().trackerLayersWithMeasurement());
  
  
                // Save trigger matching for the 4 tracks
                tmpString1.clear(); tmpString2.clear(); tmpString3.clear(); tmpString4.clear();
                const pat::Muon* patMuonM = &(*iMuonM);
                const pat::Muon* patMuonP = &(*iMuonP);
                const pat::GenericParticle* patTrkm = &(*iTrackM);
                const pat::GenericParticle* patTrkp = &(*iTrackP);
                for (unsigned int i = 0; i < TrigTable_.size(); i++)
                  {
                    myString.clear(); myString.str(""); myString << TrigTable_[i].c_str() << "*";
                    if (patMuonM->triggerObjectMatchesByPath(myString.str().c_str()).empty() == false) tmpString1.append(TrigTable_[i]+" ");
                    if (patMuonP->triggerObjectMatchesByPath(myString.str().c_str()).empty() == false) tmpString2.append(TrigTable_[i]+" ");
                    if (patTrkm ->triggerObjectMatchesByPath(myString.str().c_str()).empty() == false) tmpString3.append(TrigTable_[i]+" ");
                    if (patTrkp ->triggerObjectMatchesByPath(myString.str().c_str()).empty() == false) tmpString4.append(TrigTable_[i]+" ");
                  }
                if (tmpString1.size() == 0) tmpString1.append("NotInTable");
                if (tmpString2.size() == 0) tmpString2.append("NotInTable");
                if (tmpString3.size() == 0) tmpString3.append("NotInTable");
                if (tmpString4.size() == 0) tmpString4.append("NotInTable");
                NTuple->mumTrig    ->push_back(tmpString1);
                NTuple->mupTrig    ->push_back(tmpString2);
                NTuple->kstTrkmTrig->push_back(tmpString3);
                NTuple->kstTrkpTrig->push_back(tmpString4);
  
                // save minimum IP from any PV (2D and 3D)
                mumMind0  = mupMind0 = 100;
                TrkmMind0 = 100;
                TrkpMind0 = 100;
                mumMinIP  = 100;
                mupMinIP  = 100;
                TrkmMinIP = 100;
                TrkpMinIP = 100;
                std::pair<double,double>  IPPair;

                for (std::vector<reco::Vertex>::const_iterator ipv = recVtx->begin(); ipv != recVtx->end(); ipv++) { 
                    vert = GlobalPoint(ipv->x(), ipv->y(), ipv->z());

                    traj = muTrackmTT.trajectoryStateClosestToPoint(vert );
                    if (fabs(traj.perigeeParameters().transverseImpactParameter()) < mumMind0){
                      mumMind0  = fabs(traj.perigeeParameters().transverseImpactParameter());
                      mumMind0E = traj.perigeeError().transverseImpactParameterError();
                    }  

                    traj = muTrackpTT.trajectoryStateClosestToPoint(vert );
                    if (fabs(traj.perigeeParameters().transverseImpactParameter()) < mupMind0){
                      mupMind0  = fabs(traj.perigeeParameters().transverseImpactParameter());
                      mupMind0E = traj.perigeeError().transverseImpactParameterError();
                    }  

                    traj = TrackmTT.trajectoryStateClosestToPoint(vert );
                    if (fabs(traj.perigeeParameters().transverseImpactParameter()) < TrkmMind0){
                      TrkmMind0  = fabs(traj.perigeeParameters().transverseImpactParameter());
                      TrkmMind0E = traj.perigeeError().transverseImpactParameterError();
                    }  

                    traj = TrackpTT.trajectoryStateClosestToPoint(vert );
                    if (fabs(traj.perigeeParameters().transverseImpactParameter()) < TrkpMind0){
                      TrkpMind0  = fabs(traj.perigeeParameters().transverseImpactParameter());
                      TrkpMind0E = traj.perigeeError().transverseImpactParameterError();
                    }  
                    
                    IPPair = pionImpactParameter(muTrackmTT,*ipv);
                    if (IPPair.first < mumMinIP){
                      mumMinIP  = IPPair.first;
                      mumMinIPE = IPPair.second;
                    }
                    IPPair = pionImpactParameter(muTrackpTT,*ipv);
                    if (IPPair.first < mupMinIP){
                      mupMinIP  = IPPair.first;
                      mupMinIPE = IPPair.second;
                    }
                    IPPair = pionImpactParameter(TrackmTT,*ipv);
                    if (IPPair.first < TrkmMinIP){
                      TrkmMinIP  = IPPair.first;
                      TrkmMinIPE = IPPair.second;
                    }
                    IPPair = pionImpactParameter(TrackpTT,*ipv);
                    if (IPPair.first < TrkpMinIP){
                      TrkpMinIP  = IPPair.first;
                      TrkpMinIPE = IPPair.second;
                    }
                }

                // isolation
                double iso;
                chi = 0; ndf = 0;
                mum_isovec.clear(); mup_isovec.clear(); trkm_isovec.clear(); trkp_isovec.clear();
                float bestVtx_x = bestVtx.x();
                float bestVtx_y = bestVtx.y();
                float bestVtx_z = bestVtx.z();

                 for (std::vector<pat::GenericParticle>::const_iterator iTrackIso = thePATTrackHandle->begin(); iTrackIso != thePATTrackHandle->end(); iTrackIso++)
                 {
                   if (iTrackIso == iTrackM || iTrackIso == iTrackP) continue;
                   if (muTrackm == iTrackIso->track() || muTrackp == iTrackIso->track()) continue;
                  const reco::TransientTrack TrackIsoTT(iTrackIso->track(), &(*bFieldHandle));
                  
                  iso = calculateIsolation(ClosestApp, muTrackmTT, TrackIsoTT, muonMass, bestVtx_x, bestVtx_y, bestVtx_z);
                  if (iso > 0 ) mum_isovec.push_back(iso);

                  iso = calculateIsolation(ClosestApp, muTrackpTT, TrackIsoTT, muonMass, bestVtx_x, bestVtx_y, bestVtx_z);
                  if (iso > 0 ) mup_isovec.push_back(iso);

                  iso = calculateIsolation(ClosestApp, TrackmTT, TrackIsoTT, muonMass, bestVtx_x, bestVtx_y, bestVtx_z);
                  if (iso > 0 ) trkm_isovec.push_back(iso);

                  iso = calculateIsolation(ClosestApp, TrackpTT, TrackIsoTT, muonMass, bestVtx_x, bestVtx_y, bestVtx_z);
                  if (iso > 0 ) trkp_isovec.push_back(iso);
                 } 

                NTuple->mumIso      -> push_back(mum_isovec);
                NTuple->mupIso      -> push_back(mup_isovec);
                NTuple->kstTrkmIso  -> push_back(trkm_isovec);
                NTuple->kstTrkpIso  -> push_back(trkp_isovec);


                NTuple->mumMinIP2D      -> push_back(mumMind0);
                NTuple->mumMinIP2DE     -> push_back(mumMind0E);
                NTuple->mupMinIP2D      -> push_back(mupMind0);
                NTuple->mupMinIP2DE     -> push_back(mupMind0E);
                NTuple->kstTrkmMinIP2D  -> push_back(TrkmMind0);
                NTuple->kstTrkmMinIP2DE -> push_back(TrkmMind0E);
                NTuple->kstTrkpMinIP2D  -> push_back(TrkpMind0);
                NTuple->kstTrkpMinIP2DE -> push_back(TrkpMind0E);

                NTuple->mumMinIP      -> push_back(mumMinIP);
                NTuple->mumMinIPE     -> push_back(mumMinIPE);
                NTuple->mupMinIP      -> push_back(mupMinIP);
                NTuple->mupMinIPE     -> push_back(mupMinIPE);
                NTuple->kstTrkmMinIP  -> push_back(TrkmMinIP);
                NTuple->kstTrkmMinIPE -> push_back(TrkmMinIPE);
                NTuple->kstTrkpMinIP  -> push_back(TrkpMinIP);
                NTuple->kstTrkpMinIPE -> push_back(TrkpMinIPE);
 
  
                // #####################
                // # Clear all vectors #
                // #####################
                vertexTracks.clear();
                bParticles.clear();
                bBarParticles.clear();
                kstParticles.clear();
                kstBarParticles.clear();
              } // End for Track+
              if (skip == true) continue;
            } // End for Track-
          muonParticles.clear(); 
        } // End for mu+
        if (skip == true) continue;
      } // End for mu-
    } // End if bestVtx is true
    else if (printMsg) std::cout << __LINE__ << " : continue --> invalid Pri.Vtx" << std::endl;
//     if (B0MassConstraint != NULL) delete B0MassConstraint;


    NTuple->runN   = iEvent.id().run();
    NTuple->eventN = iEvent.id().event();
    NTuple->bsX    = beamSpot.position().x();
    NTuple->bsY    = beamSpot.position().y();  
    
    if (NTuple->nB > 0)
    {
      for (std::vector<reco::Vertex>::const_iterator iVertex = recVtx->begin(); iVertex != recVtx->end(); iVertex++)
        {
          if(iVertex->ndof() < PRIVTXNDOF)                 continue;
          if(fabs(iVertex->z()) > PRIVTXMAXZ)              continue;
          if(fabs(iVertex->position().rho()) > PRIVTXMAXR) continue;
          NTuple->recoVtxN++;
        }

      // #####################################
      // # Save: Primary Vertex and BeamSpot #
      // #####################################
      NTuple->priVtxCL = TMath::Prob(bestVtx.chi2(), static_cast<int>(rint(bestVtx.ndof())));
      NTuple->priVtxX  = bestVtx.x();
      NTuple->priVtxY  = bestVtx.y();
      NTuple->priVtxZ  = bestVtx.z();
    }
    else if (printMsg) std::cout << __LINE__ << " : @@@ No B0 --> K*0 (K pi) mu+ mu- candidate were found in the event @@@" << std::endl;
  } // End if doGenReco_ == 1 || doGenReco_ == 2






  // ###########################################
  // # Get truth information from genParticles #
  // ###########################################


  // Particle IDs:
  //   11 = e-
  //   12 = nu_e
  //   13 = mu-
  //   14 = nu_mu
  //   16 = nu_tau
  //   22 = gamma
  //  130 = KL
  //  211 = pi+
  //  310 = Ks
  //  311 = K0
  //  313 = K*0
  //  321 = K+
  //  443 = J/psi; 100443 psi(2S)
  //  511 = B0
  //  521 = B+
  //  531 = Bs

  // 2212 = p
  // 5122 = Lambda_b


  if (doGenReco_ == 2)
  {
    // #################################
    // # Save pileup information in MC #
    // #################################
    edm::Handle< std::vector<PileupSummaryInfo> > PupInfo;
    iEvent.getByToken(puToken_, PupInfo);
    for (std::vector<PileupSummaryInfo>::const_iterator PVI = PupInfo->begin(); PVI != PupInfo->end(); PVI++)
    {
      NTuple->bunchXingMC->push_back(PVI->getBunchCrossing());
      NTuple->numInteractionsMC->push_back(PVI->getPU_NumInteractions());
      NTuple->trueNumInteractionsMC->push_back(PVI->getTrueNumInteractions());
    }
  }


  if ((doGenReco_ == 2) || (doGenReco_ == 3))
  {
    std::vector<std::vector<unsigned int>* > posDau;
    std::vector<std::vector<unsigned int>* > negDau;

    for (unsigned int itGen = 0; itGen < genParticles->size(); itGen++)
    {
      const reco::Candidate* FirstPart = &(*genParticles)[itGen];


      // ##########################
      // # Check for:             #
      // # B0 / B0bar             #
      // # Bs / Bsbar             #
      // # Lambda_b / Lambda_bbar #
      // # B+ / B-                #
      // ##########################
      if ((abs(FirstPart->pdgId()) == 511) || (abs(FirstPart->pdgId()) == 531) || (abs(FirstPart->pdgId()) == 5122) || (abs(FirstPart->pdgId()) == 521))
      {         
        if (printMsg) std::cout << "\n" << __LINE__ << " : @@@ Found B0/B0bar OR Bs/Bsbar OR Lambda_b/Lambda_bbar OR B+/B- in MC @@@" << std::endl;

        // ########################################################
        // # Search for oscillations                              #
        // # If there are, then move to the non oscillating stage #
        // ########################################################
        FirstPart = skipOscillations(FirstPart);

        int imum      = -1;
        int imup      = -1;
        int ikst      = -1;
        int ikst_trkm = -1;
        int ikst_trkp = -1;
        int ipsi      = -1;
        int ipsi_mum  = -1;
        int ipsi_mup  = -1;

        bool PhotonB0 = false;
        bool PhotonKst = false;
        bool PhotonPsi = false;

        bool isSignal = true;
        bool isPsi2SnotJPsi = false;

        if (abs(FirstPart->pdgId()) == 511)
        {    
          for (unsigned int i = 0; i < FirstPart->numberOfDaughters(); i++)
          {
            if (isSignal == false) break;
              
            const reco::Candidate* genDau = FirstPart->daughter(i);

            if (genDau->pdgId() == 22)
            {
              PhotonB0 = true;
              if (printMsg) std::cout << __LINE__ << " : found B0/B0bar photon" << std::endl;
              continue;
            }
            if (genDau->pdgId() == 13)
            {
              imum = i;
              if (printMsg) std::cout << __LINE__ << " : found mu-" << std::endl;
              continue;
            }
            if (genDau->pdgId() == -13)
            {
              imup = i;
              if (printMsg) std::cout << __LINE__ << " : found mu+" << std::endl;
              continue;
            }
            if (abs(genDau->pdgId()) == 313)
            {
              ikst = i;

              if (printMsg) std::cout << __LINE__ << " : found K*0/K*0bar" << std::endl;

              for (unsigned int j = 0; j < genDau->numberOfDaughters(); j++)
              {
                const reco::Candidate* genDau2 = genDau->daughter(j);
              
                if (genDau2->pdgId() == 22)
                {
                  PhotonKst = true;
                  if (printMsg) std::cout << __LINE__ << " : found K*0/K*0bar photon" << std::endl;
                  continue;
                }
                if ((genDau2->pdgId() == -211) || (genDau2->pdgId() == -321))
                {
                  ikst_trkm = j;
                  if (printMsg) std::cout << __LINE__ << " : found K*0/K*0bar track-" << std::endl;
                  continue;
                }
                if ((genDau2->pdgId() == 211) || (genDau2->pdgId() == 321))
                {
                  ikst_trkp = j;
                  if (printMsg) std::cout << __LINE__ << " : found K*0/K*0bar track+" << std::endl;
                  continue;
                }
                isSignal = false;
                break;
              }
              continue;
            }
            if (abs(genDau->pdgId() == 443) || abs(genDau->pdgId() == 100443))
            {
              if (abs(genDau->pdgId() == 443)) isPsi2SnotJPsi = false;
              else isPsi2SnotJPsi = true;

              ipsi = i;
              
              if (printMsg) std::cout << __LINE__ << " : found J/psi or psi(2S)" << std::endl;

              for (unsigned int j = 0; j < genDau->numberOfDaughters(); j++)
              {
                const reco::Candidate* genDau2 = genDau->daughter(j);
              
                if (genDau2->pdgId() == 22)
                {
                  PhotonPsi = true;
                  if (printMsg) std::cout << __LINE__ << " : found J/psi or psi(2S) photon" << std::endl;
                  continue;
                }
                if (genDau2->pdgId() == 13)
                {
                  ipsi_mum = j;
                  if (printMsg) std::cout << __LINE__ << " : found J/psi or psi(2S) mu-" << std::endl;
                  continue;
                }
                if (genDau2->pdgId() == -13)
                {
                  ipsi_mup = j;
                  if (printMsg) std::cout << __LINE__ << " : found J/psi or psi(2S) mu+" << std::endl;
                  continue;
                }
                isSignal = false;
                break;
              }
            continue;
            }
            isSignal = false;
            break;
          } // End for B0 / B0bar daughters
        } // End if B0 / B0bar


        const reco::Candidate* genPsi = NULL;
        const reco::Candidate* genMum = NULL;
        const reco::Candidate* genMup = NULL;
        const reco::Candidate* genKst = NULL;
        const reco::Candidate* genKst_trkm = NULL;
        const reco::Candidate* genKst_trkp = NULL;


        bool foundSomething = false;
        if ((abs(FirstPart->pdgId()) == 511) && (isSignal == true) &&
           (((imum != -1) && (imup != -1) && (ikst != -1) && (ikst_trkm != -1) && (ikst_trkp != -1)) ||
           ((NTuple->genSignal != 1) && (NTuple->genSignal != 2) &&
            (ikst != -1) && (ikst_trkm != -1) && (ikst_trkp != -1) && (ipsi != -1) && (ipsi_mum != -1) && (ipsi_mup != -1))))
        {
          // ################################################################
          // # Found Signal B0/B0bar --> K*0/K*0bar (K pi) mu+ mu-          #
          // # Found Signal B0/B0bar --> K*0/K*0bar (K pi) J/psi (mu+mu-)   #
          // # Found Signal B0/B0bar --> K*0/K*0bar (K pi) psi(2S) (mu+mu-) #
          // ################################################################

          NTuple->ClearMonteCarlo();

          if ((imum != -1) && (imup != -1) && (ikst != -1) && (ikst_trkm != -1) && (ikst_trkp != -1))
          {
            if (printMsg) std::cout << __LINE__ << " : @@@ Found B0/B0bar --> K*0/K*0bar (K pi) mu+ mu- @@@" << std::endl;
              
            if      (FirstPart->daughter(ikst)->pdgId() == 313)  NTuple->genSignal = 1;
            else if (FirstPart->daughter(ikst)->pdgId() == -313) NTuple->genSignal = 2;
              
            genMum = FirstPart->daughter(imum);
            genMup = FirstPart->daughter(imup);
          }
          else
          {
            if (isPsi2SnotJPsi == false)
            {
              if (printMsg) std::cout << __LINE__ << " : @@@ Found B0/B0bar --> K*0/K*0bar (K pi) J/psi (mu+mu-) @@@" << std::endl;
              
              if      (FirstPart->daughter(ikst)->pdgId() == 313)  NTuple->genSignal = 3;
              else if (FirstPart->daughter(ikst)->pdgId() == -313) NTuple->genSignal = 4;
            }
            else
            {
              if (printMsg) std::cout << __LINE__ << " : @@@ Found B0/B0bar --> K*0/K*0bar (K pi) psi(2S) (mu+mu-) @@@" << std::endl;
              
              if      (FirstPart->daughter(ikst)->pdgId() == 313)  NTuple->genSignal = 5;
              else if (FirstPart->daughter(ikst)->pdgId() == -313) NTuple->genSignal = 6;
            }
              
            genPsi = FirstPart->daughter(ipsi);
            genMum = genPsi->daughter(ipsi_mum);
            genMup = genPsi->daughter(ipsi_mup);

            NTuple->genSignPsiHasFSR = PhotonPsi;
          }

          genKst = FirstPart->daughter(ikst);
          genKst_trkm = FirstPart->daughter(ikst)->daughter(ikst_trkm);
          genKst_trkp = FirstPart->daughter(ikst)->daughter(ikst_trkp);

 
          // ##################
          // # Save info. FSR #
          // ##################
          NTuple->genSignHasFSR = PhotonB0;
          NTuple->genSignKstHasFSR = PhotonKst;


          foundSomething = true;
        } // End if B0/B0bar signal
        else
        {
          if (printMsg) std::cout << __LINE__ << " : @@@ Start particle decay-tree scan for opposite charged stable particles for background search @@@" << std::endl;
          searchForStableChargedDaughters (FirstPart, itGen, &posDau, &negDau);


          if ((NTuple->genSignal == 0) && (posDau.size() != 0) && (negDau.size() != 0))
          {
            // ###############################
            // # Search for Background from: #
            // # B0 / B0bar                  #
            // # Bs / Bsbar                  #
            // # Lambda_b / Lambda_bbar      #
            // # B+ / B-                     #
            // ###############################

            for (unsigned int i = 0; i < negDau.size(); i++)
            {
              genMum = findDaughter (genParticles, &negDau, i);
              if (genMum->pdgId() == 13)
              {
                for (unsigned int j = 0; j < posDau.size(); j++)
                {
                  genMup = findDaughter (genParticles, &posDau, j);
                  if (genMup->pdgId() == -13)
                  {
                    if (printMsg) std::cout << __LINE__ << " : found dimuons for possible background" << std::endl;
                    foundSomething = true;
                    break;
                  }
                }
                if (foundSomething == true) break;
              }
            }

            if ((foundSomething == true) && (posDau.size() == 2) && (negDau.size() == 1))
            {
              foundSomething = false;
              for (unsigned int i = 0; i < posDau.size(); i++)
              {
                genKst_trkp = findDaughter (genParticles, &posDau, i);
                if (genKst_trkp != genMup)
                {
                  NTuple->ClearMonteCarlo();

                  NTuple->genMuMuBG = FirstPart->pdgId();
                  NTuple->genMuMuBGnTrksp = 1;
                  foundSomething = true;
                  if (printMsg) std::cout << __LINE__ << " : get background positive track: " << genKst_trkp->pdgId() << "\tfrom mother: " << genKst_trkp->mother()->pdgId() << std::endl;
                }
              }
            }
            else if ((foundSomething == true) && (posDau.size() == 1) && (negDau.size() == 2))
            {
              foundSomething = false;
              for (unsigned int i = 0; i < negDau.size(); i++)
              {
                genKst_trkm = findDaughter (genParticles, &negDau, i);
                if (genKst_trkm != genMum)
                {
                  NTuple->ClearMonteCarlo();

                  NTuple->genMuMuBG = FirstPart->pdgId();              
                  NTuple->genMuMuBGnTrksm = 1;
                  foundSomething = true;
                  if (printMsg) std::cout << __LINE__ << " : get background negative track: " << genKst_trkm->pdgId() << "\tfrom mother: " << genKst_trkm->mother()->pdgId() << std::endl;
                }
              }
            }
            else if ((foundSomething == true) && (posDau.size() == 2) && (negDau.size() == 2))
            {
              foundSomething = false;
              for (unsigned int i = 0; i < negDau.size(); i++)
              {
                genKst_trkm = findDaughter (genParticles, &negDau, i);
                if (genKst_trkm != genMum)
                {
                  for (unsigned int j = 0; j < posDau.size(); j++)
                  {
                    genKst_trkp = findDaughter (genParticles, &posDau, j);
                    if (genKst_trkp != genMup)
                    {
                      NTuple->ClearMonteCarlo();

                      NTuple->genMuMuBG = FirstPart->pdgId();
                      NTuple->genMuMuBGnTrksp = 1;
                      NTuple->genMuMuBGnTrksm = 1;
                      foundSomething = true;
                      if (printMsg)
                      {
                        std::cout << __LINE__ << " : get background negative track: " << genKst_trkm->pdgId() << "\tfrom mother: " << genKst_trkm->mother()->pdgId();
                        std::cout << "\tand positive track: " << genKst_trkp->pdgId() << "\tfrom mother: " << genKst_trkp->mother()->pdgId() << std::endl;
                      }
                    }
                  }
                }
              }
            }
            else if ((foundSomething == true) && (posDau.size() >= 2) && (negDau.size() >= 2))
            {
              foundSomething = false;

              double bestMass = 0.;
              unsigned int negDauIndx = 0;
              unsigned int posDauIndx = 0;

              for (unsigned int i = 0; i < negDau.size(); i++)
              {
                genKst_trkm = findDaughter (genParticles, &negDau, i);
                if (genKst_trkm != genMum)
                {
                  for (unsigned int j = 0; j < posDau.size(); j++)
                  {
                    genKst_trkp = findDaughter (genParticles, &posDau, j);
                    if (genKst_trkp != genMup)
                    {
                      double invMass = Utility->computeInvMass (genKst_trkm->px(),genKst_trkm->py(),genKst_trkm->pz(),genKst_trkm->mass(),
                                                                genKst_trkp->px(),genKst_trkp->py(),genKst_trkp->pz(),genKst_trkp->mass());
                      
                      if (fabs(invMass - Utility->kstMass) < fabs(bestMass - Utility->kstMass))
                      {
                        bestMass = invMass;
                        negDauIndx = i;
                        posDauIndx = j;
                      }
                    }
                  }
                }
              }

              if (bestMass > 0.)
              {
                genKst_trkm = findDaughter (genParticles, &negDau, negDauIndx);
                genKst_trkp = findDaughter (genParticles, &posDau, posDauIndx);

                NTuple->ClearMonteCarlo();
                NTuple->genMuMuBG = FirstPart->pdgId();
                NTuple->genMuMuBGnTrksp = 1;
                NTuple->genMuMuBGnTrksm = 1;
                foundSomething = true;
                if (printMsg)
                {
                  std::cout << __LINE__ << " : get background negative track: " << genKst_trkm->pdgId() << "\tfrom mother: " << genKst_trkm->mother()->pdgId();
                  std::cout << "\tand positive track: " << genKst_trkp->pdgId() << "\tfrom mother: " << genKst_trkp->mother()->pdgId() << std::endl;
                }
              }
            }
            else
            {
              foundSomething = false;
              if (printMsg) std::cout << __LINE__ << " : @@@ No background found @@@" << std::endl;
            }
          }
          else if (printMsg) std::cout << __LINE__ << " : @@@ No possible background found @@@" << std::endl;
        } // End else background


        // ########################
        // # Save gen information #
        // ########################
        if (foundSomething == true)
        {
          if (printMsg) std::cout << __LINE__ << " : @@@ Saving signal OR background compatible with signal @@@" << std::endl;


          // #############################
          // # Search for primary vertex #
          // #############################
          const reco::Candidate* PVtx = FirstPart;
          while (PVtx->numberOfMothers() > 0) PVtx = PVtx->mother(0);


          // ############################
          // # Generated Primary Vertex #
          // ############################
          NTuple->genPriVtxX = PVtx->vx();
          NTuple->genPriVtxY = PVtx->vy();
          NTuple->genPriVtxZ = PVtx->vz();

          
          // #####################
          // # Generated B0 Mass #
          // #####################
          NTuple->genB0Mass = FirstPart->mass();
          NTuple->genB0Px = FirstPart->px();
          NTuple->genB0Py = FirstPart->py();
          NTuple->genB0Pz = FirstPart->pz();


          // ####################
          // # Generated B0 Vtx #
          // ####################
          NTuple->genB0VtxX = FirstPart->vx();
          NTuple->genB0VtxY = FirstPart->vy();
          NTuple->genB0VtxZ = FirstPart->vz();
          
    
          if (NTuple->genSignal != 0)
          {
            // ######################
            // # Generated K*0 Mass #
            // ######################
            NTuple->genKstMass = genKst->mass();
            NTuple->genKstPx = genKst->px();
            NTuple->genKstPy = genKst->py();
            NTuple->genKstPz = genKst->pz();
              
            // #####################
            // # Generated K*0 Vtx #
            // #####################
            NTuple->genKstVtxX = genKst->vx();
            NTuple->genKstVtxY = genKst->vy();
            NTuple->genKstVtxZ = genKst->vz();
          }
          else if ((NTuple->genMuMuBGnTrksm != 0) && (NTuple->genMuMuBGnTrksp != 0) &&
               (genKst_trkm->vx() == genKst_trkp->vx()) &&
               (genKst_trkm->vy() == genKst_trkp->vy()) &&
               (genKst_trkm->vz() == genKst_trkp->vz()))
          {
            double invMass = Utility->computeInvMass (genKst_trkm->momentum().x(),genKst_trkm->momentum().y(),genKst_trkm->momentum().z(),genKst_trkm->mass(),
                                                      genKst_trkp->momentum().x(),genKst_trkp->momentum().y(),genKst_trkp->momentum().z(),genKst_trkp->mass());

            // ######################
            // # Generated K*0 Mass #
            // ######################
            NTuple->genKstMass = invMass;
            NTuple->genKstPx = genKst_trkm->momentum().x() + genKst_trkp->momentum().x();
            NTuple->genKstPy = genKst_trkm->momentum().y() + genKst_trkp->momentum().y();
            NTuple->genKstPz = genKst_trkm->momentum().z() + genKst_trkp->momentum().z();

            // #####################
            // # Generated K*0 Vtx #
            // #####################
            NTuple->genKstVtxX = genKst_trkm->vx();
            NTuple->genKstVtxY = genKst_trkm->vy();
            NTuple->genKstVtxZ = genKst_trkm->vz();
          }
       
          // ###########################################
          // # Generated J/psi or psi(2S) Mass and Vtx #
          // ###########################################
          if ((NTuple->genSignal == 3) || (NTuple->genSignal == 4) || (NTuple->genSignal == 5) || (NTuple->genSignal == 6))
          {
            NTuple->genPsiMass = genPsi->mass();
            NTuple->genPsiVtxX = genPsi->vx();
            NTuple->genPsiVtxY = genPsi->vy();
            NTuple->genPsiVtxZ = genPsi->vz();
          }
          else if ((NTuple->genSignal == 0) &&
               (genMum->vx() == genMup->vx()) &&
               (genMum->vy() == genMup->vy()) &&
               (genMum->vz() == genMup->vz()))
          {
            double invMass = Utility->computeInvMass (genMum->momentum().x(),genMum->momentum().y(),genMum->momentum().z(),genMum->mass(),
                                                      genMup->momentum().x(),genMup->momentum().y(),genMup->momentum().z(),genMup->mass());
              

            NTuple->genPsiMass = invMass;
            NTuple->genPsiVtxX = genMum->vx();
            NTuple->genPsiVtxY = genMum->vy();
            NTuple->genPsiVtxZ = genMum->vz();
          }

          // #################
          // # Generated mu- #
          // #################
          NTuple->genMumMother = genMum->mother()->pdgId();
          NTuple->genMumPx = genMum->px();
          NTuple->genMumPy = genMum->py();
          NTuple->genMumPz = genMum->pz();
          
          // #################
          // # Generated mu+ #
          // #################
          NTuple->genMupMother = genMup->mother()->pdgId();
          NTuple->genMupPx = genMup->px();
          NTuple->genMupPy = genMup->py();
          NTuple->genMupPz = genMup->pz();

          // ########################
          // # Generated K*0 track- #
          // ########################
          if ((NTuple->genSignal != 0) || (NTuple->genMuMuBGnTrksm != 0))
          {
            NTuple->genKstTrkmMother = genKst_trkm->mother()->pdgId();
            NTuple->genKstTrkmID = genKst_trkm->pdgId();
            NTuple->genKstTrkmPx = genKst_trkm->px();
            NTuple->genKstTrkmPy = genKst_trkm->py();
            NTuple->genKstTrkmPz = genKst_trkm->pz();
          }

          // ########################
          // # Generated K*0 track+ #
          // ########################
          if ((NTuple->genSignal != 0) || (NTuple->genMuMuBGnTrksp != 0))
          {
            NTuple->genKstTrkpMother = genKst_trkp->mother()->pdgId();
            NTuple->genKstTrkpID = genKst_trkp->pdgId();
            NTuple->genKstTrkpPx = genKst_trkp->px();
            NTuple->genKstTrkpPy = genKst_trkp->py();
            NTuple->genKstTrkpPz = genKst_trkp->pz();
          }
        } // End if foundSomething
      } // End if B0/B0bar OR Bs/Bsbar OR Lambda_b/Lambda_bbar OR B+/B-


      // ##############################################
      // # Check to see if J/psi or psi(2S) is prompt #
      // ##############################################
      bool isPrompt = false;
      const reco::Candidate& PsiCand = (*genParticles)[itGen];
      
      if ((abs(PsiCand.pdgId()) == 443) || (abs(PsiCand.pdgId()) == 100443))
      {
        isPrompt = true;

        for (unsigned int i = 0; i < PsiCand.numberOfMothers(); i++)
        {
          const reco::Candidate* psiMom = PsiCand.mother(i);
          
          if (((abs(psiMom->pdgId()) < 600) && (abs(psiMom->pdgId()) > 500)) || ((abs(psiMom->pdgId()) < 6000) && (abs(psiMom->pdgId()) > 5000)))
          {
            isPrompt = false;
            continue;
          }
          else
          {
            for (unsigned int i = 0; i < psiMom->numberOfMothers(); i++)
            {
              const reco::Candidate* psiMom2 = psiMom->mother(i);

              if (((abs(psiMom2->pdgId()) < 600) && (abs(psiMom2->pdgId()) > 500)) || ((abs(psiMom2->pdgId()) < 6000) && (abs(psiMom2->pdgId()) > 5000)))
              {
                isPrompt = false;
                continue;
              }
              else
              {
                for (unsigned int i = 0; i < psiMom2->numberOfMothers(); i++)
                {
                  const reco::Candidate* psiMom3 = psiMom2->mother(i);

                  if (((abs(psiMom3->pdgId()) < 600) && (abs(psiMom3->pdgId()) > 500)) || ((abs(psiMom3->pdgId()) < 6000) && (abs(psiMom3->pdgId()) > 5000)))
                  {
                    isPrompt = false;
                    continue;
                  }
                }
              }
            }
          }
        }
          
        if (isPrompt == true)
        {
          NTuple->genPsiPrompt = 1;
          if (printMsg) std::cout << __LINE__ << " : found prompt J/psi or psi(2S)" << std::endl;
        }
        continue;
      }

      for (unsigned int i = 0; i < posDau.size(); i++)
      {
        posDau[i]->clear();
        delete posDau[i];
      }
      posDau.clear();
      for (unsigned int i = 0; i < negDau.size(); i++)
      {
        negDau[i]->clear();
        delete negDau[i];
      }
      negDau.clear();
    } // End for genParticles
  } // End if doGenReco_ == 2 || doGenReco_ == 3




  // ####################################
  // # Perform matching with candidates #
  // ####################################
  if ((NTuple->genSignal != 0) || (NTuple->genMuMuBG != 0))
    for (unsigned int i = 0; i < NTuple->nB; i++)
    {
    // ###########################
    // # Check matching with mu- #
    // ###########################
    deltaEtaPhi = Utility->computeEtaPhiDistance (NTuple->genMumPx,NTuple->genMumPy,NTuple->genMumPz, NTuple->mumPx->at(i),NTuple->mumPy->at(i),NTuple->mumPz->at(i));
    NTuple->mumDeltaRwithMC->push_back(deltaEtaPhi);
    if (deltaEtaPhi < RCUTMU)
    {
      NTuple->truthMatchMum->push_back(1);
      if (printMsg) std::cout << __LINE__ << " : found matched mu-" << std::endl;
    }
    else NTuple->truthMatchMum->push_back(0);
              

    // ###########################
    // # Check matching with mu+ #
    // ###########################
    deltaEtaPhi = Utility->computeEtaPhiDistance (NTuple->genMupPx,NTuple->genMupPy,NTuple->genMupPz, NTuple->mupPx->at(i),NTuple->mupPy->at(i),NTuple->mupPz->at(i));
    NTuple->mupDeltaRwithMC->push_back(deltaEtaPhi);
    if (deltaEtaPhi < RCUTMU)
    {
      NTuple->truthMatchMup->push_back(1);
      if (printMsg) std::cout << __LINE__ << " : found matched mu+" << std::endl;
    }
    else NTuple->truthMatchMup->push_back(0);
              

    // ##################################
    // # Check matching with K*0 track- #
    // ##################################
    if ((NTuple->genSignal != 0) || (NTuple->genMuMuBGnTrksm != 0))
    {
      deltaEtaPhi = Utility->computeEtaPhiDistance (NTuple->genKstTrkmPx,NTuple->genKstTrkmPy,NTuple->genKstTrkmPz, NTuple->kstTrkmPx->at(i),NTuple->kstTrkmPy->at(i),NTuple->kstTrkmPz->at(i));
      NTuple->kstTrkmDeltaRwithMC->push_back(deltaEtaPhi);
      if (deltaEtaPhi < RCUTTRK)
      {
        NTuple->truthMatchTrkm->push_back(1);
        if (printMsg) std::cout << __LINE__ << " : found matched track-" << std::endl;
      }
      else NTuple->truthMatchTrkm->push_back(0);
    }
    else
    {
      NTuple->truthMatchTrkm->push_back(0);
      NTuple->kstTrkmDeltaRwithMC->push_back(-1.0);
    }


    // ##################################
    // # Check matching with K*0 track+ #
    // ##################################
    if ((NTuple->genSignal != 0) || (NTuple->genMuMuBGnTrksp != 0))
    {
      deltaEtaPhi = Utility->computeEtaPhiDistance (NTuple->genKstTrkpPx,NTuple->genKstTrkpPy,NTuple->genKstTrkpPz, NTuple->kstTrkpPx->at(i),NTuple->kstTrkpPy->at(i),NTuple->kstTrkpPz->at(i));
      NTuple->kstTrkpDeltaRwithMC->push_back(deltaEtaPhi);
      if (deltaEtaPhi < RCUTTRK) 
      {
        NTuple->truthMatchTrkp->push_back(1);
        if (printMsg) std::cout << __LINE__ << " : found matched track+" << std::endl;
      }
      else NTuple->truthMatchTrkp->push_back(0);
    }
    else
    {
      NTuple->truthMatchTrkp->push_back(0);
      NTuple->kstTrkpDeltaRwithMC->push_back(-1.0);
    }


    // ####################################################
    // # Check matching with B0 --> track+ track- mu+ mu- #
    // ####################################################
    if ((NTuple->truthMatchTrkm->back() == 1) && (NTuple->truthMatchTrkp->back() == 1) &&
        (NTuple->truthMatchMum->back() == 1) && (NTuple->truthMatchMup->back() == 1))
    {
      NTuple->truthMatchSignal->push_back(1);
      if (printMsg) std::cout << __LINE__ << " : @@@ Found matched B0 --> track+ track- mu+ mu- @@@" << std::endl;
    }
    else NTuple->truthMatchSignal->push_back(0);
  }
  else
    for (unsigned int i = 0; i < NTuple->nB; i++)
    {
      NTuple->mumDeltaRwithMC->push_back(-1.0);
      NTuple->mupDeltaRwithMC->push_back(-1.0);
      NTuple->kstTrkmDeltaRwithMC->push_back(-1.0);
      NTuple->kstTrkpDeltaRwithMC->push_back(-1.0);
    
      NTuple->truthMatchMum->push_back(0);
      NTuple->truthMatchMup->push_back(0);
      NTuple->truthMatchTrkm->push_back(0);
      NTuple->truthMatchTrkp->push_back(0);
      NTuple->truthMatchSignal->push_back(0);
    }


  // ###################
  // # Fill the ntuple #
  // ###################
  if (printMsg) std::cout << __LINE__ << " : @@@ Filling the tree @@@" << std::endl;
  theTree->Fill();
  NTuple->ClearNTuple();
}


std::string B0KstMuMu::getMuCat (reco::Muon const& muon)
{
  std::stringstream muCat;
  muCat.str("");

  if (muon.isGlobalMuon() == true)
  {
    muCat << " GlobalMuon";
    if (muon::isGoodMuon(muon, muon::GlobalMuonPromptTight) == true) muCat << " GlobalMuonPromptTight";
  }
  if (muon.isTrackerMuon() == true)
  {
    muCat << " TrackerMuon";
    if (muon::isGoodMuon(muon, muon::TrackerMuonArbitrated) == true) muCat << " TrackerMuonArbitrated";
    if (muon::isGoodMuon(muon, muon::TMLastStationTight)     == true) muCat << " TMLastStationTight";
    if (muon::isGoodMuon(muon, muon::TMLastStationLoose)     == true) muCat << " TMLastStationLoose";
    if (muon::isGoodMuon(muon, muon::TM2DCompatibilityTight) == true) muCat << " TM2DCompatibilityTight";
    if (muon::isGoodMuon(muon, muon::TM2DCompatibilityLoose) == true) muCat << " TM2DCompatibilityLoose";
    if (muon::isGoodMuon(muon, muon::TMOneStationTight)      == true) muCat << " TMOneStationTight";
    if (muon::isGoodMuon(muon, muon::TMOneStationLoose)      == true) muCat << " TMOneStationLoose";
    if (muon::isGoodMuon(muon, muon::TMLastStationAngTight)  == true) muCat << " TMLastStationAngTight";
    if (muon::isGoodMuon(muon, muon::TMLastStationAngLoose)  == true) muCat << " TMLastStationAngLoose";
    if (muon::isGoodMuon(muon, muon::TMOneStationAngTight)   == true) muCat << " TMOneStationAngTight";
    if (muon::isGoodMuon(muon, muon::TMOneStationAngLoose)   == true) muCat << " TMOneStationAngLoose";
  }
  if (muon.isStandAloneMuon() == true) muCat << " StandAloneMuon";
  if (muon.isCaloMuon()       == true) muCat << " CaloMuon";
  if ((muon.isGlobalMuon() == false) && (muon.isTrackerMuon() == false) && (muon.isStandAloneMuon() == false) && (muon.isCaloMuon() == false)) muCat << " NotInTable";

  return muCat.str();
}


const reco::Candidate* B0KstMuMu::skipOscillations (const reco::Candidate* Mother)
{
  if (abs(Mother->pdgId()) != 521)
    for (unsigned int i = 0; i < Mother->numberOfDaughters(); i++)
      if ((abs(Mother->daughter(i)->pdgId()) == 511) || (abs(Mother->daughter(i)->pdgId()) == 531) || (abs(Mother->daughter(i)->pdgId()) == 5122))
      {
        if (printMsg) std::cout << __LINE__ << " : @@@ Found oscillating B0/B0bar OR Bs/Bsbar OR Lambda_b/Lambda_bbar @@@" << std::endl;
        Mother = Mother->daughter(i);
      }

  return Mother;
}


const reco::Candidate* B0KstMuMu::findDaughter (edm::Handle<reco::GenParticleCollection> genParticles,
                        std::vector<std::vector<unsigned int>* >* Dau,
                        unsigned int it)
{
  const reco::Candidate* gen;
  if (Dau != NULL)
  {
    gen = skipOscillations(&(*genParticles)[(*Dau)[it]->operator[](0)]);
    for (unsigned int i = 1; i < (*Dau)[it]->size(); i++) gen = gen->daughter((*Dau)[it]->operator[](i));
    return gen;
  }
  else return NULL;
}


void B0KstMuMu::searchForStableChargedDaughters (const reco::Candidate* Mother,
                         unsigned int motherIndex,
                         std::vector<std::vector<unsigned int>* >* posDau,
                         std::vector<std::vector<unsigned int>* >* negDau)
{
  const reco::Candidate* genDau;
  unsigned int sizeNegVec = negDau->size();
  unsigned int sizePosVec = posDau->size();

  for (unsigned int i = 0; i < Mother->numberOfDaughters(); i++)
  {
    genDau = Mother->daughter(i);

    if (printMsg) std::cout << __LINE__ << " : start exploring daughter: " << genDau->pdgId() << "\tindex: " << i << "\tof mother: " << Mother->pdgId() << std::endl;

    if ((genDau->pdgId() == 11)   ||
        (genDau->pdgId() == 13)   ||
        (genDau->pdgId() == -211) ||
        (genDau->pdgId() == -321) ||
        (genDau->pdgId() == -2212))
    {
      std::vector<unsigned int>* vecDau = new std::vector<unsigned int>;
      vecDau->push_back(i);
      negDau->push_back(vecDau);
      if (printMsg) std::cout << __LINE__ << " : found possible background negative track: " << genDau->pdgId() << "\tfrom mother: " << genDau->mother()->pdgId() << std::endl;
    }
    else if ((genDau->pdgId() == -11) ||
             (genDau->pdgId() == -13) ||
             (genDau->pdgId() == 211) ||
             (genDau->pdgId() == 321) ||
             (genDau->pdgId() == 2212))
    {
      std::vector<unsigned int>* vecDau = new std::vector<unsigned int>;
      vecDau->push_back(i);
      posDau->push_back(vecDau);
      if (printMsg) std::cout << __LINE__ << " : found possible background positive track: " << genDau->pdgId() << "\tfrom mother: " << genDau->mother()->pdgId() << std::endl;
    }
    else if ((abs(genDau->pdgId()) != 311) &&
             (abs(genDau->pdgId()) != 310) &&
             (abs(genDau->pdgId()) != 130) &&
             (abs(genDau->pdgId()) != 22) &&
             (abs(genDau->pdgId()) != 12) &&
             (abs(genDau->pdgId()) != 14) &&
             (abs(genDau->pdgId()) != 16))
      searchForStableChargedDaughters (genDau, i, posDau, negDau);
    else if (printMsg) std::cout << __LINE__ << " : found long living neutral particle: " << genDau->pdgId() << std::endl;
    }

  for (unsigned int it = negDau->size(); it > sizeNegVec; it--) negDau->operator[](it-1)->insert(negDau->operator[](it-1)->begin(), motherIndex);
  for (unsigned int it = posDau->size(); it > sizePosVec; it--) posDau->operator[](it-1)->insert(posDau->operator[](it-1)->begin(), motherIndex);
}


void B0KstMuMu::beginJob ()
{
  edm::Service<TFileService> fs;
  theTree = fs->make<TTree>("B0KstMuMuNTuple","B0KstMuMuNTuple");
  NTuple->MakeTreeBranches(theTree);

  // ############################
  // # Print out HLT-trigger cuts #
  // ############################
  std::cout << "\n@@@ Pre-selection cuts @@@" << std::endl;
  std::cout << __LINE__ << " : CLMUMUVTX      = " << CLMUMUVTX      << std::endl;
  std::cout << __LINE__ << " : LSMUMUBS       = " << LSMUMUBS       << std::endl;
  std::cout << __LINE__ << " : DCAMUMU        = " << DCAMUMU        << std::endl;
  std::cout << __LINE__ << " : DCAMUBS        = " << DCAMUBS        << std::endl;
  std::cout << __LINE__ << " : COSALPHAMUMUBS = " << COSALPHAMUMUBS << std::endl;
  std::cout << __LINE__ << " : MUMINPT        = " << MUMINPT        << std::endl;
  std::cout << __LINE__ << " : MUMAXETA       = " << MUMAXETA       << std::endl;
  std::cout << __LINE__ << " : MINMUMUPT      = " << MINMUMUPT      << std::endl;
  std::cout << __LINE__ << " : MINMUMUINVMASS = " << MINMUMUINVMASS << std::endl;
  std::cout << __LINE__ << " : MAXMUMUINVMASS = " << MAXMUMUINVMASS << std::endl;
  
  // ##############################
  // # Print out pre-selection cuts #
  // ##############################
  std::cout << __LINE__ << " : B0MASSLOWLIMIT = " << B0MASSLOWLIMIT << std::endl;
  std::cout << __LINE__ << " : B0MASSUPLIMIT  = " << B0MASSUPLIMIT  << std::endl;
  std::cout << __LINE__ << " : CLB0VTX        = " << CLB0VTX        << std::endl;
  std::cout << __LINE__ << " : KSTMASSWINDOW  = " << KSTMASSWINDOW  << std::endl;
  std::cout << __LINE__ << " : HADDCASBS      = " << HADDCASBS      << std::endl;
  std::cout << __LINE__ << " : MINHADPT       = " << MINHADPT       << std::endl;


  std::cout << "\n@@@ Global constants @@@" << std::endl;
  std::cout << __LINE__ << " : TRKMAXR    = " << TRKMAXR << std::endl;
  
  std::cout << __LINE__ << " : PRIVTXNDOF = " << PRIVTXNDOF << std::endl;
  std::cout << __LINE__ << " : PRIVTXMAXZ = " << PRIVTXMAXZ << std::endl;
  std::cout << __LINE__ << " : PRIVTXMAXR = " << PRIVTXMAXR << std::endl;
  
  std::cout << __LINE__ << " : MUVARTOLE  = " << MUVARTOLE << std::endl;
  std::cout << __LINE__ << " : HADVARTOLE = " << HADVARTOLE << std::endl;

  std::cout << __LINE__ << " : RCUTMU     = " << RCUTMU << std::endl;
  std::cout << __LINE__ << " : RCUTTRK    = " << RCUTTRK << std::endl;
}


void B0KstMuMu::endJob ()
{
  theTree->GetDirectory()->cd();
  theTree->Write();
}


void B0KstMuMu::endLuminosityBlock(const edm::LuminosityBlock& lumiBlock, const edm::EventSetup& iSetup)
{
  if ((doGenReco_ == 2) || (doGenReco_ == 3))
    {
      edm::Handle<GenFilterInfo> genFilter;
      lumiBlock.getByToken(genFilterToken_, genFilter);

      NTuple->numEventsTried  = genFilter->numEventsTried();
      NTuple->numEventsPassed = genFilter->numEventsPassed();

      if (printMsg)
    {
      std::cout << "\n@@@ End of a luminosity block @@@" << std::endl;
      std::cout << __LINE__ << " : number of events tried  = " << NTuple->numEventsTried << std::endl;
      std::cout << __LINE__ << " : number of events passed = " << NTuple->numEventsPassed << std::endl;
    }
    }
}

float B0KstMuMu::calculateIsolation( ClosestApproachInRPhi ClosestApp, 
                                     reco::TransientTrack muTrackmTT, 
                                     reco::TransientTrack TrackIsoTT,
                                     const ParticleMass muonMass,
                                     float vtx_x,
                                     float vtx_y,
                                     float vtx_z
                                    )
{
  float iso = 0;

  KinematicParticleVertexFitter KPVtxFitter; 
  KinematicParticleFactoryFromTransientTrack partFactory;

  float chi, ndf, alpha;
  TLorentzVector tmp_trk_lv, tmp_cand_lv, tmp_lv;
  ClosestApp.calculate(TrackIsoTT.initialFreeState(), muTrackmTT.initialFreeState());

  if ( ClosestApp.status() && ClosestApp.distance() < 0.1 )
  {
  
    chi = 0.;
    ndf = 0.;
    float mumasserr = 3.5e-9;
    std::vector<RefCountedKinematicParticle> tmpParticles;
    tmpParticles.push_back(partFactory.particle(muTrackmTT, muonMass,chi,ndf, mumasserr) );
    tmpParticles.push_back(partFactory.particle(TrackIsoTT, muonMass,chi,ndf, mumasserr) );
    
    RefCountedKinematicTree tmpVertexFitTree = KPVtxFitter.fit(tmpParticles);
    if ( ! tmpVertexFitTree->isValid()) return 0;
    tmpVertexFitTree->movePointerToTheTop();
    RefCountedKinematicVertex tmpVertex   = tmpVertexFitTree->currentDecayVertex();
    if (tmpVertexFitTree->isValid() && tmpVertex -> vertexIsValid() )
    {
      tmp_trk_lv .SetPtEtaPhiM( TrackIsoTT.track().pt(), TrackIsoTT.track().eta(), TrackIsoTT.track().phi(), Utility -> kaonMass);
      tmp_cand_lv.SetPtEtaPhiM( muTrackmTT.track().pt(), muTrackmTT.track().eta(), muTrackmTT.track().phi(), Utility -> kaonMass);
      tmp_lv =  tmp_trk_lv + tmp_cand_lv;
      
      alpha = tmp_lv.Angle(TVector3 ( tmpVertex->position().x() - vtx_x, 
                                      tmpVertex->position().y() - vtx_y, 
                                      tmpVertex->position().z() - vtx_z )
                                    );
      iso = (tmp_lv).P() * alpha / ( (tmp_lv).P() + tmp_cand_lv.Pt() + tmp_trk_lv.Pt() );
    }
  }

  return iso;
}

// Compute impact parameter 3D wrt Reco Vtx
std::pair<double,double> B0KstMuMu::pionImpactParameter(reco::TransientTrack piTT, reco::Vertex myVtx)
{
    std::pair<double,double> measure;
    std::pair<bool,Measurement1D>  piIP_pair = IPTools::absoluteImpactParameter3D(piTT, myVtx);
    if (piIP_pair.first)
    {
      measure.first  = piIP_pair.second.value();
      measure.second = piIP_pair.second.significance();
    }
    else 
    {
      measure.first  = 0;
      measure.second = 0;
    } 
    return measure;
}


// Define this as a plug-in
DEFINE_FWK_MODULE(B0KstMuMu);
