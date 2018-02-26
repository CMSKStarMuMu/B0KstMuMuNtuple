#ifndef UTILS_H
#define UTILS_H

#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TF2.h>
#include <TF3.h>
#include <TMatrixTSym.h>
#include <TGraphAsymmErrors.h>

#if ROOFIT
#include <RooRealVar.h>
#include <RooFitResult.h>
#endif

#include <string>
#include <vector>

#include "B0KstMuMuTreeContent.h"


class Utils
{
 public:
  
  Utils(bool rightFlavorTag = true);
  ~Utils() {};

  double computeInvMass (double Px1,
			 double Py1,
			 double Pz1,
			 double mass1,
			 double Px2,
			 double Py2,
			 double Pz2,
			 double mass2,
			 double Px3 = 0,
			 double Py3 = 0,
			 double Pz3 = 0,
			 double mass3 = 0);

  double computeEta (double Px,
                     double Py,
                     double Pz);
  
  double computePhi (double Px, double Py);
  
  double computeEtaPhiDistance (double Px1,
				double Py1,
				double Pz1,
                                double Px2,
                                double Py2,
				double Pz2);
  
  void computeLS (double Vx,
		  double Vy,
		  double Vz,
		  double Wx,
		  double Wy,
		  double Wz,
		  double VxErr2,
		  double VyErr2,
		  double VzErr2,
		  double VxyCov,
		  double VxzCov,
		  double VyzCov,
		  double WxErr2,
		  double WyErr2,
		  double WzErr2,
		  double WxyCov,
		  double WxzCov,
		  double WyzCov,
		  double* deltaD,
		  double* deltaDErr);

  void computeCosAlpha (double Vx,
			double Vy,
			double Vz,
			double Wx,
			double Wy,
			double Wz,
			double VxErr2,
			double VyErr2,
			double VzErr2,
			double VxyCov,
			double VxzCov,
			double VyzCov,
			double WxErr2,
			double WyErr2,
			double WzErr2,
			double WxyCov,
			double WxzCov,
			double WyzCov,
			double* cosAlpha,
			double* cosAlphaErr);

  unsigned int IsInTriggerTable     (B0KstMuMuTreeContent* NTupleIn, double* HLTCutVar1, double* HLTCutVar2, int index = 0, double evtFrac = -1.0);
  double muonMass;
  double pionMass;
  double kaonMass;
  double protonMass;
  double kstMass;
  double B0Mass;
  double JPsiMass;
  double PsiPMass;
  double D0Mass;
  
  double JPsiBF;
  double JPsiKpiBF;
  double KstMuMuBF;
  double KstKpiMuMuBF;
  double PsiPBF;
  double PsiPKpiBF;

  double muonMassErr;
  double pionMassErr;
  double kaonMassErr;
  double B0MassErr;
  double kstSigma;

  double PI;

  bool RIGHTflavorTAG;

  int B0ToKstMuMu;
  int B0ToJPsiKst;
  int B0ToPsi2SKst;


 private:

  TF1* KstMassShape;

  std::vector<std::string> TrigTable;

  double ProbThreshold;
  double scrambleFraction;


};

#endif
