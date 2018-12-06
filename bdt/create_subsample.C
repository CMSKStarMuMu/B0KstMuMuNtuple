#include "TROOT.h"
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"

#include <iostream>
#include <algorithm>
#include <vector>
#include <math.h> 

int n_subsamples = 11;

void create_subsample(){
  
  TChain* tree = new TChain("ntuple");
  tree -> Add("../../test/flat_ntuples/ntuples/ntuple_flat_LMNR_PostRefitMomenta_2016BCDEFGH_skimSoftMu.root");
  
  std::vector<TFile*> out_files;
  for (int i=0; i < n_subsamples; i++){
      TFile* ifile = new TFile(Form("sub_samples/sample_data_LMNR_events_%d.root", i),"RECREATE");
      out_files.push_back(ifile);
      std::cout << "output file : " << out_files.at(i) -> GetName() << std::endl;
  }
  
  tree -> SetBranchStatus("*",1);
  double eventN;
  tree -> SetBranchAddress( "eventN", &eventN);
  int nentries = tree->GetEntries();
  std::cout << "original n entries: " << nentries << std::endl;

  std::vector<TTree*> out_trees;
  std::vector<bool> keepentries;
  for (int i = 0; i < n_subsamples; i++){
      out_files.at(i) -> cd();
      TTree* itree = (TTree*) tree->CloneTree(0);
      out_trees.push_back(itree);
      out_trees.at(i)->SetAutoSave(1000000);
      keepentries.push_back(false);
  }

  double q   = 1.;
  double mod, fractpart, intpart;
  for (Int_t eventNo=0; eventNo < nentries; eventNo++)
  {
    Int_t IgetEvent   = tree   -> GetEntry(eventNo);
    
    q         = fabs(eventN) / n_subsamples;
    fractpart = modf (q, &intpart);
    mod       = fabs(eventN) - n_subsamples*intpart;

    out_trees.at(mod) -> Fill();
  }

  for (int i = 0; i < n_subsamples; i++){
      out_files.at(i) -> cd();
      out_trees.at(i) -> Write();
  }  

  return;
}