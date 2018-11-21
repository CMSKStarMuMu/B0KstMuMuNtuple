Ntuple production for 2016 data  


cmsrel CMSSW_8_0_24  
cd CMSSW_8_0_24/src  
git cms-addpkg PhysicsTools/PatAlgos/  
git remote add patTrk git@github.com:sarafiorendi/cmssw.git  
git fetch patTrk  
git checkout patTrk/pat4tracks_Kmumu PhysicsTools/PatAlgos/plugins/PATTriggerMatchEmbedder.cc  
git checkout patTrk/pat4tracks_Kmumu PhysicsTools/PatAlgos/python/tools/trigTools.py  
scram b -j 6  

mkdir B0KstarMM  
cd B0KstarMM  
git clone git@github.com:CMSKStarMuMu/B0KstMuMu.git .  
cd B0KstMuMu  
scram b -j 6  


