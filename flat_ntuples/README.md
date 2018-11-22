* **flatNtuples.py** takes as input the ntuples produced by the crab and prepares a flattened version (one entry per each reco B0 candidate)
* **addNvtxWeight.py** adds the variable "weight" to the MC, based on the data distribution. Usage is   
```
python addNVtxWeight.py "fileData" "fileMC" outFileName -t ntuple>
```
* **flatNtuplesGENMC.py** only reads and write the GEN quantities. Also adds angular variables
