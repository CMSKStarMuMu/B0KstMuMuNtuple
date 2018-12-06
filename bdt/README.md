**create_subsample.C**  
Create N exclusive subsamples to be used for the BDT training/evaluation

**final_bdt_sub_samples.py**  
Train the BDT over N-1 subsamples created with the above macro, and save the results as classifiers/tag.pkl  
Performance plots are saved into classifiers/plots/

**add_bdt_subsamples.py**  
Add to each origina subsample a branch  for the BDT-score (named "bdt_prob") and a boolean for the pre-selection ("pass_preselection").  
The BDT score for sample i is evaluated on the classifier trained on all the other subsamples (excluding subsample i).   

