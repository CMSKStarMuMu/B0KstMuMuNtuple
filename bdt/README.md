create_subsample.C
create N exclusive subsamples to be used for the BDT training/evaluation

final_bdt_sub_samples.py
train the BDT over N-1 subsamples and save the results as classifiers/tag.pkl
performance plots are saved into classifiers/plots/

add_bdt_subsamples.py
add the BDT-score branch (named "bdt_prob") and a boolean for the pre-selection ("pass_preselection") to the original subsamples.
The BDT score for sample i is evaluated on the classifier trained on all the other subsamples (excluding subsample i). 

