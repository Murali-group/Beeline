## Introduce the hyper-parameter used in DeepSEM

- n_epochs: Number of Epochs for training DeepSEM.
- task: Determine which task to run. Select from (non_celltype_GRN,celltype_GRN,embedding,simulation)
- setting: Determine whether or not to use the default hyper-parameter
- batch_size: The batch size used in the training process.
- data_file: The input scRNA-seq gene expression file.
- net_file: The ground truth of GRN. Only used in GRN inference task if available. 
- alpha: The loss coefficient for L1 norm of W, which is same as \alpha used in our paper.
- beta: The loss coefficient for KL term (beta-VAE), which is same as \beta used in our paper.
- lr: The learning rate of used for RMSprop.
- lr_step_size: The step size of LR decay.
- gamma: The decay factor of learning rate.
- n_hidden: The Number of hidden neural used in MLP
- K: The number of Gaussian kernel in GMM, default =1
- K1: The Number of epoch for optimize MLP. Notes that we optimize MLP and W alternately. The default setting denotes
 to optimize MLP for one epoch then optimize W for two epochs.
- K2: The Number of epoch for optimize W. Notes that we optimize MLP and W alternately. The default setting denotes
 to optimize MLP for one epoch then optimize W for two epochs.
- save_name: The directory to save the output file of DeepSEM


## Default hyper-parameter selection

We provide default hyper-parameter in main.py. Therefore, for most situations, directly run main.py with default 
hyper-parameter can perform quite well.

Note that, the hyper-parameters used in Figure 2 and Supplementary Figure 1 is **not default hyper-parameters**, which is 
selected by cross validation between cell line for fair comparison. For exactly reproducing performance shown in 
Figure 2 and Supplementary Figure 1. Please used the following hyper-parameters:

- **Celltype non-specific GRN inference alpha=100,beta=1,n_epochs=90**
- **Celltype specific GRN inference alpha=0.1 beta=0.01 n_epochs=150 for hESC and mHSC-E and alpha=1,beta=0.01,
n_epochs=150 for others**.

You can also use the default hyper-parameter, the overall average performance between the default hyper-parameters and 
hyper-parameter determined by cross-validation should be less than 0.5%

## Important hyper-parameter and suggestion for hyper-parameter tuning

Alpha and beta are most important hyper-parameters among all hyper-parameters. Hyper-parameter beta balance the 
novelty and similarity when comparision to real cells. More specific hyper-parameter sweeping results are shown in our manuscript (Supplementary Section 7). 
 
Suggestion for hyper-parameter tuning: If the MSE loss does not decrease significantly for a long time during the 
training process. Slightly decrease the hyper-parameter beta can solve this problem (hardly ever happened in our 
experiments). 
Suggestion for large dataset: Note that the cell number of all benchmark datasets is less than 1500, so we suggest to subsample cells for each epoch and rational increase number of epoch when for large dataset. We STRONGLY recommend to select TF+highly variable gene for larget dataset.
