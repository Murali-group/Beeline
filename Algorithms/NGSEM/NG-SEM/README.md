# NG-SEM

## Description
This repository is an implementation for **NG-SEM** (short for, an effective **N**on-**G**aussian **S**tructural **E**quation **M**odeling framework for gene regulatory network inference from single cell RNA-seq data).


## Installation

`NG-SEM` can be installed from github directly as follows:

  ```R
  if (!require("devtools")){
        install.packages("devtools")
  }
  install_github("jiayingzhao/NG-SEM")
  ```
  
## Quick Start
`ng_sem()` is the wrapper function to perform the GRN inference.

To perform this method, simply run:

  ```R
  library(NGSEM)
  
  res <- ng_sem(dat, nk=3, param.miter=500, param.error=1e-6, n.cores = 24, seed_list = NULL)

  ```
**ng_sem() parameters**
  
- `dat`: a gene-by-cell matrix (or data.frame, **required**) for single-cell RNA-seq expression data, which is not restricted to experiment measurements and is applicable to RPM (reads per million reads), TPM (transcripts per kilobase per millions reads) or RPKM (reads per kilobase per millions reads).
- `nk`: the number of Gaussian components for the Gaussian-mixture noise model (`default: 3`).
- `param.miter`: the maximum iteration to implement the EM algorithm (`default: 500`).
- `param.error`: the error bound for the termination criteria to implement the EM algorithm (`default: 1e-6`).
- `n.cores`: the number of working CPU cores for parallelization (`default: 24`).
- `seed_list`: seeds for reproducibility (`default: NULL`). The seed list for reproducing can be found as an *rds* file under the *./data* folder of this repository 



`ng_sem()` returns the following list of values:

- `w.mat`: a G-by-G GRN adjacency matrix.
- `likelihood`: a list records the likelihood.
- `iteration`: a vector of length G that records the iterations.
- `Comp_Weight`: a G-by-nk matrix that contains the Gaussian component weights.
- `Comp_Sigma2`: a G-by-nk matrix that contains the Gaussian component variances.
- `time`: CPU time to implement **NG-SEM**.


## Contact

**Jiaying Zhao** (HKU) jyzhao@connect.hku.hk
