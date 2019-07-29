# Parameter Search

## Run the methods
This branch contains the code for performing the parameter search. 
The main difference between this branch and the main branch is that the parameters used when running a given method are appended to the output file (e.g., "nbins8" for SINCERITES).

To re-run the parameter estimation for the dyn-LI model, use this command:
```
python masterscript.py --config config-files/2019-07-dyn-LI/dyn-LI-500.yaml --alg LEAP --alg SCODE --alg SINCERITIES  --alg SCRIBE
```
The script will start a screen session for each method, and run them on the 10 datasets present in yaml file, using the parameters in the yaml file.  

Since GRISLI and SCINGE have many more parameters, I setup `masterscript.py` as well as `grisliRunner.py` and `scingeRunner.py` to submit them to our cluster named baobab using `qsub`, and run them without docker. 
```
ssh baobab
python masterscript.py --config config-files/2019-07-dyn-LI/dyn-LI-500.yaml --alg GRISLI --alg SCINGE --qsub
```

## Visualize the AUROC for each method
To visualize the results and make a table of the best parameters for each method, use the jupyter notebook located here:
```
src/jupyter-notebooks/benchmark-methods.ipynb
```
