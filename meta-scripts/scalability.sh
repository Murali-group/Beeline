
# script to evaluate the time and memory (using /usr/bin/time --verbose) usage of methods when subsetting 100 to 10,000 genes

#config_file="config-files/2019-07-datasets/Camp2017.yaml"
#declare -a datasets=("Camp2017" "Chu2016" "Nestorowa2016")
declare -a datasets=("Chu2016" "Nestorowa2016")
#declare -a num_genes_list=(100 500 1000 2000 5000)
#declare -a num_genes_list=(2000 5000)
declare -a num_genes_list=(100 500 1000)
#algs="--alg SINCERITIES --alg GRNBOOST2 --alg PPCOR --alg PIDC"
#algs="--alg GRNBOOST2 --alg PPCOR --alg PIDC"
algs="--alg GENIE3 --alg LEAP --alg SCINGE"
#algs="--alg GENIE3 --alg LEAP"
#algs="--alg SCODE --alg SCRIBE --alg GRISLI"
#algs="--alg SCRIBE"

for dataset in ${datasets[*]}; do
    config_file="config-files/2019-07-datasets/$dataset.yaml"
    for num_genes in ${num_genes_list[*]}; do
        echo ""
        echo "----------------------------------------------------------------------------------------------------" 
        date
        echo "python -u subset_genes.py --config $config_file --most-variable-genes --num-genes $num_genes $algs --run-algs "
        python -u subset_genes.py --config $config_file --most-variable-genes --num-genes $num_genes $algs --run-algs 
    done
done

echo "DONE DONE DONE"
date

# TODO send an email when done?

