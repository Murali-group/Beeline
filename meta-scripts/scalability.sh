
# script to evaluate the time and memory (using /usr/bin/time --verbose) usage of methods when subsetting 100 to 10,000 genes
#python=/data/jeff-law/tools/anaconda3/bin/python
python=python
#python="/home/adyprat/anaconda3/bin/python"

# human datasets
declare -a datasets=("Camp2017" "Chu2016")
tfs_file="/data/inputs/single-cell/datasets/TFs/human/human-tfs.csv"
# mouse datasets
#declare -a datasets=("Nestorowa2016-E" "Nestorowa2016-GM" "Nestorowa2016-L" "Shalek2014" "Hayashi2018")
#tfs_file="/data/inputs/single-cell/datasets/TFs/mouse/mouse-tfs.csv"

#declare -a num_genes_list=(100 500 1000 2000 5000)
#declare -a num_genes_list=(1000)
#declare -a num_genes_list=(500 1000)
declare -a num_genes_list=(500 1000)
#algs="--alg SINCERITIES --alg GRNBOOST2 --alg PPCOR --alg SCRIBE --alg GENIE3 --alg PIDC --alg LEAP"

# this will run each of the methods sequentially to ensure a fiar comparison of running time and memory usage
for dataset in ${datasets[*]}; do
    #config_file="config-files/2019-08-datasets/$dataset.yaml"
    config_file="config-files/2019-08-datasets-updated/$dataset.yaml"
    for num_genes in ${num_genes_list[*]}; do
        echo ""
        echo "----------------------------------------------------------------------------------------------------" 
        date
        #echo "python -u subset_genes.py --config $config_file --most-variable-genes --num-genes $num_genes $algs --run-algs "
        #python -u subset_genes.py --config $config_file --most-variable-genes --num-genes $num_genes $algs --run-algs 
        #python -u subset_genes.py --config $config_file --most-variable-genes --num-genes $num_genes $algs 
        cmd="""$python subset_genes.py \
            --config $config_file \
            --most-variable-genes \
            --num-genes $num_genes \
            --pval-cutoff 0.01 \
            --bf-corr \
            --sort-by-variance \
            --include-tfs $tfs_file \

        """
            #--forced
            # don't run the algs here. Use run_parallel.py to run them in parallel by starting a screen session for each one
            #--run-algs $algs \
        echo "$cmd"
        $cmd
    done
done

echo "DONE DONE DONE"
date

# TODO send an email when done?

