
# script to evaluate the time and memory (using /usr/bin/time --verbose) usage of methods when subsetting 100 to 10,000 genes
#python=/data/jeff-law/tools/anaconda3/bin/python
python=python

declare -a datasets=("Camp2017" "Chu2016" "Nestorowa2016-E" "Nestorowa2016-GM" "Nestorowa2016-L" "Shalek2014" "Hayashi2018")

declare -a num_genes_list=(
#"500_varsort" 
#"500_varsort_tfs" 
"1000_varsort" 
"1000_varsort_tfs" 
)
#algs="--alg SINCERITIES --alg GRNBOOST2 --alg PPCOR --alg SCRIBE --alg GENIE3 --alg PIDC --alg LEAP"
#algs="--alg SINCERITIES --alg LEAP"
algs="--alg SINCERITIES"
#algs="--alg SCRIBE"
#algs="--alg SCODE"

# this will run each of the methods sequentially to ensure a fiar comparison of running time and memory usage
for dataset in ${datasets[*]}; do
    for num_genes in ${num_genes_list[*]}; do
        config_file="config-files/2019-08-datasets-updated/${dataset}/${dataset}-0.01_BF_${num_genes}.yaml"
        echo ""
        echo "----------------------------------------------------------------------------------------------------" 
        date
        #echo "python -u subset_genes.py --config $config_file --most-variable-genes --num-genes $num_genes $algs --run-algs "
        #python -u subset_genes.py --config $config_file --most-variable-genes --num-genes $num_genes $algs --run-algs 
        #python -u subset_genes.py --config $config_file --most-variable-genes --num-genes $num_genes $algs 
        cmd="""$python masterscript.py \
            --config $config_file \
            $algs \
            --param-per-job

        """
            #--eval-only \
        echo "$cmd"
        $cmd
    done
done

echo "DONE DONE DONE"
date

# TODO send an email when done?


