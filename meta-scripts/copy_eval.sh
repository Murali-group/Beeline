
# script to evaluate the time and memory (using /usr/bin/time --verbose) usage of methods when subsetting 100 to 10,000 genes
#python=/data/jeff-law/tools/anaconda3/bin/python
python=python
#python="/home/adyprat/anaconda3/bin/python"

# human datasets
#declare -a datasets=("Camp2017" "Chu2016")
#tfs_file="/data/inputs/single-cell/datasets/TFs/human/human-tfs.csv"
# mouse datasets
declare -a datasets=("Camp2017" "Chu2016" "Nestorowa2016-E" "Nestorowa2016-GM" "Nestorowa2016-L" "Shalek2014" "Hayashi2018")
#tfs_file="/data/inputs/single-cell/datasets/TFs/mouse/mouse-tfs.csv"

#declare -a num_genes_list=(100 500 1000 2000 5000)
#declare -a num_genes_list=(100 500)
declare -a num_genes_list=(
"500_varsort" 
"500_varsort_tfs" 
"1000_varsort" 
"1000_varsort_tfs" 
)
#algs="--alg SINCERITIES --alg GRNBOOST2 --alg PPCOR --alg SCRIBE --alg GENIE3 --alg PIDC --alg LEAP"
#algs="--alg SINCERITIES --alg LEAP"
#algs="--alg SINCERITIES"
#algs=""
best_params_file="outputs/2019-08-datasets-updated/parameter-search-2019-08-real-data-updated.csv"

# this will run each of the methods sequentially to ensure a fiar comparison of running time and memory usage
for dataset in ${datasets[*]}; do
    for num_genes in ${num_genes_list[*]}; do
        config_file="config-files/2019-08-datasets-updated/${dataset}/${dataset}-0.01_BF_${num_genes}.yaml"
        #echo ""
        echo "----------------------------------------------------------------------------------------------------" 
        #date
        #echo "python -u subset_genes.py --config $config_file --most-variable-genes --num-genes $num_genes $algs --run-algs "
        #python -u subset_genes.py --config $config_file --most-variable-genes --num-genes $num_genes $algs --run-algs 
        #python -u subset_genes.py --config $config_file --most-variable-genes --num-genes $num_genes $algs 
        #cmd="""$python masterscript.py \
        #    --config $config_file \
        #    $algs \
        #    --eval-only \
#
        #"""
            #--param-per-job
        new_params_file="${best_params_file}-$num_genes.csv"
        head -n 1 $best_params_file > $new_params_file
        #grep "$num_genes\$" $best_params_file | grep "SINCERITIES" >> $new_params_file
        grep "$num_genes\$" $best_params_file >> $new_params_file
        echo """grep \"$num_genes\$\" $best_params_file >> $new_params_file"""
        # limit the best params file to teh current num genes
        cmd="""python setup_best_params.py \
            --config $config_file \
            --best-params-file $new_params_file \

        """
            #--copy-out-files /data/jeff-law/projects/2019-04-single-cell/RNMethods/outputs/2019-08-datasets-best-params/
        echo "$cmd"
        $cmd
        #exit
    done
done

echo "DONE DONE DONE"
date

# TODO send an email when done?



