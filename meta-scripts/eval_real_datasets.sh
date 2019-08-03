
# script to evaluate the methods using a variety of different networks
# if running on cuthbert
python=/data/jeff-law/tools/anaconda3/bin/python
# if running on the csb machines
#python=python

base_net_path="/data/inputs/single-cell/datasets"
# whichever algorithms have finished can be set here
algs="--alg SINCERITIES --alg GRNBOOST2 --alg PPCOR --alg PIDC --alg GENIE3 --alg LEAP --alg SCINGE --alg SCODE --alg SCRIBE --alg GRISLI"

#declare -a num_genes_list=(100 500 1000 2000 5000)
declare -a num_genes_list=(100 500)

# human datasets with gene names
declare -a datasets=("Camp2017" "Chu2016")
#declare -a datasets=("Chu2016")
# non-cell-type-specific networks
declare -a networks=(
   # "${base_net_path}/string/human/9606.protein.links.v11.0-gene-names-c400.csv"
   # "${base_net_path}/string/human/9606.protein.links.v11.0-gene-names-c700.csv"
    #"${base_net_path}/TRRUST/processed/trrust_data.human.csv"
    #"${base_net_path}/DoRothEA/ref-network-A.csv"
    #"${base_net_path}/DoRothEA/ref-network-A-B.csv"
    #"${base_net_path}/DoRothEA/ref-network-A-B-C.csv"
    "${base_net_path}/RegNetwork/human-net.csv"
    )

for dataset in ${datasets[*]}; do
    for num_genes in ${num_genes_list[*]}; do
    for net in ${networks[*]}; do
        config_file="config-files/2019-07-datasets/$dataset/${dataset}-${num_genes}.yaml"
        echo "----------------------------------------------------------------------------------------------------" 
  echo "$python -u eval_net.py --config $config_file --ref-net-file $net $algs --auc  --time --epr"
        #$python -u eval_net.py --config $config_file --ref-net-file $net $algs --auc  --time --epr
        $python -u eval_net.py --config $config_file --ref-net-file $net $algs --auc  --time --epr --stats-only
    done
    done
done
echo "skipping mouse"
exit

# ------------------------------------------------------
# ------------------------------------------------------
# now the mouse datasets with Ensembl IDs 
declare -a datasets=("Nestorowa2016")
declare -a networks=(
    #"${base_net_path}/string/mouse/10090.protein.links.v11.0-ensembl-ids-c400.csv"
    #"${base_net_path}/string/mouse/10090.protein.links.v11.0-ensembl-ids-c700.csv"
    "${base_net_path}/TRRUST/processed/trrust_data.mouse.csv"
    )

for dataset in ${datasets[*]}; do
    for num_genes in ${num_genes_list[*]}; do
    for net in ${networks[*]}; do
        config_file="config-files/2019-07-datasets/$dataset/${dataset}-${num_genes}.yaml"
        echo "----------------------------------------------------------------------------------------------------" 
  echo "$python -u eval_net.py --config $config_file --ref-net-file $net $algs --auc  --time --epr"
        $python -u eval_net.py --config $config_file --ref-net-file $net $algs --auc  --time --epr
    done
    done
done
