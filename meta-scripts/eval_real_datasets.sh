
# script to evaluate the methods using a variety of different networks
# if running on cuthbert
#python=/data/jeff-law/tools/anaconda3/bin/python
# if running on the csb machines
python=python
# use adyprat's python for now
#python=/home/adyprat/anaconda3/bin/python

base_net_path="/data/inputs/single-cell/datasets"
# whichever algorithms have finished can be set here
#algs="--alg SINCERITIES --alg GRNBOOST2 --alg PPCOR --alg PIDC --alg GENIE3 --alg LEAP --alg SCINGE --alg SCODE --alg SCRIBE --alg GRISLI"
#algs="--alg SINCERITIES --alg LEAP --alg SCODE --alg SCRIBE"
algs="--alg SINCERITIES --alg LEAP --alg SCRIBE"

#declare -a num_genes_list=(100 500 1000 2000 5000)
#declare -a num_genes_list=("0.01_BF_100_varsort_tfs")
declare -a num_genes_list=(
"0.01_BF_500_varsort"
#"0.01_BF_500_varsort_tfs"
"0.01_BF_1000_varsort"
#"0.01_BF_1000_varsort_tfs"
)

# human datasets with gene names
declare -a datasets=("Camp2017" "Chu2016")
# non-cell-type-specific networks
declare -a networks=(
   # "${base_net_path}/string/human/9606.protein.links.v11.0-gene-names-c400.csv"
   # "${base_net_path}/string/human/9606.protein.links.v11.0-gene-names-c700.csv"
    #"${base_net_path}/TRRUST/processed/trrust_data.human.csv"
    #"${base_net_path}/DoRothEA/ref-network-A.csv"
    #"${base_net_path}/DoRothEA/ref-network-A-B.csv"
    #"${base_net_path}/DoRothEA/ref-network-A-B-C.csv"
    #"${base_net_path}/RegNetwork/human-net.csv"
    "${base_net_path}/Final/human/human-net.csv"
    )

for dataset in ${datasets[*]}; do
    for num_genes in ${num_genes_list[*]}; do
    for net in ${networks[*]}; do
        config_file="config-files/2019-08-datasets-updated/$dataset/${dataset}-${num_genes}.yaml"
        echo "----------------------------------------------------------------------------------------------------" 
  echo "$python -u eval_net.py --config $config_file --ref-net-file $net $algs --auc  --time --epr --stats-only"
        $python -u eval_net.py --config $config_file --ref-net-file $net $algs --auc  --time --epr --stats-only 
    done
    done
done
#echo "skipping mouse"
#exit

# ------------------------------------------------------
# ------------------------------------------------------
# mouse datasets with gene names
declare -a datasets=("Nestorowa2016-E" "Nestorowa2016-GM" "Nestorowa2016-L" "Shalek2014" "Hayashi2018")
declare -a networks=(
    #"${base_net_path}/string/mouse/10090.protein.links.v11.0-ensembl-ids-c400.csv"
    #"${base_net_path}/string/mouse/10090.protein.links.v11.0-ensembl-ids-c700.csv"
    #"${base_net_path}/TRRUST/processed/trrust_data.mouse.csv"
    #"${base_net_path}/RegNetwork/mouse-net.csv"
    #"${base_net_path}/ESCAPE/chip-x-net.csv"
    #"${base_net_path}/ESCAPE/logof-net.csv"
    "${base_net_path}/Final/mouse/mouse-net.csv"
    )

for dataset in ${datasets[*]}; do
    for num_genes in ${num_genes_list[*]}; do
    for net in ${networks[*]}; do
        config_file="config-files/2019-08-datasets-updated/$dataset/${dataset}-${num_genes}.yaml"
        echo "----------------------------------------------------------------------------------------------------" 
  echo "$python -u eval_net.py --config $config_file --ref-net-file $net $algs --auc  --time --epr --stats-only"
        $python -u eval_net.py --config $config_file --ref-net-file $net $algs --auc  --time --epr --stats-only 
    done
    done
done
