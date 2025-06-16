# This script is to generate expression inputs and corresponding reference network for real data by ground truth networks.

current_net_path="inputs/BEELINE-Networks/Networks"
python="python" 

declare -a num_genes_list=(
"TFs-500"
"TFs-1000"
)

# algs="--alg PIDC" #--alg GRNBOOST2 --alg GENIE3 --alg SCSGL --alg SCTENIFOLDNET --alg SCORPION"

# #-------------------------------- Non-specific ground-truth network --------------------------------#
# # human datasets with gene names
# declare -a datasets=("hESC" "hHep")
# declare -a networks=(
#     "${current_net_path}/human/Non-Specific-ChIP-seq-network.csv"
#     )

# for dataset in ${datasets[*]}; do
#     for num_genes in ${num_genes_list[*]}; do
#     for net in ${networks[*]}; do
#         config_file="config-files/scRNA-seq/Nonspecific/$dataset/${dataset}-${num_genes}-Generator.yaml"
#         echo "----------------------------------------------------------------------------------------------------" 
#   echo "$python -u generateRefNetworks.py --config $config_file --ref-net-file $net"
#         $python -u generateRefNetworks.py --config $config_file --ref-net-file $net 
#     done
#     done
# done

# # mouse datasets with gene names
# declare -a datasets=("mDC" "mESC" "mHSC-E" "mHSC-GM" "mHSC-L")
# declare -a networks=(
#     "${current_net_path}/mouse/Non-Specific-ChIP-seq-network.csv"
#     )

# for dataset in ${datasets[*]}; do
#     for num_genes in ${num_genes_list[*]}; do
#     for net in ${networks[*]}; do
#         config_file="config-files/scRNA-seq/Nonspecific/$dataset/${dataset}-${num_genes}-Generator.yaml"
#         echo "----------------------------------------------------------------------------------------------------" 
#   echo "$python -u generateRefNetworks.py --config $config_file --ref-net-file $net"
#         $python -u generateRefNetworks.py --config $config_file --ref-net-file $net 
#     done
#     done
# done

# #-------------------------------- String ground-truth network --------------------------------#
# # human datasets with gene names
# declare -a datasets=("hESC" "hHep")
# declare -a networks=(
#     "${current_net_path}/human/STRING-network.csv"
#     )

# for dataset in ${datasets[*]}; do
#     for num_genes in ${num_genes_list[*]}; do
#     for net in ${networks[*]}; do
#         config_file="config-files/scRNA-seq/String/$dataset/${dataset}-${num_genes}-Generator.yaml"
#         echo "----------------------------------------------------------------------------------------------------" 
#   echo "$python -u generateRefNetworks.py --config $config_file --ref-net-file $net"
#         $python -u generateRefNetworks.py --config $config_file --ref-net-file $net 
#     done
#     done
# done

# # mouse datasets with gene names
# declare -a datasets=("mDC" "mESC" "mHSC-E" "mHSC-GM" "mHSC-L")
# declare -a networks=(
#     "${current_net_path}/mouse/STRING-network.csv"
#     )

# for dataset in ${datasets[*]}; do
#     for num_genes in ${num_genes_list[*]}; do
#     for net in ${networks[*]}; do
#         config_file="config-files/scRNA-seq/String/$dataset/${dataset}-${num_genes}-Generator.yaml"
#         echo "----------------------------------------------------------------------------------------------------" 
#   echo "$python -u generateRefNetworks.py --config $config_file --ref-net-file $net"
#         $python -u generateRefNetworks.py --config $config_file --ref-net-file $net 
#     done
#     done
# done


#-------------------------------- Cell type-specific ground-truth network --------------------------------#
# # human datasets with gene names
# declare -a datasets=("hESC")
# declare -a networks=(
#     "${current_net_path}/human/hESC-ChIP-seq-network.csv"
#     )

# for dataset in ${datasets[*]}; do
#     for num_genes in ${num_genes_list[*]}; do
#     for net in ${networks[*]}; do
#         config_file="config-files/scRNA-seq/Specific/$dataset/${dataset}-${num_genes}-Generator.yaml"
#         echo "----------------------------------------------------------------------------------------------------" 
#   echo "$python -u generateRefNetworks.py --config $config_file --ref-net-file $net"
#         $python -u generateRefNetworks.py --config $config_file --ref-net-file $net 
#     done
#     done
# done

# declare -a datasets=("hHep")
# declare -a networks=(
#     "${current_net_path}/human/HepG2-ChIP-seq-network.csv"
#     )

# for dataset in ${datasets[*]}; do
#     for num_genes in ${num_genes_list[*]}; do
#     for net in ${networks[*]}; do
#         config_file="config-files/scRNA-seq/Specific/$dataset/${dataset}-${num_genes}-Generator.yaml"
#         echo "----------------------------------------------------------------------------------------------------" 
#   echo "$python -u generateRefNetworks.py --config $config_file --ref-net-file $net"
#         $python -u generateRefNetworks.py --config $config_file --ref-net-file $net 
#     done
#     done
# done

# # mouse datasets with gene names
# declare -a datasets=("mDC")
# declare -a networks=(
#     "${current_net_path}/mouse/mDC-ChIP-seq-network.csv"
#     )

# for dataset in ${datasets[*]}; do
#     for num_genes in ${num_genes_list[*]}; do
#     for net in ${networks[*]}; do
#         config_file="config-files/scRNA-seq/Specific/$dataset/${dataset}-${num_genes}-Generator.yaml"
#         echo "----------------------------------------------------------------------------------------------------" 
#   echo "$python -u generateRefNetworks.py --config $config_file --ref-net-file $net"
#         $python -u generateRefNetworks.py --config $config_file --ref-net-file $net 
#     done
#     done
# done

# declare -a datasets=("mESC")
# declare -a networks=(
#     "${current_net_path}/mouse/mESC-ChIP-seq-network.csv"
#     )

# for dataset in ${datasets[*]}; do
#     for num_genes in ${num_genes_list[*]}; do
#     for net in ${networks[*]}; do
#         config_file="config-files/scRNA-seq/Specific/$dataset/${dataset}-${num_genes}-Generator.yaml"
#         echo "----------------------------------------------------------------------------------------------------" 
#   echo "$python -u generateRefNetworks.py --config $config_file --ref-net-file $net"
#         $python -u generateRefNetworks.py --config $config_file --ref-net-file $net 
#     done
#     done
# done


declare -a datasets=("mESC")
declare -a networks=(
    "${current_net_path}/mouse/mESC-lofgof-network.csv"
    )

for dataset in ${datasets[*]}; do
    for num_genes in ${num_genes_list[*]}; do
    for net in ${networks[*]}; do
        config_file="config-files/scRNA-seq/lofgof/$dataset/${dataset}-${num_genes}-Generator.yaml"
        echo "----------------------------------------------------------------------------------------------------" 
  echo "$python -u generateRefNetworks.py --config $config_file --ref-net-file $net"
        $python -u generateRefNetworks.py --config $config_file --ref-net-file $net 
    done
    done
done

# declare -a datasets=("mHSC-E")
# declare -a networks=(
#     "${current_net_path}/mouse/mHSC-ChIP-seq-network.csv"
#     )

# for dataset in ${datasets[*]}; do
#     for num_genes in ${num_genes_list[*]}; do
#     for net in ${networks[*]}; do
#         config_file="config-files/scRNA-seq/Specific/$dataset/${dataset}-${num_genes}-Generator.yaml"
#         echo "----------------------------------------------------------------------------------------------------" 
#   echo "$python -u generateRefNetworks.py --config $config_file --ref-net-file $net"
#         $python -u generateRefNetworks.py --config $config_file --ref-net-file $net 
#     done
#     done
# done

# declare -a datasets=("mHSC-GM")
# declare -a networks=(
#     "${current_net_path}/mouse/mHSC-ChIP-seq-network.csv"
#     )

# for dataset in ${datasets[*]}; do
#     for num_genes in ${num_genes_list[*]}; do
#     for net in ${networks[*]}; do
#         config_file="config-files/scRNA-seq/Specific/$dataset/${dataset}-${num_genes}-Generator.yaml"
#         echo "----------------------------------------------------------------------------------------------------" 
#   echo "$python -u generateRefNetworks.py --config $config_file --ref-net-file $net"
#         $python -u generateRefNetworks.py --config $config_file --ref-net-file $net 
#     done
#     done
# done

# declare -a datasets=("mHSC-L")
# declare -a networks=(
#     "${current_net_path}/mouse/mHSC-ChIP-seq-network.csv"
#     )

# for dataset in ${datasets[*]}; do
#     for num_genes in ${num_genes_list[*]}; do
#     for net in ${networks[*]}; do
#         config_file="config-files/scRNA-seq/Specific/$dataset/${dataset}-${num_genes}-Generator.yaml"
#         echo "----------------------------------------------------------------------------------------------------" 
#   echo "$python -u generateRefNetworks.py --config $config_file --ref-net-file $net"
#         $python -u generateRefNetworks.py --config $config_file --ref-net-file $net 
#     done
#     done
# done