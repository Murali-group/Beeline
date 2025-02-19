#!/bin/bash
## specify $1 for input *.csv file, $2 for number of parallel jobs, $3 for trajectory_file , $4 for branch_file.

source activate py366


## transpose input matrix from cell*gene to gene*cell, and generate list of all pairs of genes

cat $1 | cut -d ',' -f 2- | tail -n +2 | sed 's/,/ /g' > cell_gene.tsv
cat $1 | head -n 1 | cut -d ',' -f 2- | tr ',' '\n' > gene_names
num_gene=`cat cell_gene.tsv | wc -l | sed -e 's/^[ \t]*//'`
python ./PreProcessScript.py

#if [ -d "genes" ]; then
#	rm -rf genes
#fi
#mkdir genes
#cp cell_gene_trsps.csv genes/cell_gene_trsps.csv

## split by genes
#cd genes
#split -a 3 -l 1 cell_gene_trsps.csv cell_gene_
#rm -f cell_gene_trsps.csv
#cd ..
#ls -1 genes/ > list_genefiles


# split pair list into # of jobs

num_job=$2
if [ -d "pair_jobs" ]; then
	rm -rf pair_jobs
fi
mkdir pair_jobs
mv all_pairs.csv pair_jobs/all_pairs.csv
cd pair_jobs
num_pair=`cat all_pairs.csv | wc -l | sed -e 's/^[ \t]*//'`
num_line=`expr $(expr ${num_pair} / ${num_job}) + 1`
split -a 3 -l ${num_line} all_pairs.csv pair_list_
rm -f all_pairs.csv
cd ..
ls -1 pair_jobs/ > list_jobfiles


## add loops to submit jobs to the PBS cluster system

if [ -d "outputs" ]; then
	rm -rf outputs
fi
mkdir outputs
num_job=`grep -cv '^[[:space:]]*$' list_jobfiles`
#module load moab torque
for ((job_loop=1;job_loop<=${num_job};job_loop++))
do
	cat PrePBSScript > pbs_batch_job_${job_loop}
	input_file=`cat list_jobfiles | head -n ${job_loop} | tail -n 1`
	output_id=`cat list_jobfiles | head -n ${job_loop} | tail -n 1 | cut -d '_' -f 3`
	echo -e "python runTEpbsV2.py pair_jobs/${input_file} outputs/TE_out_${output_id}.csv $3 $4" >> pbs_batch_job_${job_loop}
	cat PostPBSScript >> pbs_batch_job_${job_loop}
	#qsub pbs_batch_job_${job_loop}
	echo ${job_loop}
	#sleep 0.05
done
