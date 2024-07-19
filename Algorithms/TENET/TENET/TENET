#!/bin/bash
## specify $1 for input *.csv file, $2 for number of parallel jobs, $3 for trajectory_file , $4 for cell_select_file, $5 for historyLength.


# transpose input matrix from cell*gene to gene*cell, and generate list of all pairs of genes
cat $1 | cut -d ',' -f 2- | tail -n +2 | sed 's/,/ /g' > cell_gene.tsv
cat $1 | head -n 1 | cut -d ',' -f 2- | tr ',' '\n' > gene_names
num_gene=`cat cell_gene.tsv | wc -l | sed -e 's/^[ \t]*//'`
python ./PreProcessScript.py

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

## mpirun, by given number of jobs, check number of available cpu cores first
mpirun_cmd='time mpirun'
if [ -d "outputs" ]; then
	rm -rf outputs
fi
mkdir outputs
num_job=`grep -cv '^[[:space:]]*$' list_jobfiles`

for ((loop=1;loop<=${num_job};loop++))
do
	input_file=`cat list_jobfiles | head -n $loop | tail -n 1`
	output_id=`cat list_jobfiles | head -n $loop | tail -n 1 | cut -d '_' -f 3`
	#echo $mpirun_cmd
	mpirun_cmd=`echo $mpirun_cmd | sed -e "s/$/ -np 1 python runTE.py pair_jobs\/${input_file} outputs\/TE_out_${output_id}.csv $3 $4 $5 :/g"`
	#echo $mpirun_cmd
done
echo $mpirun_cmd | sed 's/ :$//g' > mpirun_script.sh
chmod a+x mpirun_script.sh
./mpirun_script.sh
sleep 5

cat outputs/*.csv > TE_result_all.csv
chmod a+x makeTEasMatrix.py
python makeTEasMatrix.py
rm TE_result_all.csv