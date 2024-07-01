# TENET
A tool for reconstructing Transfer Entropy-based causal gene NETwork from pseudo-time ordered single cell transcriptomic data 

<div>
<img src="https://user-images.githubusercontent.com/33104430/83049138-fe23bb00-a04a-11ea-9d8a-59ca7582e759.png" width="90%"></img>
</div>

## Citation

Nucleic Acids Research, gkaa1014, https://doi.org/10.1093/nar/gkaa1014

## Dependency

	python3
	openmpi (>4.0)
	JPype

## 1. Run TENET using expression data in a csv file and pseudotime result in a text file
#### Usage

	./TENET [expression_file_name] [number_of_threads] [trajectory_file_name] [cell_select_file_name] [history_length]

#### example

	./TENET expression_data.csv 10 trajectory.txt cell_select.txt 1


#### Input 

###### (1) expression_file (raw count is recommended) - a csv file with N cells in the rows and M genes in the columns (same format with wishbone pseudotime package).

		GENE_1	GENE_2	GENE_3	...	GENE_M

	CELL_1	

	CELL_2

	CELL_3

	.
	.
	.

	CELL_N

###### (2) number_of_threads - You can use this multi-threads option. This will take lots of memory depending on the squared number of genes * the number of cells. If the program fail, you need to reduce this.

###### (3) trajectory_file - a text file of pseudotime data with N time points in the same order as the N cells of the expression file.

	0.098
	0.040
	0.023
	.
	.
	.
	0.565

###### (4) cell_select_file - a text file of cell selection data with N Boolean (1 for select and 0 for non-select) data in the same order as the N cells of the expression file.

	1
	1
	0
	.
	.
	.
	1

###### (5) history_length - the length of history. In the benchmark data TENET provides best result when the length of history set to 1.

#### Output

	TE_result_matrix.txt - TEij, M genes x M genes matrix representing the causal relationship from GENEi to GENEj.

	TE	GENE_1	GENE_2	GENE_3	...	GENE_M
	GENE_1	0	0.05	0.02	...	0.004
	GENE_2	0.01	0	0.04	...	0.12
	GENE_3	0.003	0.003	0	...	0.001
	.
	.
	.
	GENE_M	0.34	0.012	0.032	...	0

## 2. Run TENET with hdf5 file including PAGA pseudotime result
#### Usage

	./TENET4PAGAhdf5 [hdf5_file_name] [number_of_threads] [history_length] [variable_in_adata]

#### example

	./TENET4PAGAhdf5 Data.Tuck/Tuck_PAGA510genes.h5ad 10 1 X
	./TENET4PAGAhdf5 Data.Tuck/Tuck_PAGA510genes.h5ad 10 1 raw

#### Input

###### (1) hdf5 file stored after running PAGA.

###### (2) [variable_in_adata]
	If the expression matrix stored in adata.X, then choose X. If it is adata.raw.X, then choose raw.

## 3. Run TENET from TF to target using expression data in a csv file and pseudotime result in a text file
#### Usage

        ./TENET_TF [expression_file_name] [number_of_threads] [trajectory_file_name] [cell_select_file_name] [history_length] [species]

#### example

        ./TENET_TF expression_data.csv 10 trajectory.txt cell_select.txt 1 mouse

#### Input

###### (6) species - [human/mouse/rat]

#### Output

        TE_result_matrix.txt

## 4. Run TENET single core version
#### Usage

	python TENETsinglecore [expression_file_name] [trajectory_file_name] [cell_select_file_name] [history_length]

#### example

	python TENETsinglecore expression_data.csv trajectory.txt cell_select.txt 1

#### Output

	TE_result_matrix.txt

## 5. Downstream analysis

#### (1) Reconstructing GRN
###### Usage
	python makeGRN.py [cutoff for FDR]
	python makeGRNsameNumberOfLinks.py [number of links]
	python makeGRNbyTF.py [species] [cutoff for FDR]
	python makeGRNbyTFsameNumberOfLinks.py [species] [number of links]
	** Note that "TE_result_matrix.txt" should be in the same folder.

###### Example
	python makeGRN.py 0.01
	python makeGRNsameNumberOfLinks.py 1000
	python makeGRNbyTF.py human 0.01
	python makeGRNbyTFsameNumberOfLinks.py human 1000

###### Output file
	TE_result_matrix.fdr0.01.sif
	TE_result_matrix.NumberOfLinks1000.sif
	TE_result_matrix.byGRN.fdr0.01.sif
	TE_result_matrix.byGRN.NumberOflinks1000.sif

###### Parameter
	[cutoff for fdr] - A cutoff value for FDR by z-test
	[number of links] - The number of links of the GRN
	[species] - User can choose [human/mouse/rat]

#### (2) Trimming indirect edges
###### Usage
	python trim_indirect.py [name of GRN] [cutoff]
###### Example
	python trim_indirect.py TE_result_matrix.fdr0.01.sif 0
###### Output file
	TE_result_matrix.fdr0.01.trimIndirect0.0.sif
###### Parameter
	[cutoff] - A cutoff value for trimming indirect edges. Recommended range is -0.1 to 0.1

#### (3) Counting out-degree of a given GRN
###### Usage
	python countOutdegree.py [name of GRN]
###### Example
	python countOutdegree.py TE_result_matrix.fdr0.01.sif
###### Output file
	TE_result_matrix.fdr0.01.sif.outdegree.txt
