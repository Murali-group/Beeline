import scanpy
import sys

adata=scanpy.read_h5ad(sys.argv[1])
if sys.argv[2]=="raw":
    data_matrix=adata.raw.X
elif sys.argv[2]=="X":
    data_matrix=adata.X

ofile = open("temp.csv","w")
for gene_name in adata.var_names:
    ofile.write(","+gene_name)
for i in range(len(adata.obs_names)):
    ofile.write("\n"+adata.obs_names[i])
    for j in range(len(adata.var_names)):
        ofile.write(","+str(data_matrix[i][j]))
ofile.close()

ofile = open("PAGApseudotime.txt","w")
for pseudotime in adata.obs['dpt_pseudotime']:
    ofile.write(str(pseudotime)+'\n')
ofile.close()

ofile = open("PAGAcell_select.txt","w")
for i in range(len(adata.obs_names)):
    ofile.write('1\n')
ofile.close()
