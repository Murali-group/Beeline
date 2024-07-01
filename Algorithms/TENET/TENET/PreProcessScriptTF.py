import csv
import os
import numpy
import sys
#import time
#start_time = time.time()
abspath = os.path.abspath(sys.argv[0])
dname = os.path.dirname(abspath)
os.chdir(dname)
#os.chdir("/Users/senli/Documents/Junil")
reader=csv.reader(open("cell_gene.tsv", "r"), delimiter=" ")
x=list(reader)
expression_data=numpy.array(x).astype("float")
expression_data=expression_data.T
with open("cell_gene_trsps.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerows(expression_data)
num_gene=len(expression_data)

ifile = open("gene_names")
gene_names=[]
for line in ifile:
    gene_names.append(line.replace("\n",""))
ifile.close()
ifile = open("GO_symbol_"+sys.argv[1]+"_regulation_of_transcription+sequence-specific_DNA_binding_list.txt")
TFs=[]
for line in ifile:
    TFs.append(line.replace("\n",""))
ifile.close()
TFindex=[]
for i in range(num_gene):
    if gene_names[i] in TFs:
        TFindex.append(i+1)

indx=[]
for i in TFindex:
    for j in range(num_gene):
        if i!=j+1:
            indx.append([i,j+1])
indx=numpy.array(indx)
indx=indx.astype(int)

with open("all_pairs.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerows(indx)

#print("--- %s seconds ---" % (time.time() - start_time))
