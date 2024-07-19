import numpy
import sys

ifile = open("gene_names")
gene_names=[]
for line in ifile:
    gene_names.append(line.replace("\n",""))
ifile.close()

file_name="TE_result_all.csv"
ifile = open(file_name)
cutOff=0
source=[]
TE=[]
target=[]
for line in ifile:
    temp = line.replace("\n","").split(",")
    if float(temp[1])>cutOff:
        source.append(int(temp[0]))
        TE.append(float(temp[2]))
        target.append(int(temp[1]))
ifile.close()

TE=numpy.array(TE)
TEsortIndex=numpy.argsort(TE)

NumberOfLinks=sys.argv[1]
ofile = open(file_name.replace(".csv",".NumberOfLinks")+NumberOfLinks+".sif","w")
for i in range(int(NumberOfLinks)):
    ofile.write(gene_names[source[TEsortIndex[-i-1]]-1]+"\t"+str(TE[-i-1])+"\t"+gene_names[target[TEsortIndex[-i-1]]-1]+"\n")
ofile.close()



