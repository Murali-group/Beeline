import numpy
import sys

file_name="TE_result_matrix.txt"

ifile = open(file_name)
line = ifile.readline()
temp = line.split()
gene_name=[]
for i in range(len(temp)-1):
    gene_name.append(temp[i+1])

cutOff=0
sourceIndex=0
TEnetwork=[]
source=[]
TE=[]
target=[]
for line in ifile:
    temp = line.split()
    for targetIndex in range(len(temp)-1):
        if float(temp[targetIndex+1])>cutOff:            
            source.append(gene_name[sourceIndex])
            TE.append(float(temp[targetIndex+1]))
            target.append(gene_name[targetIndex])
    sourceIndex=sourceIndex+1
ifile.close()

TE=numpy.array(TE)

TEsortIndex=numpy.argsort(TE)

NumberOfLinks=sys.argv[1]
ofile = open(file_name.replace(".txt",".NumberOfLinks")+NumberOfLinks+".sif","w")
for i in range(int(NumberOfLinks)):
    ofile.write(source[TEsortIndex[-i-1]]+"\t"+str(TE[-i-1])+"\t"+target[TEsortIndex[-i-1]]+"\n")
ofile.close()



