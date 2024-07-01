import numpy
import statsmodels.sandbox.stats.multicomp
import scipy.stats
import sys

TFlist=[]
species=sys.argv[1]
ifile = open("GO_symbol_"+species+"_regulation_of_transcription+sequence-specific_DNA_binding_list.txt")
for line in ifile:
    TFlist.append(line.replace("\n","").replace("\r",""))

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
    if gene_name[sourceIndex] in TFlist:
        temp = line.split()
        for targetIndex in range(len(temp)-1):
            if float(temp[targetIndex+1])>cutOff:
                source.append(gene_name[sourceIndex])
                TE.append(float(temp[targetIndex+1]))
                target.append(gene_name[targetIndex])
    sourceIndex=sourceIndex+1
ifile.close()

TEzscore=(TE-numpy.mean(TE))/numpy.std(TE)
TEpvalue=1-scipy.stats.norm.cdf(TEzscore)
TEfdr=statsmodels.sandbox.stats.multicomp.multipletests(TEpvalue,alpha=0.05,method='fdr_bh')

fdrCutoff=float(sys.argv[2])
ofile = open(file_name.replace(".txt",".byGRN.fdr")+str(fdrCutoff)+".sif","w")
for i in range(len(source)):
    if TEfdr[1][i]<fdrCutoff:
        ofile.write(source[i]+"\t"+str(TE[i])+"\t"+target[i]+"\n")
ofile.close()
