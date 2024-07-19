import numpy
import statsmodels.sandbox.stats.multicomp
import scipy.stats
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

TEzscore=(TE-numpy.mean(TE))/numpy.std(TE)
TEpvalue=1-scipy.stats.norm.cdf(TEzscore)
TEfdr=statsmodels.sandbox.stats.multicomp.multipletests(TEpvalue,alpha=0.05,method='fdr_bh')

fdrCutoff=float(sys.argv[1])
ofile = open(file_name.replace(".csv",".fdr")+str(fdrCutoff)+".sif","w")
for i in range(len(source)):
    if TEfdr[1][i]<fdrCutoff:
        ofile.write(gene_names[source[i]-1]+"\t"+str(TE[i])+"\t"+gene_names[target[i]-1]+"\n")
ofile.close()
