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
X=numpy.arange(num_gene-1,1,-1)
X_1=numpy.insert(X,0,1)
X_2=-numpy.insert(X,0,0)+1
bp=X_1.cumsum()-1
indx_1=numpy.zeros(shape=(bp[-1]+1,1))
indx_2=numpy.zeros(shape=(bp[-1]+1,1))
for i in range(len(bp)):
    indx_1[bp[i]]=1
    indx_2[bp[i]]=X_2[i]
indx_1=indx_1.cumsum(axis=0)
indx_2=indx_2+1
indx_2=indx_2.cumsum(axis=0)
indx=numpy.concatenate((indx_1,indx_2),axis=1)
indx=indx.astype(int)
with open("all_pairs.csv", "w") as f:
    writer = csv.writer(f)
    writer.writerows(indx)

#print("--- %s seconds ---" % (time.time() - start_time))
