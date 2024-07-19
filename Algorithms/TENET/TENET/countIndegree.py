import numpy
import os, glob
import sys

file_name=sys.argv[1]
ifile = open(file_name)
TFlist=[];TFlistIndegree=[]
for line in ifile:
    temp = line.split()
    if temp[2] not in TFlist:
        TFlist.append(temp[2])
        TFlistIndegree.append(1)
    else:
        TFlistIndegree[TFlist.index(temp[2])]=TFlistIndegree[TFlist.index(temp[2])]+1
TFlistIndegreeIndex=numpy.argsort(TFlistIndegree)

ofile = open(file_name+".indegree.txt","w")
for i in range(len(TFlist)):
    ofile.write(TFlist[TFlistIndegreeIndex[-i-1]]+"\t"+str(TFlistIndegree[TFlistIndegreeIndex[-i-1]])+"\n")
ofile.close()
