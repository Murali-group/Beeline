import os,glob

ifile = open("gene_names")
gene_name=[]
for line in ifile:
    gene_name.append(line.replace("\n","").replace("\r",""))
ifile.close()

TEmatrix=[]
for i in range(len(gene_name)):
    TEmatrixTemp=[]
    for j in range(len(gene_name)):
        TEmatrixTemp.append(0)
    TEmatrix.append(TEmatrixTemp)

ifile = open("TE_result_all.csv")
for line in ifile:
    temp=line.replace("\n","").replace("\r","").split(",")
    TEmatrix[int(temp[0])-1][int(temp[1])-1]=float(temp[2])
    if len(temp)>3:
        TEmatrix[int(temp[1])-1][int(temp[0])-1]=float(temp[3])

ofile = open("TE_result_matrix.txt","w")
ofile.write("TE")
for i in range(len(gene_name)):
    ofile.write("\t"+gene_name[i])
for i in range(len(gene_name)):
    ofile.write("\n"+gene_name[i])
    for j in range(len(gene_name)):
        ofile.write("\t"+str(TEmatrix[i][j]))
ofile.close()
