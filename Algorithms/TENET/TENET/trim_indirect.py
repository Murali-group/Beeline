import sys

threshold=float(sys.argv[2])

file_name=sys.argv[1]

ifile = open(file_name)
TF=[];TFtarget=[];TFtargetTE=[];TFtargetIndirect=[]
for line in ifile:
    temp = line.split()
    if temp[0] not in TF:
        TF.append(temp[0])
        TFtarget.append([temp[2]])
        TFtargetTE.append([float(temp[1])])
        TFtargetIndirect.append([0])
    else:
        indexTemp=TF.index(temp[0])
        TFtarget[indexTemp].append(temp[2])
        TFtargetTE[indexTemp].append(float(temp[1]))
        TFtargetIndirect[indexTemp].append(0)

for i in range(len(TF)):
    for j in range(len(TFtarget[i])):
        for k in range(len(TFtarget[i])):
            if j!=k and TFtarget[i][j] in TF:
                indexTemp=TF.index(TFtarget[i][j])
                if TFtarget[i][k] in TFtarget[indexTemp]:
                    if TFtargetTE[i][k]<min(TFtargetTE[i][j],TFtargetTE[indexTemp][TFtarget[indexTemp].index(TFtarget[i][k])])+threshold:
                        TFtargetIndirect[i][k]=1

ofile = open(file_name.replace(".sif",".trimIndirect"+str(threshold)+".sif"),"w")
for i in range(len(TF)):
    for j in range(len(TFtarget[i])):
        if TFtargetIndirect[i][j]==0:
            ofile.write(TF[i]+"\t"+str(TFtargetTE[i][j])+"\t"+TFtarget[i][j]+"\n")
ofile.close()

