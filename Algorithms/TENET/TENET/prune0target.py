import sys

file_name=sys.argv[1]
ifile = open(file_name)
node_unique=[];node_outdegree=[];networks=[]
for line in ifile:
    networks.append(line)
    temp = line.split()
    if temp[0] not in node_unique:
        node_unique.append(temp[0])
        node_outdegree.append(1)
    else:
        node_outdegree[node_unique.index(temp[0])]=node_outdegree[node_unique.index(temp[0])]+1
    if temp[2] not in node_unique:
        node_unique.append(temp[2])
        node_outdegree.append(0)
ifile.close()

ofile = open(".".join(file_name.split(".")[:-1])+".prune0target.sif","w")
for line in networks:
    temp = line.split()
    if node_outdegree[node_unique.index(temp[2])]>0:
        ofile.write(line)
ofile.close()
