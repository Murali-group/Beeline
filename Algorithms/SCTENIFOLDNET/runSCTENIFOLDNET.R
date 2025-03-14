library(scTenifoldNet)
library(reshape2)
args <- commandArgs(trailingOnly = T)
inFile <- args[1]
outFile <-  args[2]

# input expression data
inputExpr <- read.table(inFile, sep=",", header = 1, row.names = 1)
geneNames <- rownames(inputExpr)
rownames(inputExpr) <- c(geneNames)
inputExpr <- as.matrix(inputExpr)
#nGenes <- nrow(inputExpr)

# Run pcNet 
# Link to paper: https://doi.org/10.1101/2020.02.12.931469
set.seed(1)
pcNetResults= as.matrix(pcNet(X = inputExpr, nComp = 9))
#set.seed(1)
#pcNetResults = makeNetworks(inputExpr, nComp = round(nGenes/2), q = 0, nNet = 10)
#set.seed(1)
#pcNetResults = tensorDecomposition(pcNetResults)
#pcNetResults = as.matrix(pcNetResults$X)
diag(pcNetResults) <- 1

# Write output to a file
# https://stackoverflow.com/questions/38664241/ranking-and-counting-matrix-elements-in-r
DF = melt(pcNetResults)

#DF = data.frame(Gene1 = geneNames[c(row(pcorResults$estimate))], Gene2 = geneNames[c(col(pcorResults$estimate))]
#                , corVal = c(pcorResults$estimate), pValue =  c(pcorResults$p.value))
colnames(DF) = c('Gene1', 'Gene2', 'corVal')
outDF <- DF[order(DF$corVal, decreasing=TRUE), ]
write.table(outDF, outFile, sep = "\t", quote = FALSE, row.names = FALSE)
