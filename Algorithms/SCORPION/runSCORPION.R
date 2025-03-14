library(SCORPION)
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

# Run SCORPION

X <- SCORPION:::makeSuperCells(inputExpr, n.pc = 5)
X <- cor(t(as.matrix(X)), method = 'sp')
# Write output to a file
# https://stackoverflow.com/questions/38664241/ranking-and-counting-matrix-elements-in-r
DF = melt(X)

#DF = data.frame(Gene1 = geneNames[c(row(pcorResults$estimate))], Gene2 = geneNames[c(col(pcorResults$estimate))]
#                , corVal = c(pcorResults$estimate), pValue =  c(pcorResults$p.value))
colnames(DF) = c('Gene1', 'Gene2', 'corVal')
outDF <- DF[order(DF$corVal, decreasing=TRUE), ]
write.table(outDF, outFile, sep = "\t", quote = FALSE, row.names = FALSE)
