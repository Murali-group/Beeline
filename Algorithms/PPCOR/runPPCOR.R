library(ppcor)
args <- commandArgs(trailingOnly = T)
inFile <- args[1]
outFile <-  args[2]

# input expression data
inputExpr <- read.table(inFile, sep=",", header = 1, row.names = 1)
geneNames <- rownames(inputExpr)
rownames(inputExpr) <- c(geneNames)

# Run pcor using spearman's correlation as mentioned in the PNI paper 
# Link to paper: https://www.pnas.org/content/114/23/5822

pcorResults=  pcor(x= t(as.matrix(inputExpr)), method = "spearman")

# Write output to a file
# https://stackoverflow.com/questions/38664241/ranking-and-counting-matrix-elements-in-r
DF = data.frame(Gene1 = geneNames[c(row(pcorResults$estimate))], Gene2 = geneNames[c(col(pcorResults$estimate))]
                , corVal = c(pcorResults$estimate), pValue =  c(pcorResults$p.value))
outDF <- DF[order(DF$corVal, decreasing=TRUE), ]
write.table(outDF, outFile, sep = "\t", quote = FALSE, row.names = FALSE)
