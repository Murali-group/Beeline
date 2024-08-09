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

nGenes <- nrow(inputExpr)
n.pc <- min(5, nGenes - 1)  # Ensure n.pc is less than number of genes

# X <- SCORPION:::makeSuperCells(inputExpr, n.pc = n.pc)
# Example adjustment for choosing between irlba and svd
# if (n.pc / nGenes > 0.1) {
#     # Use SVD
#     svd_results <- svd(inputExpr)
#     X <- svd_results$u[, 1:n.pc] %*% diag(svd_results$d[1:n.pc])
# } else {
#     # Use SCORPION makeSuperCells
#     X <- SCORPION:::makeSuperCells(inputExpr, n.pc = n.pc)
# }
if (n.pc / nGenes > 0.1) {
    # Use SVD
    svd_results <- svd(inputExpr)
    X <- svd_results$u[, 1:n.pc] %*% diag(svd_results$d[1:n.pc])
    # Retain the original row names
    rownames(X) <- rownames(inputExpr)
} else {
    # Use SCORPION makeSuperCells
    X <- SCORPION:::makeSuperCells(inputExpr, n.pc = n.pc)
    # Assuming makeSuperCells retains row names, if not, add them similarly
    rownames(X) <- rownames(inputExpr)
}


X <- cor(t(as.matrix(X)), method = 'sp')
# Write output to a file
# https://stackoverflow.com/questions/38664241/ranking-and-counting-matrix-elements-in-r
DF = melt(X)

#DF = data.frame(Gene1 = geneNames[c(row(pcorResults$estimate))], Gene2 = geneNames[c(col(pcorResults$estimate))]
#                , corVal = c(pcorResults$estimate), pValue =  c(pcorResults$p.value))
colnames(DF) = c('Gene1', 'Gene2', 'corVal')
outDF <- DF[order(DF$corVal, decreasing=TRUE), ]
write.table(outDF, outFile, sep = "\t", quote = FALSE, row.names = FALSE)
