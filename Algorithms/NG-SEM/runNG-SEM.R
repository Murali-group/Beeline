library(NGSEM)

args = commandArgs(trailingOnly = T)
inFile <- args[1]
outFile <- args[2]
nk <- args[3]
miter <- args[4]
error <- args[5]
cores <- args[6]

inputExpr <- read.table(inFile, sep=",", header = 1, row.names = 1)
geneNames <- rownames(inputExpr)
rownames(inputExpr) <- c(geneNames)

results = ng_sem(as.matrix(inputExpr), nk, miter, error, cores, NULL)

