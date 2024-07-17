source("ng-sem/R/loggaussian.R")
source("ng-sem/R/ng_sem_estimation.R")
source("ng-sem/R/ng_sem_expectation.R")
source("ng-sem/R/ng_sem_maximization.R")
source("ng-sem/R/ng_sem.R")

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

DF = data.frame(Gene1 = geneNames[c(row(results$w.mat))], Gene2 = geneNames[c(col(results$w.mat))]
                , weight = c(results$w.mat))
outDF <- DF[order(DF$weight, decreasing=TRUE), ]

write.table(outDF, outFile, sep = "\t", quote = FALSE, row.names = FALSE)

