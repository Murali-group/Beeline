suppressPackageStartupMessages(library(BTR))
suppressPackageStartupMessages(library(doParallel))


num_cores = 32 #specify the number of cores to be used.
doParallel::registerDoParallel(cores=num_cores)

args <- commandArgs(trailingOnly = T)
inFile <- args[1]
inMdl <- args[2]
maxGenesPerRule <- as.numeric(args[3])
allow_self <- as.logical(args[4])
outFile <-  args[5]
Rseed <-  args[6]

set.seed(Rseed)
cat("Random seed set to ", Rseed,"\n \n")

cat("Running BTR with input file = ", inFile, ", bmodel = ", inMdl, 
    ", max_varperrule = ", maxGenesPerRule, 
    ", self_loop = ", allow_self,"\n \n")
# input expression data
inputExpr <- read.table(inFile, sep=",", header = 1, row.names = 1)
inputExpr <- t(inputExpr)

#initial boolean model
boolInit <- read.table(inMdl, sep=",", header = 1, row.names = NULL)
boolInit = initialise_model(boolInit)


# run BTR
final_model = model_train(cdata = inputExpr, bmodel = boolInit, max_varperrule = maxGenesPerRule,
                          verbose = T, self_loop = allow_self)
#write output
writeBM(final_model, outFile, gene.names = T, rownames = T)
cat("Output written to ", outFile, "\n")