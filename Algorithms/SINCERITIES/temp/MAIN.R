###########################################################################
###                      SINCERITIES MAIN SCRIPT                        ###
###########################################################################

# Please prepare the data in an Excel sheet and then save it as csv file 
# using the format below 
#
# Data: s-by-m+1 matrix, where s is the total number of observations/single 
# cells and m is the number of genes. The first m columns contain 
# the expression level of each m genes, and the last column contains the
# the time-stamps.
#
# Two data formats are accepted:
# 
# A) with row header
# ------------------------------------
# Gene1  Gene2  Gene3  ... Genej  Time
#  27     80     56    ...  69      0
#  73     20     90    ...  45      0
#   .     .      .     ...  .       .
#   .     .      .     ...  .       .
#   .     .      .     ...  .       .
# ------------------------------------
# 
# B) without row header
# ------------------------------------
#  27     80     56    ...  69      0
#  73     20     90    ...  45      0
#   .     .      .     ...  .       .
#   .     .      .     ...  .       .
#   .     .      .     ...  .       .
# ------------------------------------
# 
# Please specify row header in funcion uploading. If row header is absent
# set the header argument in uploading() to false.
# 
# PACKAGES required:
# kSamples
# glmnet
# ppcor
###########################################################################

library(kSamples)
library(glmnet)
library(ppcor)

args = commandArgs(trailingOnly=TRUE)
print(args)
# *** Data loading ***
uploading <- dget("SINCERITIES functions/uploading.R")
DATA <- uploading(args[1])

# *** SINCERITIES ***

SINCERITITES <- dget("SINCERITIES functions/SINCERITIES.R")
result <- SINCERITITES(DATA,distance=1,method = 1,noDIAG = 0,SIGN = 1)

#Final ranked list of regulatory edges
adj_matrix <- result$adj_matrix/max(result$adj_matrix)
final_ranked_predictions <- dget("SINCERITIES functions/final_ranked_predictions.R")
table <- final_ranked_predictions(adj_matrix,DATA$genes,SIGN=1,saveFile = TRUE, fileNAME = args[2])
