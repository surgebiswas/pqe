rm(list=ls()) # clean environment

library(destiny)
library(Biobase)
library(MASS)

setwd("~/GitHub/pqe/data")
USE_IMPUTED = FALSE

# Read in the data
if (USE_IMPUTED) {
  exp_matfile = "expression_filtered_and_DE_genes_imputed_expression_mat.txt"
  impute_word = "imputed"
} else {
  exp_matfile = "expression_filtered_and_DE_genes_expression_mat.txt"
  impute_word = "original"
}

# Read in the expression matrix
exp_mat_log_fpkm <- (as.matrix(read.table(file = exp_matfile, header=TRUE, row.names=1))) # log2(RPKM + 0.001), genes x cell
head(exp_mat_log_fpkm[,1:10])

# Phenotype and feature data
pdata <- read.table(file="expression_filtered_and_DE_genes_design_mat.txt", header=TRUE, row.names=1)
pdata$EStage <- as.double(gsub("E", "",pdata$EStage))

fdata <- data.frame(genes=rownames(exp_mat_log_fpkm))
rownames(fdata) <- rownames(exp_mat_log_fpkm)

pd <- new("AnnotatedDataFrame", data = pdata)
fd <- new("AnnotatedDataFrame", data = fdata)


exp_data <- ExpressionSet(exp_mat_log_fpkm, phenoData=pd, featureData=fd)

dm <- DiffusionMap(exp_data)

# There is some error related to "object 'linepad' not found". 
# Let's export the diffusion components directly and plot them in matlab.
ev <- eigenvectors(dm)[,1:3]
rownames(ev) <- rownames(pdata)
write.table(ev, file=paste("~/GitHub/pqe/data/diffusion_components_", impute_word, ".txt", sep=""),
            quote=FALSE,sep="\t", row.names=TRUE, col.names=FALSE)