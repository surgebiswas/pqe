rm(list=ls()) # clean environment

library(monocle)
setwd("~/GitHub/pqe/data")
USE_IMPUTED = FALSE
USE_MONOCLE2 = TRUE

# Read in the data
if (USE_IMPUTED) {
  exp_matfile = "expression_filtered_and_DE_genes_imputed_expression_mat.txt"
  impute_word = "imputed"
} else {
  exp_matfile = "expression_filtered_and_DE_genes_expression_mat.txt"
  impute_word = "original"
}

exp_mat_log_fpkm <- read.table(file = exp_matfile, header=TRUE, row.names=1)
head(exp_mat_log_fpkm[,1:10])

exp_mat_fpkm = 2^(as.matrix(exp_mat_log_fpkm)) - 0.001 # convert to FPKM


pdata <- read.table(file="expression_filtered_and_DE_genes_design_mat.txt", header=TRUE, row.names=1)
pdata$EStage <- as.double(gsub("E", "",pdata$EStage))

fdata <- data.frame(genes=rownames(exp_mat_log_fpkm))
rownames(fdata) <- rownames(exp_mat_log_fpkm)

# Created annotated data frames
# We are following a workflow that is a combination of 
# - https://bioconductor.org/packages/release/bioc/vignettes/monocle/inst/doc/monocle-vignette.pdf and
# - https://hemberg-lab.github.io/scRNA.seq.course/pseudotime-analysis.html
pd <- new("AnnotatedDataFrame", data = pdata)
fd <- new("AnnotatedDataFrame", data = fdata)


mat2use <-  exp_mat_fpkm #t ( scale(t(log2(exp_mat_fpkm + 0.001))) ); #
expFamily = tobit()
nmeth = 'log'

dCellDataSet <- monocle::newCellDataSet(mat2use, phenoData = pd, featureData = fd, 
                                     expressionFamily=expFamily)

if (USE_MONOCLE2) {
  redmeth = "DDRTree"
  method_ver = "Monocle2"
} else {
  redmeth = "ICA"
  method_ver = "Monocle"
}

dCellDataSet <- monocle::reduceDimension(dCellDataSet, norm_method=nmeth, pseudo_expr=0.001,
                                         reduction_method=redmeth)

dCellDataSet <- monocle::orderCells(dCellDataSet, reverse = FALSE)

# Plot
outplot = paste("~/GitHub/pqe/figures/trajectory_", method_ver, "_", impute_word,".png", sep="")
png(filename=outplot)
monocle::plot_cell_trajectory(dCellDataSet,color_by="EStage")
dev.off()

