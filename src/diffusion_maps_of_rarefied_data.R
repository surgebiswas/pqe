rm(list=ls()) # clean environment

library(destiny)
library(Biobase)
library(MASS)

setwd("~/GitHub/pqe/data/rarefaction")

vals = c("1", "5", "10", "30", "50", "80", "90", "100")
for (i in vals) {
  infile = paste("rarefied_", i, ".txt", sep="")
  cat(infile)
  exp_mat_log_fpkm <- log2(as.matrix(read.table(file = infile, header=TRUE, row.names=1)) + 0.001) # log2(RPKM + 0.001), genes x cell
  dm <- DiffusionMap(t(exp_mat_log_fpkm))
  ev <- eigenvectors(dm)[,1:3]
  write.table(ev, file=paste("~/GitHub/pqe/data/rarefaction/DC_", i, ".txt", sep=""),
              quote=FALSE,sep="\t", row.names=FALSE, col.names=FALSE)
  
} 


