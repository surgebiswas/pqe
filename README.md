# pqe

Code and analysis to answer questions for my preliminary qualifying exam, chaired by Peter Karchenko.

## Progress log and workflow
### April 21, 2017
1. Read the following papers so far:
- Single-Cell Analysis Reveals a Close Relationship between Differentiating Dopamine and Subthalamic Nucleus Neuronal Lineages. Kee et al. -> Paper this exam is centered on.
- The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. Trapnell et al. -> Original Monocle paper.

2. Began analyzing Kee et al. data. 
- Downloaded the relevant count files and metadata from GEO (Accession: GSE87069). 
- Performed quality/expression filtering and DEG analysis to obtain a working gene-set for downstream analysis and modeling. See src/pdf/quality_filtering.pdf. 

### April 22, 2017
3. Read the following papers:
- Scaling single-cell genomics from phenomenology to mechanism. Tanay et al. -> Single cell review paper from Regev lab.
- Diffusion maps for high-dimensional single-cell analysis of differentiation data. Hagverdi et al. -> Diffusion map method for pseudotime estimation
- destiny: diffusion maps for large-scale single- cell data in R. Angerer et al. -> Package to do diffusion based pseudotime estimation
- MAGIC: A diffusion-based imputation method reveals gene-gene interactions in single-cell RNA-sequencing data. Dijk et al. -> Expression imputation for single-cell RNA-Seq data.

4. Started exploratory analysis.
- Examining the effects of expression imputation
- PCA on imputed vs not
- Hierarchical clustering on imputed vs not
- Expression plots of genes highlighted in Kee et al.
- These simple analyses are suggesting a different picture than what Kee et al. has proposed. For example, GFP expression is strongly anti-correlated with developmental time, as is Lmx1a expression. In the paper however, they suggest Lmx1a expression is consitent in a sub-population of cells across time, and that within the Lmx1a+ population, there is a non dopaminergic set of neurons. I don't quite see that yet, but need to do further analysis. It seems there is a bimodal population of cells at E12.5 (see Dbx1 expression plot, and PCA plots). Let's run the pseudotime algorithms to see what they reveal.

