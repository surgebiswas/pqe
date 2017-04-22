# pqe

Code and analysis to answer my preliminary qualifying exam, chaired by Peter Karchenko.

## Progress log and workflow
### April 21, 2017
Read the following papers so far:
- Single-Cell Analysis Reveals a Close Relationship between Differentiating Dopamine and Subthalamic Nucleus Neuronal Lineages. Kee et al. -> Paper this exam is centered on.
- The dynamics and regulators of cell fate decisions are revealed by pseudotemporal ordering of single cells. Trapnell et al. -> Original Monocle paper.

Began analyzing Kee et al. data. 
- Downloaded the relevant count files and metadata from GEO (Accession: GSE87069). 
- Performed quality/expression filtering and DEG analysis to obtain a working gene-set for downstream analysis and modeling. See src/pdf/quality_filtering.pdf. 
