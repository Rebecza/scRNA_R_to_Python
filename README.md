# scRNA_R_to_Python
Moving objects preproccesed in R (Seurat-based) to Python (AnnData/Scanpy-based).

This repository contains two simple R functions for transferring your Seurat object to Python: 
- By saving as separate matrices.
- By saving with SeuratDisk to a .h5ad file (a HDF5-based AnnData containing file).

And a python notebook for loading your Seurat-derived object to perform SCEPIA analysis.

- Further Python transfer explanation: [here](https://github.com/Rebecza/scRNA-seq/wiki/X.-Python-transfer#transfer-seurat-object)

![](https://github.com/Rebecza/scRNA-seq/blob/main/doc/SeuratS4-AnnDataH5AD_sq.png)
