seurat_to_matrices <- function(seurat_object, rds_path, dataset_name, output_folder, integration_object) {
  if(missing(rds_path)) rds_path <- NULL
  if(missing(seurat_object)) seurat_object <- NULL
  if(missing(dataset_name)) dataset_name <- NULL
  if(missing(integration_object)) integration_object <- FALSE

  library(Matrix)
  library(Seurat)

  if(!is.null(rds_path)) {
    seuset <- readRDS(rds_path)
  } else if (!is.null(seurat_object)) {
    seuset <- seurat_object
  } else {
    print("Provide a Seurat object to the argument seurat_object or a path with rds_path.")
  }
  
  # Create output folder
  system(paste("mkdir -p ", output_folder))
  
  # Make meta.data characters:
  seuset@meta.data <- data.frame(lapply(seuset@meta.data, as.character), stringsAsFactors = FALSE)

  ############################################################
  ## Saving the matrices needed for building a adata object ##
  ############################################################

  assay_type <- Assays(seuset)

  # Save sparse matrix for raw counts
  for (i in assay_type){
    output_name <- paste0(output_folder,"/","seurat_",as.character(i))
    print(paste0("Generating files at: ", output_folder,"/"," for assay: ", as.character(i), ", a counts.mtx and scaled_counts.csv"))
    writeMM(seuset[[i]]@counts, paste0(output_folder,"/","seurat_",i,"_counts.mtx"))
    write.csv(seuset[[i]]@scale.data, paste0(output_name,"_scaled_counts.csv"), quote = FALSE)
    # save genes and cells names
    write(x = rownames(seuset[[i]]), file = paste0(output_name, "_genes.tsv"))
    write(x = colnames(seuset[[i]]), file = paste0(output_name, "_barcodes.tsv"))
  }

  # Save metadata table, umap embeddings & HVG info
  write.csv(seuset@meta.data, paste0(output_folder, "/", "seurat_metadata.csv"), quote = FALSE)
  write.csv(seuset@reductions$umap@cell.embeddings, paste0(output_folder, "/", "seurat_umap_embedding.csv"), quote = FALSE)
  if (integration_object != TRUE){
    write.csv(seuset[[assay_type[1]]]@meta.features[,c("vst.variance.standardized", "vst.variable")], paste0(output_folder, "/", assay_type[1], "seurat_HVG.csv"), quote = FALSE)
  }
}
