seurat_to_python <- function(seurat_object, rds_path, dataset_name, output_folder, velocity_object, integration_object, h5ad_object, h5ad_scaled_counts) {
  if(missing(rds_path)) rds_path <- NULL
  if(missing(seurat_object)) seurat_object <- NULL
  if(missing(dataset_name)) dataset_name <- NULL
  if(missing(h5ad_object)) h5ad_object <- FALSE
  if(missing(velocity_object)) velocity_object <- FALSE
  if(missing(integration_object)) integration_object <- FALSE
  if(missing(h5ad_scaled_counts)) h5ad_scaled_counts <- TRUE

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

  # Set assay names per Seurat object type
  if (velocity_object == TRUE) {
    assay_type <- c("sf", "uf")
  } else if (integration_object == TRUE) {
    assay_type <- Assays(seuset)
  } else {
    assay_type <- c("RNA")
  }

  ############################################################
  ## Saving the matrices needed for building a adata object ##
  ############################################################
  # Save sparse matrix for raw counts
  for(i in assay_type){
    print(paste0("Generating files at: ", output_folder,"/"," for assay: ",i, ", a counts.mtx and scaled_counts.csv"))
    writeMM(seuset[[i]]@counts, paste0(output_folder,"/","seurat_",i,"_counts.mtx"))
    write.csv(seuset[[i]]@scale.data, paste0(output_folder,"/","seurat_",i,"_scaled_counts.csv"), quote = FALSE)
    # save genes and cells names
    write(x = rownames(seuset[[i]]), file = paste0(output_folder, "/", i, "_genes.tsv"))
    write(x = colnames(seuset[[i]]), file = paste0(output_folder, "/", i, "_barcodes.tsv"))
  }
  

  # Save metadata table, umap embeddings & HVG info
  write.csv(seuset@meta.data, paste0(output_folder, "/","seurat_metadata.csv"), quote = FALSE)
  write.csv(seuset@reductions$umap@cell.embeddings, paste0(output_folder,"/","seurat_umap_embedding.csv"), quote = FALSE)
  if(integration_object != TRUE){
    write.csv(seuset[[i[1]]]@meta.features[,c("vst.variance.standardized", "vst.variable")], paste0(output_folder,"/","seurat_HVG.csv"), quote = FALSE)
  }
  
  
  #############################################################
  ## Transfer to H5AD format -- Seurat-Disk package needed!! ##
  #############################################################
  # only runs when h5ad_object == TRUE

  if(h5ad_object == TRUE){

    ## Check if SeuratDisk is installed before running:
    if(!require("SeuratDisk")){

      print("Make sure you install SeuratDisk to perform the H5AD transfer!")
      print("Install mojaveazure/loomR from devtools, remotes and mojaveazure/seurat-disk. For instructions: https://github.com/mojaveazure/seurat-disk.")

    }else{

      library(SeuratDisk)

      ## Obtain name for dataset
      if (is.null(dataset_name)){ 
        # If no dataset_name is specified, take from rds_path:
        if (!is.null(rds_path)){
          dataset_name <- gsub(".rds", "", gsub("^.*\\/", "", rds_path))
        } else {
          dataset_name <- "seuset"
        }
        
      }

      ## Velocity object assay name adjustments
      if (velocity_object == TRUE){
        # Rename assays for scvelo:
        seuset <- RenameAssays(seuset, sf = "spliced")
        seuset <- RenameAssays(seuset, uf = "unspliced")

        # Save an RNA assay for adata.X:
        seuset[["RNA"]] <- seuset[["spliced"]]
        DefaultAssay(seuset) <- "RNA"
      }
      
      ## Object after Seurat integration
      if (integration_object == TRUE){
        DefaultAssay(seuset) <- "integrated"
      }

      ## Optional: Delete the scaled data for transfer to Python (normalized counts will be transferred):
      if (h5ad_scaled_counts == FALSE){
        seuset <- DietSeurat(seuset, counts = TRUE,
                 data = TRUE,
                 scale.data = FALSE,
                 dimreducs = c("pca", "umap"),
                 graphs = c(paste0(assay_type[1],"_nn"), paste0(assay_type[1],"_snn")))
      } else if (h5ad_scaled_counts == TRUE) {
        dataset_name <- paste0(dataset_name, "_scaled")
      }

      print("Conversion from .h5Seurat to .h5ad.")
      ## Convert to h5ad:
      h5Seurat_path <- paste0(output_folder, "/", dataset_name, ".h5Seurat" )
      SaveH5Seurat(seuset, filename = h5Seurat_path)
      Convert(h5Seurat_path, dest = "h5ad")
      
      print(paste0("H5AD object saved at: ", h5Seurat_path))
      
      ## Print structure of h5ad object:
      hfile <- Connect(h5Seurat_path)
      print(hfile$index())
      print(hfile[["assays"]])
      print(hfile[["reductions"]])
      hfile$close_all()

    }

  }

}


