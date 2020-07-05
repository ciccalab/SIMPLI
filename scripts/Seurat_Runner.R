####### SETUP #######
rm(list = ls())
library(dplyr)
library(data.table)
source(file = "/opt/RunMultiCCA_Changed.R")
library(Seurat)

arguments <- commandArgs(trailingOnly = TRUE)
print(arguments)
in_file_name <- arguments[[1]]
cell_population <- arguments[[2]]
markers <- strsplit(arguments[[3]], "@")[[1]]
resolutions <- as.numeric(strsplit(arguments[[4]], "@")[[1]])
output_prefix <- arguments[[5]]
output_folder <- arguments[[6]]

####### Get the data for clustering #######
cat("Reading the clustering data by sample\n")
clustering_data <- fread(in_file_name)
clustering_data <- clustering_data[cell_type == cell_population]
sample_names <- unique(clustering_data$Metadata_sample_name)
clustering_data_by_sample <- lapply(sample_names, function(name){
	mat <- as.matrix(transpose(clustering_data[Metadata_sample_name == name, ..markers]))
	row.names(mat) <- markers
	colnames(mat) <- clustering_data[Metadata_sample_name == name, CellName]
	return(mat)
})
names(clustering_data_by_sample) <- sample_names 

####### Create a seurat object, scale and center the data for each sample #######
cat("Scaling and centering clustering data by sample.\n")
seurats_by_sample <- lapply(sample_names, function(name){
	rat <- CreateSeuratObject(raw.data = clustering_data_by_sample[[name]], min.cells = -1, min.genes = -1, is.expr = -1)
	rat@meta.data$Metadata_sample_name <- name
	rat <- ScaleData(rat, do.scale = TRUE, do.center = TRUE,  min.cells.to.block = 10000, check.for.norm = FALSE,
		display.progress = FALSE)
	return(rat)
})
names(seurats_by_sample) <- sample_names

####### Run a CCA analysis and combine all samples #######
total_dimensions <- length(markers) - 1
if (length(sample_names) > 1){
  if(length(sample_names) > 2){
    cat("Running the multiple CCA analysis.\n")
    seurats_integrated <- try(RunMultiCCA_Changed(seurats_by_sample, genes.use = markers,
      num.cc = total_dimensions, check_duplicates = FALSE))
	while(class(seurats_integrated) == "try-error" & total_dimensions > 0)
	{
		total_dimensions <- total_dimensions - 1	
		cat(paste0("Reduced dimensions: ", total_dimensions, "\n"))
		seurats_integrated <- try(RunMultiCCA_Changed(seurats_by_sample, genes.use = markers,
		num.cc = total_dimensions))
	}
  }
  if(length(sample_names) == 2){
    cat("Running the CCA analysis.\n")
    seurats_integrated <- try(RunCCA(seurats_by_sample[[1]], seurats_by_sample[[2]],
      genes.use = markers, num.cc = total_dimensions))
	while(class(seurats_integrated) == "try-error" & total_dimensions > 0)
	{
		total_dimensions <- total_dimensions - 1	
		cat(paste0("Reduced dimensions: ", total_dimensions, "\n"))
		seurats_integrated <- try(RunCCA(seurats_by_sample[[1]], seurats_by_sample[[2]],
		genes.use = markers, num.cc = total_dimensions))
	}
  }
  reduction.type = "cca"
  cat(paste0("Aligning Subspaces, dimensions: ", total_dimensions, "\n"))
  seurats_integrated <- AlignSubspace(seurats_integrated, reduction.type = reduction.type,
    grouping.var = "Metadata_sample_name", verbose = F, dims.align = 1:total_dimensions,
    num.genes = max(total_dimensions / 2, 2))
}else{
  reduction.type = "pca"
  cat("Performing PCA\n")
  seurats_integrated <- RunPCA(object = seurats_by_sample[[1]], pc.genes = markers, do.print = F)
}  
rm(seurats_by_sample)

####### Align the CCA subspaces #######
cat("Finding Clusters\n")
seurats_integrated@meta.data <- as.data.table(seurats_integrated@meta.data, keep.rownames = "CellName")
seurats_integrated@meta.data[, orig.ident := NULL]
seurats_integrated <- FindClusters(seurats_integrated, reduction.type = reduction.type,
  dims.use = 1:total_dimensions, resolution = resolutions[[1]], save.SNN = TRUE, reuse.SNN = FALSE)
if (length(resolutions) > 1){
  identities <- lapply(resolutions[2:length(resolutions)], function(res){
  	seurats_integrated <- FindClusters(seurats_integrated, reduction.type = reduction.type,
      dims.use = 1:total_dimensions, resolution = res,	save.SNN = FALSE, reuse.SNN = TRUE)@ident
  })
  names(identities) <- paste0("res_", resolutions[2:length(resolutions)], "_ids")
  seurats_integrated@meta.data <- cbind(seurats_integrated@meta.data, as.data.table(identities))
}
setnames(seurats_integrated@meta.data, paste0("res.", resolutions[[1]]), paste0("res_", resolutions[[1]], "_ids"))
seurats_integrated@meta.data <- merge(seurats_integrated@meta.data, clustering_data,
 by = c("CellName", "Metadata_sample_name"))

####### Align the CCA subspaces #######
dir.create(output_folder, recursive = T, showWarnings = F)
out_name <- paste0(output_folder, "/", output_prefix, "_clusters")
save(seurats_integrated, file = paste0(out_name, ".RData"))
fwrite(seurats_integrated@meta.data, file = paste0(out_name, ".csv"))
cat("Finished!\n")
