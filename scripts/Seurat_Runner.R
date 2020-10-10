####### SETUP #######
rm(list = ls())
library(dplyr)
library(data.table)
source(file = "/opt/RunMultiCCA_Changed.R")
library(Seurat)

arguments <- commandArgs(trailingOnly = TRUE)
print(arguments)
cell_file_name <- arguments[[1]]
sample_file_name <- arguments[[2]]
cell_type_target <- arguments[[3]]
markers <- strsplit(arguments[[4]], "@")[[1]]
resolutions <- as.numeric(strsplit(arguments[[5]], "@")[[1]])
output_prefix <- arguments[[6]]
output_folder <- arguments[[7]]

####### Get the data for clustering #######
cat("Reading the clustering data by sample\n")
clustering_data <- fread(cell_file_name)
clustering_data <- clustering_data[cell_type == cell_type_target]

samples <- fread(sample_file_name)
samples[is.na(comparison), comparison :=  "NA"]
sample_names <- samples[comparison != "NA", sample_name]
clustering_data <- clustering_data[Metadata_sample_name %in% sample_names, ]
sample_names <- sample_names[sample_names %in% unique(clustering_data$Metadata_sample_name)]

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
      num.cc = total_dimensions))
	while(class(seurats_integrated) == "try-error" & total_dimensions > 1)
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
	while(class(seurats_integrated) == "try-error" & total_dimensions > 1)
	{
		total_dimensions <- total_dimensions - 1	
		cat(paste0("Reduced dimensions: ", total_dimensions, "\n"))
		seurats_integrated <- try(RunCCA(seurats_by_sample[[1]], seurats_by_sample[[2]],
		genes.use = markers, num.cc = total_dimensions))
	}
  }
  reduction.type = "cca"
  cat(paste0("Aligning Subspaces, dimensions: ", total_dimensions, "\n"))
  seurats_aligned <- try(AlignSubspace(seurats_integrated, reduction.type = reduction.type,
    grouping.var = "Metadata_sample_name", verbose = F, dims.align = 1:total_dimensions,
    num.genes = max(total_dimensions / 2, 2)))
	while(class(seurats_aligned) == "try-error" & total_dimensions > 1)
	{
		total_dimensions <- total_dimensions - 1	
		cat(paste0("Aligning Subspaces, dimensions: ", total_dimensions, "\n"))
		seurats_aligned <- try(AlignSubspace(seurats_integrated, reduction.type = reduction.type,
			grouping.var = "Metadata_sample_name", verbose = F, dims.align = 1:total_dimensions,
			num.genes = max(total_dimensions / 2, 2)))
	}
}else{
  reduction.type = "pca"
  cat("Performing PCA\n")
  seurats_aligned <- RunPCA(object = seurats_by_sample[[1]], pc.genes = markers, do.print = F)
}  
rm(seurats_by_sample)

####### Align the CCA subspaces #######
cat("Finding Clusters\n")
seurats_aligned@meta.data <- as.data.table(seurats_aligned@meta.data, keep.rownames = "CellName")
seurats_aligned@meta.data[, orig.ident := NULL]
seurats_aligned <- FindClusters(seurats_aligned, reduction.type = reduction.type,
	dims.use = 1:total_dimensions, resolution = resolutions[[1]], save.SNN = TRUE, reuse.SNN = FALSE,
	edge.file.name = "edge_file_seurat.txt")
if (length(resolutions) > 1){
  identities <- lapply(resolutions[2:length(resolutions)], function(res){
  	seurats_aligned <- FindClusters(seurats_aligned, reduction.type = reduction.type,
      dims.use = 1:total_dimensions, resolution = res, save.SNN = FALSE, reuse.SNN = TRUE,
	  edge.file.name = "edge_file_seurat.txt")@ident
  })
  names(identities) <- paste0("res_", resolutions[2:length(resolutions)], "_ids")
  seurats_aligned@meta.data <- cbind(seurats_aligned@meta.data, as.data.table(identities))
}
setnames(seurats_aligned@meta.data, paste0("res.", resolutions[[1]]), paste0("res_", resolutions[[1]], "_ids"))
seurats_aligned@meta.data <- merge(seurats_aligned@meta.data, clustering_data,
	by = c("CellName", "Metadata_sample_name"))

####### Align the CCA subspaces #######
dir.create(output_folder, recursive = T, showWarnings = F)
out_name <- paste0(output_folder, "/", output_prefix, "-clusters")
save(seurats_aligned, file = paste0(out_name, ".RData"))
fwrite(seurats_aligned@meta.data, file = paste0(out_name, ".csv"))
cat("Finished!\n")
