####### SETUP #######
rm(list = ls())
library(data.table)
library(uwot)

source("/opt/Plot_Functions.R")

arguments <- commandArgs(trailingOnly = TRUE)                                                                              
print(arguments)                                                                                                           
cell_file_name <- arguments[[1]]                                                                                             
sample_metadata_file_name <- arguments[[2]]
category_column <- arguments[[3]]   
cell_population <- arguments[[4]]                                                                                          
markers <- strsplit(arguments[[5]], "@")[[1]]                                                                              
resolutions <- as.numeric(strsplit(arguments[[6]], "@")[[1]])                                                              
high_color <- arguments[[7]]                                                                                               
mid_color <- arguments[[8]]                                                                                                
low_color <- arguments[[9]]                                                                                                
output_folder <- arguments[[10]]   

resolutions <- paste0("res_", resolutions, "_ids")

######## Cell data #######
Cells <- fread(cell_file_name)
Cells <- Cells[cell_type == cell_population, ]
Samples <- fread(sample_metadata_file_name)
Cells <- merge(suppressWarnings(Cells[, -category_column, with = F]), Samples,
	by.x = "Metadata_sample_name", by.y = "sample_name")

######## Heatmaps #######
heatmaps <- lapply(resolutions, function(res){
  heatmapper(Cells[, c(res, markers), with = F], res_column = res, cols = markers, high_color,
	mid_color, low_color)
})
names(heatmaps) <- resolutions

################## UMAPS #####################x
set.seed(666)
UMAPS <- lapply(resolutions, function(res){
  umap_coords <- umap(Cells, n_neighbors = 40, min_dist = 0.9, learning_rate = 0.5, init = "random",
	metric = list("euclidean" = markers, "categorical" = res), n_sgd_threads = 1, n_threads = 1)
  umap_table <- cbind(umap_coords, Cells[, c(markers, res), with = F])
  setnames(umap_table, c("V1", "V2"), c("umap_x", "umap_y"))
  return(umap_table)
})
names(UMAPS) <- resolutions

###################### UMAP Plots by Cluster ##########################
UMAP_cluster_plots <- lapply(resolutions, function(res){
  label_dot_plotter(UMAPS[[res]][, c("umap_x", "umap_y", res), with = F], "umap_x", "umap_y", res)
})
names(UMAP_cluster_plots) <- resolutions

###################### UMAP Plots by Marker ##########################
umap_marker_plots <- lapply(resolutions, function(res){
  plots <- lapply(markers, function(marker){
    gradient_dot_plotter(UMAPS[[res]][, c("umap_x", "umap_y", res, marker), with = F], "umap_x", "umap_y", marker, high_color, low_color)  
  })
  names(plots) <- markers
  return(plots)
})
names(umap_marker_plots) <- resolutions


###################### Boxplots by Cluster ##########################
n_categories <- length(unique(Samples[[category_column]]))

if(n_categories == 2){
	cluster_boxplots <- lapply(resolutions, function(res){
		plot_dataset <- copy(Cells)		
		plot_dataset[, n_cells := .N, by = c("Metadata_sample_name", res)]
		plot_dataset[, total := .N, by = "Metadata_sample_name"]
		plot_dataset[, frac := n_cells / total * 100]
		cluster_boxplots <- list_boxplotter(plot_dataset, "Metadata_sample_name", res, "category", "frac")
		cluster_boxplots <- lapply(cluster_boxplots, function(x){x$Plot})
		return(cluster_boxplots)
	})
	names(cluster_boxplots) <- resolutions
}

lapply(resolutions, function(res){
	output_folder <- file.path(output_folder)
	# UMAPs by marker
	umap_marker_output_folder <- file.path(output_folder, "UMAPs_by_Marker", res)
	dir.create(umap_marker_output_folder, recursive = T, showWarnings = F)
	lapply(markers, function(marker){
		pdf_plotter(umap_marker_plots[[res]][[marker]], filename = paste0(umap_marker_output_folder,
			"/UMAP_", marker,"_", res, ".pdf"))
	})
	# UMAPs by cluster
	umap_cluster_output_folder <- file.path(output_folder, "UMAPs_by_Cluster")
	dir.create(umap_cluster_output_folder, recursive = T, showWarnings = F)
	pdf_plotter(UMAP_cluster_plots[[res]], filename = paste0(umap_cluster_output_folder, "/UMAP_cluster_", res, ".pdf"))
	# Heatmaps by cluster
	heatmap_output_folder <- file.path(output_folder, "Heatmaps")
	dir.create(heatmap_output_folder, recursive = T, showWarnings = F)
	pdf_plotter(heatmaps[[res]], filename = paste0(heatmap_output_folder, "/Heatmap_", res, ".pdf"))
	# Boxplots by cluster
	if(n_categories == 2){
		boxplot_output_folder <- file.path(output_folder, "Boxplots", res)  
		dir.create(boxplot_output_folder, recursive = T, showWarnings = F)
		clusters <- as.character(sort(unique(Cells[[res]])))
		lapply(clusters, function(cluster){
			pdf_plotter(cluster_boxplots[[res]][[cluster]], filename =  paste0(boxplot_output_folder,
				"/boxplot_", cluster,"_", res, ".pdf"))
		})
	}	
})
