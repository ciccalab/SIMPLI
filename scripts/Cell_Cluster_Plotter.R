####### SETUP #######
rm(list = ls())
library(data.table)
library(uwot)

source("/opt/Plot_Functions.R")

arguments <- commandArgs(trailingOnly = TRUE)                                                                              
print(arguments)                                                                                                           
cell_file_name <- arguments[[1]]                                                                                             
sample_metadata_file_name <- arguments[[2]]
comparison_columns <- strsplit(arguments[[3]], ",")[[1]]   
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
Samples <- fread(sample_metadata_file_name)
Samples[, c(comparison_columns) :=  lapply(.SD, function(x){ifelse(is.na(x), "NA", x)}), .SDcols = comparison_columns]
suppressWarnings(Cells[, color := NULL])
suppressWarnings(Cells[, comparison_columns := NULL])
Cells <- merge(Cells, Samples, by.x = "Metadata_sample_name", by.y = "sample_name")

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
umap_table <- cbind(umap_coords, Cells[, c("Metadata_sample_name", markers, res), with = F])
setnames(umap_table, c("V1", "V2"), c("umap_x", "umap_y"))
return(umap_table)
})
names(UMAPS) <- resolutions

###################### UMAP Plots by Cluster ##########################
umap_cluster_plots <- lapply(resolutions, function(res){
	label_dot_plotter(UMAPS[[res]][, c("umap_x", "umap_y", res), with = F], "umap_x", "umap_y", res)
})
names(umap_cluster_plots) <- resolutions

###################### UMAP Plots by Sample ##########################
umap_sample_plots <- lapply(resolutions, function(res){
	label_dot_plotter(UMAPS[[res]][, c("umap_x", "umap_y", "Metadata_sample_name"), with = F], "umap_x", "umap_y", "Metadata_sample_name")
})
names(umap_sample_plots) <- resolutions

###################### UMAP Plots by Marker ##########################
umap_marker_plots <- lapply(resolutions, function(res){
	plots <- lapply(markers, function(marker){
		gradient_dot_plotter(UMAPS[[res]][, c("umap_x", "umap_y", res, marker), with = F],
			"umap_x", "umap_y", marker, high_color, low_color)  
	})
	names(plots) <- markers
	return(plots)
})
names(umap_marker_plots) <- resolutions

###################### Boxplots by Cluster ##########################
n_categories <- sapply(comparison_columns, function(comparison_column){
	sum(unique(Samples[[comparison_column]]) != "NA")
})
names(n_categories) <- comparison_columns

sample_colors <- Samples[, color]
names(sample_colors) <- Samples[, color]

comparison_boxplots <- lapply(comparison_columns, function(comparison_column){
	all_boxplots <- list()
	if(n_categories[[comparison_column]] == 2){
		all_boxplots <- lapply(resolutions, function(res){
			plot_dataset <- copy(Cells)		
			plot_dataset <- plot_dataset[plot_dataset[[comparison_column]] != "NA"]
			plot_dataset[, n_cells := .N, by = c("Metadata_sample_name", res)]
			plot_dataset[, total := .N, by = "Metadata_sample_name"]
			plot_dataset[, frac := n_cells / total * 100]
			cluster_boxplots <- list_boxplotter(plot_dataset, "Metadata_sample_name", res, comparison_column, "frac",
				"Cells in cluster / total cells %", "color", sample_colors)
			lapply(cluster_boxplots, function(x){x$Plot})
		})
		names(all_boxplots) <- resolutions
	}
	return(all_boxplots)
})
names(comparison_boxplots) <- comparison_columns 

comparison_output_folder <- file.path(output_folder, "Cluster_Comparisons")
dir.create(comparison_output_folder, recursive = T, showWarnings = F)

for(comparison_column in comparison_columns){
	umap_output_folder <- file.path(output_folder, "UMAPs")
	dir.create(umap_output_folder, recursive = T, showWarnings = F)
	for(res in resolutions){
		multi_pdf_plotter(c(list(umap_cluster_plots[[res]]), list(umap_sample_plots[[res]]), umap_marker_plots[[res]]),
			filename = file.path(umap_output_folder, paste0("UMAPs-", res, ".pdf"))) 
		if(n_categories[[comparison_column]] == 2){
			pdf(NULL)
			arranged_grobs <- marrangeGrob(grobs = comparison_boxplots[[comparison_column]][[res]],
				nrow = 3, ncol = 2, top = NULL)
			dev.off()
			pdf(file.path(comparison_output_folder, paste0(comparison_column, "-", res, "-plots.pdf")),
				width = multi_pdf_w, height = multi_pdf_h, useDingbats = F)
			print(heatmaps[[res]])
			print(arranged_grobs)
			dev.off()
		}	
	}
}	

if (all(n_categories != 2)){
	for(res in resolutions){
		pdf_plotter(file.path(comparison_output_folder, paste0("Heatmap-", res, ".pdf")), heatmaps[[res]],
			width = single_pdf_w, height = single_pdf_h)
	}
}
