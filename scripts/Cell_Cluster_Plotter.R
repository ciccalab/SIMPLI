####### SETUP #######
rm(list = ls())
library(data.table)
library(uwot)

source("/opt/Plot_Functions.R")

arguments <- commandArgs(trailingOnly = TRUE)                                                                              
print(arguments)                                                                                                           
cell_file_name <- arguments[[1]]                                                                                             
sample_metadata_file_name <- arguments[[2]]
cell_population <- arguments[[3]]                                                                                          
markers <- strsplit(arguments[[4]], "@")[[1]]                                                                              
resolutions <- as.numeric(strsplit(arguments[[5]], "@")[[1]])                                                              
high_color <- arguments[[6]]                                                                                               
mid_color <- arguments[[7]]                                                                                                
low_color <- arguments[[8]]                                                                                                
output_folder <- arguments[[9]]   

resolutions <- paste0("res_", resolutions, "_ids")

######## Cell data #######
Cells <- fread(cell_file_name)
Samples <- fread(sample_metadata_file_name)
Samples[is.na(comparison), comparison := "NA" ]
suppressWarnings(Cells[, color := NULL])
suppressWarnings(Cells[, comparison := NULL])
Cells <- merge(Cells, Samples, by.x = "Metadata_sample_name", by.y = "sample_name")
Cells <- Cells[comparison != "NA" & cell_type == cell_population]

if(nrow(Cells) == 0)
{
	cat(paste0("No samples with comparison != NA, please check the  ", sample_metadata_file_name, " file\n."))
	quit(save = "no", status = 1)
}

######## Heatmaps #######
heatmaps <- lapply(resolutions, function(res){
	heatmapper(Cells[, c(res, markers), with = F], res_column = res, cols = markers, high_color,
		mid_color, low_color)
	})
names(heatmaps) <- resolutions

################## UMAPS #####################
set.seed(666)

checkna_warn <- function(X){  ### Override the uwot:::checkna which stops with an error!
	if (!is.null(X) && any(is.na(X))) {
		warning("Missing values found in 'X'")
    }
}	
rlang::env_unlock(env = asNamespace('uwot'))
rlang::env_binding_unlock(env = asNamespace('uwot'))
assign('checkna', checkna_warn, envir = asNamespace('uwot'))
rlang::env_binding_lock(env = asNamespace('uwot'))
rlang::env_lock(asNamespace('uwot'))

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
	label_dot_plotter(UMAPS[[res]][, c("umap_x", "umap_y", res), with = F], "umap_x", "umap_y", res, "Cluster")
})
names(umap_cluster_plots) <- resolutions

###################### UMAP Plots by Sample ##########################
umap_sample_plots <- lapply(resolutions, function(res){
	label_dot_plotter(UMAPS[[res]][, c("umap_x", "umap_y", "Metadata_sample_name"), with = F], "umap_x", "umap_y",
		"Metadata_sample_name", "Sample")
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
n_categories <- length(unique(Samples[comparison != "NA", comparison]))

sample_colors <- Samples[, color]
names(sample_colors) <- Samples[, color]

if(n_categories == 2){
	plot_dataset <- copy(Cells)		
	plot_dataset[, total := .N, by = "Metadata_sample_name"]
	comparison_boxplots <- lapply(resolutions, function(res){
		plot_dataset[, n_cells := .N, by = c("Metadata_sample_name", res)]
		plot_dataset[, frac := n_cells / total * 100]
		cluster_boxplots <- list_boxplotter(plot_dataset, "Metadata_sample_name", res, "comparison", "frac",
			"Cells in cluster / total cells %", "color", sample_colors)
		lapply(cluster_boxplots, function(x){x$Plot})
	})
	names(comparison_boxplots) <- resolutions
}

comparison_output_folder <- file.path(output_folder, "Cluster_Comparisons")
dir.create(comparison_output_folder, recursive = T, showWarnings = F)

umap_output_folder <- file.path(output_folder, "UMAPs")
dir.create(umap_output_folder, recursive = T, showWarnings = F)
for(res in resolutions){
	multi_pdf_plotter(c(list(umap_cluster_plots[[res]]), list(umap_sample_plots[[res]]), umap_marker_plots[[res]]),
		filename = file.path(umap_output_folder, paste0("UMAPs-", res, ".pdf"))) 
	pdf(file.path(comparison_output_folder, paste0(res, "-plots.pdf")),
		width = multi_pdf_w, height = multi_pdf_h, useDingbats = F)
	print(heatmaps[[res]])
	if(n_categories == 2){
		arranged_grobs <- marrangeGrob(grobs = comparison_boxplots[[res]], nrow = 2, ncol = 2, top = NULL)
		print(arranged_grobs)
	}	
	dev.off()
}

