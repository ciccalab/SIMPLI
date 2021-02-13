###### Setup ######
library(data.table)
library(fpc)
library(ggplot2)
library(ggpubr)
source("/opt/Plot_Functions.R")

arguments <- commandArgs(trailingOnly = TRUE)
dbscan_file_name <- arguments[[1]]
metadata_file_name <- arguments[[2]]
output_folder <- arguments[[3]]
cell_mask_file_list <- unlist(arguments[4:length(arguments)])

###### Read metadata file ######
metadata <- fread(metadata_file_name)
colors_per_clustered_cell_type <- unique(metadata[, .(spatial_analysis_cell_type = cell_type_to_cluster, color)])

###### Read dbscan clusters ######
dbscan_clusters <- fread(dbscan_file_name)
samples <- unique(dbscan_clusters$Metadata_sample_name)
cell_types <- colors_per_clustered_cell_type$spatial_analysis_cell_type

dbscan_clusters <- merge(dbscan_clusters, colors_per_clustered_cell_type, by = "spatial_analysis_cell_type")

###### Apply dbscan for each sample and cell type of interest and plot the maps ######
cell_mask_sizes <- lapply(samples, function(sample_name){
	cell_mask_file_name <- grep(sample_name, cell_mask_file_list, value = T)
	print(cell_mask_file_name)
	print(sample_name)
	cell_mask <- load_image(cell_mask_file_name)
	list(width = dim(cell_mask)[1], height = dim(cell_mask)[2])
})
names(cell_mask_sizes) <- samples

for (cell_type in cell_types){
	out_dir <- paste0(output_folder, "/", cell_type)
	dir.create(out_dir, recursive = T, showWarnings = F)
	for (sample_name in samples){
		my_dbscan <- dbscan_clusters[Metadata_sample_name == sample_name & spatial_analysis_cell_type == cell_type,
			.(Location_Center_X, Location_Center_Y, cluster, isseed, color)]
		outliers <- my_dbscan[cluster == 0]
		my_dbscan <- my_dbscan[cluster > 0]
		my_dbscan[, cluster := as.factor(cluster)]
		max_x <- cell_mask_sizes[[sample_name]]$width
		max_y <- cell_mask_sizes[[sample_name]]$height
		if (nrow(my_dbscan) > 0){
			my_plot <- ggscatter(data = my_dbscan, x = "Location_Center_X", y = "Location_Center_Y", merge = FALSE,
				palette = my_dbscan$color, color = "cluster", fill = "cluster", combine = FALSE, shape = 19, size = 0.3,
				point = TRUE, rug = FALSE, title = paste0(cell_type, " ", sample_name), xlab = NULL, ylab = NULL,
				ellipse = TRUE, ellipse.level = 0.95, ellipse.type = "convex", ellipse.alpha = 0.2,
				ellipse.border.remove = FALSE)
		} else {my_plot <- ggplot()}
		my_plot <- my_plot + geom_point(data = outliers, aes(Location_Center_X, Location_Center_Y), size = 0.3,
				color = "black", shape = 19) + coord_fixed() +
			scale_x_continuous(position = "top", limits = c(0, max_x), expand = c(0, 0), breaks = seq(from = 0,
				to = max_x, length.out = 4), labels = round(seq(from = 0, to = max_x, length.out = 4), 2)) +
			scale_y_continuous(limits = c(max_y, 0), expand = c(0, 0), trans = "reverse", breaks = seq(from = max_y,
				to = 0, length.out = 4), labels = round(seq(from = max_y, to = 0, length.out = 4), 2)) +
			theme_bw(base_size = 4, base_family = "sans") +
			theme(plot.title = element_text(hjust = 0.5), panel.border = element_blank(), legend.position = "none",
				axis.title.x = element_blank(),	axis.title.y = element_blank(),	panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),	panel.background = element_blank(),
				axis.line = element_line(colour = "black", size = 0.25),
				strip.background = element_rect(colour = NA, fill = NA))
		pdf_plotter(filename = paste0(out_dir, "/", cell_type, "-", sample_name, "-homotypic.pdf"), plot = my_plot)
	}
}

