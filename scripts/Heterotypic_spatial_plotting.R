###### Setup ######
library(data.table)
library(ggplot2)
source("/opt/Plot_Functions.R")

arguments <- commandArgs(trailingOnly = TRUE)
distance_file_name <- arguments[[1]]
heterotypic_metadata_file_name <- arguments[[2]]
sample_metadata_file_name <- arguments[[3]]
output_folder <- arguments[[4]]

######## Distances #######
distances <- fread(distance_file_name)
samples <- fread(sample_metadata_file_name)
metadata <- fread(heterotypic_metadata_file_name)

samples[is.na(comparison), comparison := "NA" ]
suppressWarnings(distances[, color := NULL])
suppressWarnings(distances[, comparison := NULL])
distances <- merge(distances, samples, by.x = "Metadata_sample_name", by.y = "sample_name")
distances <- distances[comparison != "NA",]
distances[, pair := paste0(spatial_analysis_cell_type1, "-", spatial_analysis_cell_type2)]

metadata[, pair := paste0(cell_type1, "-", cell_type2)]
pairs <- metadata$pair

pair_colors <- metadata[, color]
names(pair_colors) <- pairs

comps <- sort(unique(samples[!is.na(comparison), comparison]))
n_categories <- length(comps)

######## Distance plots #################
for (my_pair in pairs){
    out_dir <- paste0(output_folder, "/Distance/", my_pair)
    dir.create(out_dir, recursive = T, showWarnings = F)
    plot_data <- distances[pair == my_pair,]
    my_plot <- histogram_density(plot_data, "distance", my_pair, "Distance", "Cells", category_col = NULL,
		geom = "density", color = pair_colors[[my_pair]], fill = pair_colors[[my_pair]], alpha = 0.3, log_x = T, log_y = T,
		x_breaks = c(1, 10, 100, 500, 1000))
    pdf_plotter(filename = paste0(out_dir, "/", my_pair, "-all-heterotypic.pdf"), plot = my_plot)
    if(n_categories >= 1){
		my_plot <- histogram_density(plot_data, "distance", my_pair, "Distance", "Cells", category_col = "comparison",
			geom = "density", color = "black", fill = "grey", alpha = 0.3, log_x = T, log_y = T,
			x_breaks = c(1, 10, 100, 500, 1000))
		pdf_plotter(filename = paste0(out_dir, "/", my_pair, "-by_category-heterotypic.pdf"), plot = my_plot)
    }
}

