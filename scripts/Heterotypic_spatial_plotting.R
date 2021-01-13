###### Setup ######
library(data.table)
library(ggplot2)
source("/opt/Plot_Functions.R")

arguments <- commandArgs(trailingOnly = TRUE)
distance_file_name <- arguments[[1]]
heterotypic_metadata_file_name <- arguments[[2]]
sample_metadata_file_name <- arguments[[3]]
output_folder <- arguments[[4]]

######## Cell data #######
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
tests <- metadata$pair

test_colors <- metadata[, color]
names(test_colors) <- tests

n_categories <- length(unique(samples[comparison != "NA", comparison]))

for (test in tests){
    out_dir <- paste0(output_folder, "/", test)
    dir.create(out_dir, recursive = T, showWarnings = F)
    plot_data <- distances[pair == test ,]
    my_plot <- histogram_density(plot_data, "distance", test, "Distance", "Cells", category_col = NULL,
		geom = "density", color = test_colors[[test]],fill = test_colors[[test]], alpha = 0.3, log_x = T, log_y = T,
		x_breaks = c(1, 10, 100, 500, 1000))
    pdf_plotter(filename = paste0(out_dir, "/", test, "-all-heterotypic.pdf"), plot = my_plot)
    if(n_categories >= 1){
      my_plot <- histogram_density(plot_data, "distance", test, "Distance", "Cells", category_col = "comparison",
		geom = "density", color = "black", fill = "grey", alpha = 0.3, log_x = T, log_y = T,
		x_breaks = c(1, 10, 100, 500, 1000))
      pdf_plotter(filename = paste0(out_dir, "/", test, "-by_category-heterotypic.pdf"), plot = my_plot)
    }
}

