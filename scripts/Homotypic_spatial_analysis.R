###### Setup ######
library(data.table)
library(fpc)

arguments <- commandArgs(trailingOnly = TRUE)
coordinates_file_name <- arguments[[1]]
reachability_distance <- arguments[[2]]
min_cells <- arguments[[3]]
cell_type_column <- arguments[[4]]
cell_type_to_cluster <- arguments[[5]]
output_folder <- arguments[[6]]
output_filename <- arguments[[7]]

###### Read input coordinates ######
coordinates <- fread(coordinates_file_name)
setnames(coordinates, cell_type_column, "cell_type")
samples <- unique(coordinates$Metadata_sample_name)

###### Apply dbscan for each sample and cell type of interest ######
dbscan_output <- lapply(samples, function(sample){
    out_dbscan <- coordinates[cell_type == cell_type_to_cluster & Metadata_sample_name == sample, .(CellName,
		Metadata_sample_name, Location_Center_X, Location_Center_Y, Location_Center_Z, cell_type)]
    if (nrow(out_dbscan) > 0) {
        fpc_coordinates <- out_dbscan[, .(Location_Center_X, Location_Center_Y)]
        dbscan_results_per_sample <- fpc::dbscan(data = fpc_coordinates, e = reachability_distance, MinPts = min_cells, method = "raw")
        out_dbscan$cluster <- dbscan_results_per_sample$cluster
        return(out_dbscan)
	}
	return(NULL)
})  

output_df <- rbindlist(dbscan_output[!is.null(dbscan_output)])
dir.create(output_folder, recursive = T, showWarnings = F)
fwrite(file = paste0(output_folder, "/", output_filename), x = output_df)













