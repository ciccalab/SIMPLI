###### Setup ######
library(data.table)
library(fpc)

arguments <- commandArgs(trailingOnly = TRUE)
coordinates_file_name <- arguments[[1]]
reachability_distance <- as.numeric(arguments[[2]])
min_cells <- as.numeric(arguments[[3]])
cell_type_column <- arguments[[4]]
cell_type_to_cluster <- arguments[[5]]
output_folder <- arguments[[6]]
output_filename <- arguments[[7]]

###### Read input coordinates ######
coordinates <- fread(coordinates_file_name)
setnames(coordinates, cell_type_column, "spatial_analysis_cell_type")
coordinates[, spatial_analysis_cell_type := as.character(spatial_analysis_cell_type)]
coordinates[is.na(spatial_analysis_cell_type), spatial_analysis_cell_type := "NA"]
samples <- unique(coordinates$Metadata_sample_name)

###### Apply dbscan for each sample and cell type of interest ######
dbscan_output <- lapply(samples, function(sample_name){
    out_dbscan <- coordinates[spatial_analysis_cell_type == cell_type_to_cluster & Metadata_sample_name == sample_name,
		.(CellName,	Metadata_sample_name, Location_Center_X, Location_Center_Y, spatial_analysis_cell_type)]	
	if (nrow(out_dbscan) > 0) {
        fpc_coordinates <- out_dbscan[, .(Location_Center_X, Location_Center_Y)]
        dbscan_results_per_sample <- fpc::dbscan(data = fpc_coordinates, e = reachability_distance, MinPts = min_cells, method = "raw")
        out_dbscan$cluster <- dbscan_results_per_sample$cluster
        out_dbscan$isseed <- ifelse(is.null(dbscan_results_per_sample$isseed), FALSE, dbscan_results_per_sample$isseed)
        return(out_dbscan)
	}
	return(NULL)
})  

output_df <- rbindlist(dbscan_output[!is.null(dbscan_output)])
dir.create(output_folder, recursive = T, showWarnings = F)
fwrite(file = paste0(output_folder, "/", output_filename), x = output_df)












