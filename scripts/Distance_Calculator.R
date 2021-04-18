###### Setup ######
library(data.table)
library(fpc)

arguments <- commandArgs(trailingOnly = TRUE)
coordinates_file_name1 <- arguments[[1]]
cell_type_column1 <- arguments[[2]]
cell_type_to_cluster1 <- arguments[[3]]
coordinates_file_name2 <- arguments[[4]]
cell_type_column2 <- arguments[[5]]
cell_type_to_cluster2 <- arguments[[6]]
output_folder <- arguments[[7]]
output_filename <- arguments[[8]]

###### Read input coordinates ######
coordinates1 <- fread(coordinates_file_name1)
setnames(coordinates1, cell_type_column1, "spatial_analysis_cell_type1")
coordinates1[, spatial_analysis_cell_type1 := as.character(spatial_analysis_cell_type1)]
coordinates1[is.na(spatial_analysis_cell_type1), spatial_analysis_cell_type1 := "NA"]

coordinates2 <- fread(coordinates_file_name2)
setnames(coordinates2, cell_type_column2, "spatial_analysis_cell_type2")
coordinates2[, spatial_analysis_cell_type2 := as.character(spatial_analysis_cell_type2)]
coordinates2[is.na(spatial_analysis_cell_type2), spatial_analysis_cell_type2 := "NA"]

coordinates1 <- coordinates1[spatial_analysis_cell_type1 == cell_type_to_cluster1, .(CellName,
  Metadata_sample_name, Location_Center_X, Location_Center_Y, spatial_analysis_cell_type1,
  id1 = paste0("1-", .I))]
coordinates2 <- coordinates2[spatial_analysis_cell_type2 == cell_type_to_cluster2, .(CellName,
  Metadata_sample_name, Location_Center_X, Location_Center_Y, spatial_analysis_cell_type2,
  id2 = paste0("2-", .I))]
  
coordinates <- merge(coordinates1, coordinates2, by = "Metadata_sample_name", allow.cartesian = T, suffixes = c("1", "2"))
coordinates[, distance := ((Location_Center_X1 - Location_Center_X2)^2 + (Location_Center_Y1 - Location_Center_Y2)^2)^0.5]

coordinates <- coordinates[, .SD[distance == min(distance)], by = id1]

coordinates[, id1 := NULL]
coordinates[, id2 := NULL]

dir.create(output_folder, recursive = T, showWarnings = F)
fwrite(file = paste0(output_folder, "/", output_filename), x = coordinates)
