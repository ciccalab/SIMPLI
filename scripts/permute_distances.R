###### Setup ######
library(data.table)
library(ggplot2)

arguments <- commandArgs(trailingOnly = TRUE)
permutations <- as.integer(arguments[[1]])
distance_file_name <- arguments[[2]]
heterotypic_metadata_file_name <- arguments[[3]]
sample_metadata_file_name <- arguments[[4]]
output_folder <- arguments[[5]]

######## Cell data #######
distances <- fread(distance_file_name)
samples <- fread(sample_metadata_file_name)
metadata <- fread(heterotypic_metadata_file_name)

samples[is.na(comparison), comparison := "NA" ]
suppressWarnings(distances[, color := NULL])
suppressWarnings(distances[, comparison := NULL])
distances <- merge(distances, samples, by.x = "Metadata_sample_name", by.y = "sample_name")
distances[, pair := paste0(spatial_analysis_cell_type1, "-", spatial_analysis_cell_type2)]
distances <- distances[comparison != "NA",]
metadata[, pair := paste0(cell_type1, "-", cell_type2)]
tests <- metadata$pair

################ Shuffled Distances ###########################
shuffled <- copy(distances)
shuffled[, color := NULL]
shuffled[, comparison := NULL]

shuffled_all <- lapply(1:permutations, function(n){
    my_shuffle <- copy(shuffled)
	my_shuffle[, spatial_analysis_cell_type1 := sample(spatial_analysis_cell_type1), by = Metadata_sample_name]
    my_shuffle[, spatial_analysis_cell_type2 := sample(spatial_analysis_cell_type2), by = Metadata_sample_name]
	return(my_shuffle)
})
shuffled_all <- rbindlist(shuffled_all, idcol = "permutation")

dir.create(output_folder, recursive = T, showWarnings = F)
fwrite(file = paste0(output_folder, "/permuted_distances.csv"), x = shuffled_all)

