###### Setup ######

library(dplyr)
library(fpc)

arguments <- commandArgs(trailingOnly = TRUE)
coordinates_file_name <- arguments[[1]]
reachability_distance <- arguments[[2]]
min_cells <- arguments[[3]]
metadata_file_name <- arguments[[4]]
output_folder <- arguments[[5]]
output_filename <- arguments[[6]]

###### Read metadata file ######

metadata <- fread(metadata_file_name)
cell_types_to_be_clustered <- dplyr::filter(metadata, to_cluster == T) %>% .$cell_type %>% unique()
colours_per_clustered_cell_type <- metadata[,c("cell_type", "colour")] %>% dplyr::distinct()

###### Read input coordinates ######

coordinates <- fread(coordinates_file_name)
samples <- unique(coordinates$Metadata_sample_name)

###### Apply dbscan for each sample and cell type of interest and plot the maps ######

point_size = 0.03
stroke_val = 0.1

list_dbscan_output <- list()

for (cell_type_to_cluster in cell_types_to_be_clustered) {
  
  for (sample in samples) {
    
    out_dbscan <- dplyr::filter(coordinates, cell_type == cell_type_to_cluster & Metadata_sample_name == sample) %>% dplyr::select(c("CellName", "Metadata_sample_name", "Location_Center_X", "Location_Center_Y", "Location_Center_Z", "cell_type", "dim_x", "dim_y"))
    
    if (nrow(out_dbscan) > 0) {
      
        fpc_coordinates <- dplyr::filter(coordinates, cell_type == cell_type_to_cluster & Metadata_sample_name == sample) %>% dplyr::select(c("Location_Center_X", "Location_Center_Y"))
        dbscan_results_per_sample_per_cell_type <- fpc::dbscan(data = fpc_coordinates, e = reachability_distance, MinPts = min_cells, method = "raw")
        out_dbscan$cluster <- dbscan_results_per_sample_per_cell_type$cluster
        list_dbscan_output[[paste0(cell_type_to_cluster, "_", sample)]] <- out_dbscan
        
        out_dbscan <- merge(out_dbscan, colours_per_clustered_cell_type, by = "cell_type")
        max_x <- out_dbscan$dim_x[out_dbscan$Metadata_sample_name == sample] %>% unique()
        max_y <- out_dbscan$dim_y[out_dbscan$Metadata_sample_name == sample] %>% unique()
        
        nbr_clusters <- length(unique(out_dbscan$cluster)) - 1
        
        if (nbr_clusters > 0) {
          fviz_obj <- factoextra::fviz_cluster(object = dbscan_results_per_sample_per_cell_type, data = fpc_coordinates, ellipse.type = "convex", xlab = FALSE, ylab = FALSE, stand = FALSE, geom = "point", shape = 19, outlier.color = "black", outlier.shape = 19, pointsize = point_size, show.clust.cent = FALSE, palette = out_dbscan$colour)
        }
        else
        {
          fviz_obj <- factoextra::fviz_cluster(object = dbscan_results_per_sample_per_cell_type, data = fpc_coordinates, ellipse.type = "convex", xlab = FALSE, ylab = FALSE, stand = FALSE, geom = "point", shape = 19, outlier.color = "black", outlier.shape = 19, pointsize = point_size, show.clust.cent = FALSE)
        }
        
        fviz_obj <- fviz_obj +
          coord_fixed() +
          scale_x_continuous(position = "top", limits = c(0, max_x), expand = c(0, 0), breaks = seq(from = 0, to = max_x, length.out = 4), labels = round(seq(from = 0, to = max_x, length.out = 4), 2)) +
          scale_y_continuous(limits = c(max_y, 0), expand = c(0, 0), trans = "reverse", breaks = seq(from = max_y, to = 0, length.out = 4), labels = round(seq(from = max_y, to = 0, length.out = 4), 2)) +
          theme_bw(base_size = 4, base_family = "sans") +
          theme(plot.title = element_text(hjust = 0.5),
                panel.border = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_line(colour = "black", size = stroke_val),
                strip.background = element_rect(colour=NA, fill=NA))
        
        ggsave(filename = paste0(output_folder, "SIMPLI_DBSCAN_MAP_", sample, "_", cell_type_to_cluster, ".pdf"), device = "pdf", width = 3, height = 3, useDingbats = F)
        
    }
    
  }
  
}

output_df <- bind_rows(list_dbscan_output)
fwrite(file = output_filename, x = output_df)













