####### SETUP #######
rm(list = ls())
library(dplyr)
library(tidyr)
library(data.table)
library(EBImage)

args <- commandArgs(trailingOnly = TRUE)                                                                                
unannotated_cells_file_name <- args[[1]]                                                                                
threshold_metadata_file_name <- args[[2]]
image_metadata_file_name <- args[[3]]                                                                               
cell_mask_metadata_file_name <- args[[4]]                                                                               

out_file_name <- args[[5]]                                                                                              

############ Load the unannotated cells ##########                                                                      
cells <- fread(unannotated_cells_file_name)                                                                             
cells[, CellName := paste0(Metadata_sample_name, "-", ObjectNumber)]

############################# Marker combinationsi and Thresholds ###########################                                           
# Expects a .csv file with 3 columns:                                                                                   
# - cell_type                                                                                                           
# - threshold_marker                                                                                                    
# - threshold_value                                                                                                     
threshold_metadata <- fread(threshold_metadata_file_name)                                                               
cell_types <- threshold_metadata[, cell_type]

############################# Image Metadata ###########################                                           
# Expects a .csv file with 3 columns:                                                                                      
# - sample_name                                                                                                            
# - label                                                                                                                  
# - file_name                                                                                                              
image_metadata <- fread(image_metadata_file_name)                                                                                
all_samples <- unique(image_metadata$sample_name) 

############################# Cell Mask Metadata ###########################                                           
# Expects a .csv file with 3 columns:                                                                                      
# - sample_name                                                                                                            
# - label = Cell_Mask                                                                                                                  
# - file_name 
cell_mask_metadata <- fread(cell_mask_metadata_file_name)
image_metadata <- rbind(cell_mask_metadata, image_metadata)
rm(cell_mask_metadata)

############################# Cell Selection ###########################                                                
load_image <- function(filename)
{
  Img <- readImage(filename, "tiff", all = FALSE)
  if (length(dim(Img)) == 3) {return(Img[,,1])}
  return(Img)
}

outside_remover <- function(objs_mask, mask, fraction)
{
  objs <- imageData(objs_mask)
  to_remove <- data.table(table(as.numeric(objs[which(objs & !imageData(mask), arr.ind = TRUE)])))
  all <- data.table(table(as.numeric(objs)))
  to_remove <- merge(to_remove, all, by = "V1", all = F)
  to_remove <- to_remove[N.x / N.y > fraction]$V1
  to_remove <- which(objs %in% to_remove, arr.ind = T)
  objs[to_remove] <- 0
  as.Image(objs)
}

process_sample <- function(name, fraction, marker_expression, cell_mask)
{ 
  markers <- as.character(unlist(sapply(sapply(marker_expression, function(expr){parse(text = expr)}), all.vars)))
  Marker_images <- lapply(markers, function(marker){
    load_image(image_metadata[sample_name == name & label == marker, file_name])})
  Marker_images <- lapply(Marker_images, function(x){x & x})
  names(Marker_images) <- markers
  Marker_image <- as.Image(eval(parse(text = marker_expression), envir = Marker_images))
  Cell_image <- load_image(image_metadata[sample_name == name & label == cell_mask, file_name])
  Cell_image <- Cell_image * (2^16 - 1)
  Marker_Cell_mask <- outside_remover(Cell_image, Marker_image, fraction)
  return(list(All_cells = unique(as.numeric(Cell_image@.Data)), Marker_Cells = unique(as.numeric(Marker_Cell_mask@.Data))))
}


############# Identify cells #############
for (name in all_samples){
  for (type_n in seq_along(cell_types)){
    marker <- threshold_metadata[cell_type == cell_types[[type_n]], threshold_marker]
    fraction <- threshold_metadata[cell_type == cell_types[[type_n]], threshold_value]
    type_cells <- process_sample(name, fraction, marker, "Cell_Mask")[[2]]
    type_cells <- type_cells[type_cells != 0]
	type_cells <- paste0(name, "-", type_cells)
    cells[Metadata_sample_name == name, c(cell_types[[type_n]]) := (CellName %in% unlist(type_cells)) * type_n]
  }
}

############# Prioritize Populations #############
print(cell_types)
cells[, min_cell_type := apply(.SD, 1, function(x){ifelse(any(x) > 0, min(x[x > 0]), NA)}), .SDcols = cell_types]                                            
cells[, cell_type := ifelse(!is.na(min_cell_type), cell_types[[min_cell_type]], "UNASSIGNED"), by = CellName]                                                                 
cells[, c(cell_types, "min_cell_type") := NULL] 

############# Prioritize Populations #############
fwrite(cells, out_file_name)
