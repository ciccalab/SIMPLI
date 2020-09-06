####### SETUP #######
rm(list = ls())
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
unannotated_cells_file_name <- args[[1]]
threshold_metadata_file_name <- args[[2]]
out_file_name <- args[[3]]

############ Load the unannotated cells ##########
cells <- fread(unannotated_cells_file_name)

############################# Marker combinations ###########################
# Expects a .csv file with 3 columns:
# - cell_type
# - threshold_marker
# - threshold_value
threshold_metadata <- fread(threshold_metadata_file_name)

############################# Cell Selection ###########################
cell_types <- threshold_metadata[, cell_type]
bool_columns <- paste0(cell_types, "_bool")
expression_columns <- paste0(cell_types, "_expression")
names(expression_columns) <- cell_types
names(bool_columns) <- cell_types
cells[, (bool_columns) := lapply(threshold_metadata[, cell_type], function(type){
  get(threshold_metadata[cell_type == type, threshold_marker]) > threshold_metadata[cell_type == type, threshold_value]
})]
cells[, (expression_columns) := lapply(threshold_metadata[, cell_type], function(type){
  get(threshold_metadata[cell_type == type, threshold_marker])})]

cells[, (expression_columns) := lapply(cell_types, function(type){
  ifelse(get(bool_columns[[type]]), get(expression_columns[[type]]), NA_real_)
})]

cells[, cell_type := apply(.SD, 1, which.max), .SDcols = expression_columns]
cells[, cell_type := cell_types[as.numeric(cell_type)]]
cells[is.na(cell_type), cell_type := "UNASSIGNED"]
cells[, (c(expression_columns, bool_columns)) := NULL]

########## Output ###################
cells[, CellName := paste0(Metadata_sample_name, "-", ObjectNumber)]
fwrite(cells, out_file_name)
