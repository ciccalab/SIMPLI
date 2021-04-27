####### SETUP #######
rm(list = ls())
library(data.table)

arguments <- commandArgs(trailingOnly = TRUE)
file_type <- arguments[[1]]
outputh_path <- arguments[[2]]

filenames <- as.character(arguments[3 : length(arguments)])

cast_columns <- c("URL", "Frame")

if(file_type != "ome" & file_type != "single"){
  print(paste0("file type: " , file_type, "not 'single' or 'ome'"))
  quit(status = 99)
}

##### Prepare CellProfiler4 compatible wide format metadata
metadata_maker <- function(filename)
{
  metadata <- fread(filename)
  if(file_type == "ome" && length(filenames) == 1) metadata[, Frame := .I - 1]
  if(file_type == "single" || length(filenames) > 1) metadata[, Frame := 0]
  metadata[, URL := paste0("file:///", file_name)]
  metadata[, file_name := NULL]
  cp_metadata <- dcast(metadata, sample_name ~ label, value.var = cast_columns)
  metadata_columns <- colnames(metadata)[!(colnames(metadata) %in% cast_columns)]
  metadata_columns <- metadata_columns[!(metadata_columns %in% c("sample_name", "label", "marker"))]
  if (length(metadata_columns) > 0) {
    metadata <- unique(metadata[, c("sample_name", metadata_columns), with = F])
    cp_metadata <- merge(cp_metadata, metadata, by = "sample_name")
    setnames(cp_metadata, metadata_columns, paste0("Metadata_", metadata_columns))
  }
  setnames(cp_metadata, "sample_name", "Metadata_sample_name")
}

metadata_list <- lapply(filenames, metadata_maker)

filenames <- basename(filenames)
names(metadata_list) <- filenames

##### Output metadata in CellProfiler4 compatible wide format
lapply(names(metadata_list), function(filename){
  out_filename <- paste0(sub("-preprocessed_metadata.csv", "", filename),
    "-cp4-preprocessed_metadata.csv")
  out_filename <- file.path(outputh_path, out_filename)
  fwrite(metadata_list[[filename]], file = out_filename, sep = ",")  
})
