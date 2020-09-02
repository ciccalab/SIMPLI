####### SETUP #######
rm(list = ls())
library(dplyr)
library(tidyr)
library(data.table)
library(EBImage)

args <- commandArgs(trailingOnly = TRUE)
measurement_file_name <- args[[1]]
metadata_file_name <- args[[2]]
out_file_name <- args[[3]]

############ Load the image metadata ##########
# Expects a .csv file with 3 columns:
# - sample_name
# - label
# - preprocessed_file_name
image_metadata <- fread(metadata_file_name)
all_samples <- unique(image_metadata$sample_name)

############################# Marker combinations ###########################
# Expects a .csv file with 2 columns:
# - marker
# - main_marker
marker_combinations_table <- fread(measurement_file_name)

fix_names <- function(name_vect)
{
  name_vect <- gsub(" ", "", name_vect)
  name_vect <- gsub("&", "_", name_vect)
  name_vect <- gsub("|", "_or_", name_vect, fixed = T)
  name_vect <- gsub("!", "not_", name_vect)
}

fix_expression <- function(name_vect)
{
  name_vect <- gsub("not_", "! ", name_vect)
  name_vect <- gsub("_or_", " | ", name_vect)
  name_vect <- gsub("_", " & ", name_vect)
}

missing_main_markers <- unique(marker_combinations_table$main_marker[!(marker_combinations_table$main_marker
  %in% marker_combinations_table$marker)])

missing_main_markers <- data.table(marker = missing_main_markers, main_marker = missing_main_markers)
marker_combinations_table <- rbind(marker_combinations_table, missing_main_markers)

marker_combinations <- lapply(1:marker_combinations_table[, .N], function(n){
  c(marker_combinations_table[n, marker], marker_combinations_table[n, main_marker])
})

names(marker_combinations) <- marker_combinations_table[, marker] 

############### Helper Functions ##############
tiff_loader <- function(file_name, all = FALSE)
{
  readImage(file_name, type = "tiff", all, as.is = TRUE)
}

get_ROI_areas <- function(samples)
{
  ROI_areas <- sapply(samples, function(name){
    img <- tiff_loader(image_metadata[sample_name == name, file_name][[1]])
    dim(img)[[1]] * dim(img)[[2]]
  })
  data.table(sample_name = samples, total_ROI_area = ROI_areas)
}

process_image <- function(name_sample, expressions)
{
  to_parse <- sapply(expressions, `[[`, 1)
  to_divide <- sapply(expressions, `[[`, 2)
  parsed <- sapply(to_parse, function(expr){parse(text = expr)})
  divided <- as.character(unlist(sapply(sapply(to_divide, function(expr){parse(text = expr)}), all.vars)))
  image_names <- unique(c(as.character(unlist(sapply(parsed, all.vars))), divided))
  marker_file_names <- sapply(image_names, function(marker_name){
    image_metadata[sample_name == name_sample & label == marker_name, file_name]})
  imgs <- lapply(marker_file_names, tiff_loader)
  imgs <- lapply(imgs, function(x){x & x})
  names(imgs) <- image_names
  processed <- lapply(parsed, function(x){eval(x, envir = imgs)})
  areas <- lapply(processed, sum)
  names(areas) <- paste0(to_parse, "_ON_", to_divide)
  return(areas)
}

calculate_areas <- function(samples, expressions)
{
  areas <- lapply(samples, function(sample_name){process_image(sample_name, expressions)})
  names(areas) <- samples
  areas <- rbindlist(areas, idcol = "sample_name")
  old_cols <- colnames(areas)[colnames(areas) != "sample_name"]
  cols <- fix_names(old_cols)
  setnames(areas, old_cols, cols)
  return(areas)
}

melter <- function(annotated_table, ROI_area_table)
{
  tbl <- copy(annotated_table)
  cols <- grep("_ON_", colnames(tbl), value = T)
  tbl <- melt(tbl, measure.vars = list(cols), id.vars = "sample_name", variable.name = "marker",
              value.name = "area")
  tbl[, bad_marker := marker]
  tbl[, marker := sub("_ON_.+", "", bad_marker)]
  tbl[, main_marker := sub(".+_ON_", "", bad_marker)]
  tbl[, bad_marker := NULL]
  tbl <- unique(merge(tbl, tbl[, .(sample_name, main_marker = marker, main_marker_area = area)],
                      by = c("sample_name", "main_marker")))
  tbl <- merge(tbl, ROI_area_table, by = "sample_name")
  tbl[main_marker == marker, main_marker_area := total_ROI_area]
  tbl[main_marker == marker, main_marker := "total_ROI_area"]
  tbl[, percentage := area / main_marker_area * 100]
  tbl[main_marker_area == 0, percentage := 0]
  return(tbl)
}

################ Area Measurements ##########################
cat("Started getting ROI areas.\n")
ROI_areas <- get_ROI_areas(all_samples)
cat("Started processing the areas.\n")
raw_areas <- calculate_areas(all_samples, marker_combinations)
cat("Started formatinng the area table to long form.\n")
melted_areas <- melter(raw_areas, ROI_areas)

################ Tabular Output by Sample #################
cat("Started the table output.\n")
fwrite(melted_areas, out_file_name, na = "NA")
cat("Finished!\n")
