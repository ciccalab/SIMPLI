####### SETUP #######
rm(list = ls())
library(data.table)
library(EBImage)

args <- commandArgs(trailingOnly = TRUE)
my_sample_name <- args[[1]]
raw_metadata_file_name <- args[[2]]
file_type<- args[[3]]
output_path <- args[[4]]
output_metadata <- args[[5]]
cellprofiler_metadata <- args[[6]]

####### Helper Functions #######
tiff_loader <- function(file_name, all = FALSE)
{
  readImage(file_name, type = "tiff", all, as.is = TRUE)
}

normalizer <- function(picture_data)
{
  picture_data <- picture_data / sort(picture_data)[0.99*length(picture_data)]
  picture_data[picture_data > 1] <- 1
  return(picture_data)
}

image_writer <- function(picture, sample_name, marker, output_dir, suffix = "normalized.tiff")
{
  file_name <- file.path(output_dir, paste(sample_name, marker, suffix, sep = "-"))
  writeImage(picture, file_name, type = "tiff", quality = 100,
    bits.per.sample = 16L, compression = "none", reduce = TRUE)
  file_name <- normalizePath(file_name)
  return(file_name)
}

process_image_single <- function(image, sample_name, marker, metal, output_dir)
{
  my_image <- Image(data = image@.Data)   
  my_image@.Data <- normalizer(my_image@.Data)
  output_file_name <- image_writer(my_image, sample_name, marker, output_dir)
  return(c(sample_name = sample_name, metal = metal, label = marker,
    normalized_file = output_file_name))
}

####### Load the metatdata #######
# Expects a .csv file with these columns: sample_name,roi_name,metal,label,raw_tiff_file_name        
raw_metadata <- fread(raw_metadata_file_name, header = TRUE, sep = ",")
raw_metadata <- raw_metadata[sample_name == my_sample_name, ]

####### Image processing #######
if (file_type == "ome") {
  ome_image <- tiff_loader(raw_metadata$raw_tiff_file_name, all = T)
  ome_image@.Data <- sapply(seq(1, dim(ome_image@.Data)[[3]]), function(n){
    normalizer(ome_image@.Data[,,n])}, simplify = "array")
  output_file_name <- image_writer(ome_image, my_sample_name, "ALL", output_path,
    suffix = "normalized.ome.tiff")
  metadata <- copy(raw_metadata)
  metadata[, normalized_file := output_file_name]
  metadata[, raw_tiff_file_name := NULL]
} else if (file_type == "single") {
  metadata <- lapply(seq(1, raw_metadata[, .N]), function(n){
    meta <- raw_metadata[n]
    my_image <- tiff_loader(meta$raw_tiff_file_name)
    process_image_single(my_image, meta$sample_name,  meta$label, meta$metal,
     output_path)  
  })
  metadata <- rbindlist(lapply(metadata, function(x){as.data.table(t(x))}))
} else {
  print(paste0("file type: " , file_type, "not 'single' or 'ome'"))
  quit(status = 99)
}
fwrite(metadata, file = output_metadata, sep = ",")

file_columns <- grep("file", colnames(metadata), value = T)
cp_metatada <- dcast(metadata, sample_name ~ label, value.var = file_columns)
setnames(cp_metatada, unique(metadata$label), paste0("Image_FileName_", unique(metadata$label)))

metadata_columns <- grep("file", colnames(metadata), invert = T, value = T)
metadata_columns <- metadata_columns[!(metadata_columns %in% c("sample_name", "label", "metal"))]
if (length(metadata_columns) > 0) {
  metadata <- unique(metadata[, c("sample_name", metadata_columns), with = F])
  cp_metatada <- merge(cp_metatada, metadata, by = "sample_name")
  setnames(cp_metatada, metadata_columns, paste0("Metadata_", metadata_columns))
  
}
setnames(cp_metatada, "sample_name", "Metadata_sample_name")
fwrite(cp_metatada, cellprofiler_metadata)