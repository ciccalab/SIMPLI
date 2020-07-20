####### SETUP #######
rm(list = ls())
library(data.table)
source("/opt/Plot_Functions.R")

arguments <- commandArgs(trailingOnly = TRUE)                                                                              
print(arguments)                                                                                                           
area_file_name <- arguments[[1]]                                                                                             
sample_metadata_file_name <- arguments[[2]]
comparison_columns <- strsplit(arguments[[3]], ",")[[1]]   
output_folder <- arguments[[4]]   

######## Area Measurements Data #######
Areas <- fread(area_file_name)
Samples <- fread(sample_metadata_file_name)
Samples[, c(comparison_columns) :=  lapply(.SD, function(x){ifelse(is.na(x), "NA", x)}), .SDcols = comparison_columns]
suppressWarnings(Areas[, color := NULL])
suppressWarnings(Areas[, comparison_columns := NULL])
Areas <- merge(Areas, Samples, by = "sample_name")

Areas[, marker_combination := paste0(marker, "_ON_", main_marker)]
marker_combinations <- unique(Areas$marker_combination)

for(comparison_column in comparison_columns){
	###################### Boxplots by Marker ##########################
	n_categories <- sum(unique(Samples[[comparison_column]]) != "NA")
	# Boxplots should be made only when we have 2 groups to compare
	if(n_categories == 2){
		sample_colors <- Samples[, color]
		names(sample_colors) <- Samples[, color]
		area_boxplots <- list_boxplotter(Areas[Areas[[comparison_column]] != "NA"], "sample_name", "marker_combination",
			comparison_column, "percentage", "Marker area %", "color", sample_colors)
		area_boxplots <- lapply(area_boxplots, function(x){x$Plot})
	###################### PDF output ##########################
		boxplot_output_folder <- file.path(output_folder, "Boxplots")  
		dir.create(boxplot_output_folder, recursive = T, showWarnings = F)
		multi_pdf_plotter(plots = area_boxplots, filename = file.path(boxplot_output_folder,
			paste0("all_boxplots-", comparison_column, ".pdf")))
	}
}	
