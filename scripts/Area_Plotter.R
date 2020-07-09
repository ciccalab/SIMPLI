####### SETUP #######
rm(list = ls())
library(data.table)
source("/opt/Plot_Functions.R")

arguments <- commandArgs(trailingOnly = TRUE)                                                                              
print(arguments)                                                                                                           
area_file_name <- arguments[[1]]                                                                                             
sample_metadata_file_name <- arguments[[2]]
category_column <- arguments[[3]]   
output_folder <- arguments[[4]]   

######## Area Measurements Data #######
Areas <- fread(area_file_name)
Samples <- fread(sample_metadata_file_name)
Areas <- merge(suppressWarnings(Areas[, -category_column, with = F]), Samples,
	by = "sample_name")

Areas[, marker_combination := paste0(marker, "_ON_", main_marker)]
marker_combinations <- unique(Areas$marker_combination)

###################### Boxplots by Marker ##########################
n_categories <- length(unique(Samples[[category_column]]))
# Boxplots should be made only when we have 2 groups to compare
if(n_categories == 2){
		area_boxplots <- list_boxplotter(Areas, "sample_name", "marker_combination", "category", "percentage", "Marker area / normaliser area %")
		area_boxplots <- lapply(area_boxplots, function(x){x$Plot})
}

###################### PDF output ##########################
if(n_categories == 2){
	output_folder <- file.path(output_folder)
	lapply(marker_combinations, function(mc){
		main_marker <- unique(Areas[marker_combination == mc, main_marker])
		boxplot_output_folder <- file.path(output_folder, "Boxplots", main_marker)  
		dir.create(boxplot_output_folder, recursive = T, showWarnings = F)
		pdf_plotter(area_boxplots[[mc]], filename =  paste0(boxplot_output_folder, "/boxplot_", mc, ".pdf"))
	})
}	
