####### SETUP #######
rm(list = ls())
library(data.table)
source("/opt/Plot_Functions.R")

arguments <- commandArgs(trailingOnly = TRUE)                                                                              
area_file_name <- arguments[[1]]                                                                                             
sample_metadata_file_name <- arguments[[2]]
output_folder <- arguments[[3]]   

######## Area Measurements Data #######
Areas <- fread(area_file_name)
Samples <- fread(sample_metadata_file_name)
Samples[is.na(comparison), comparison := "NA" ]
suppressWarnings(Areas[, color := NULL])
suppressWarnings(Areas[, comparisons := NULL])
Areas <- merge(Areas, Samples, by = "sample_name")

Areas[, marker_combination := paste0(marker, "_ON_", main_marker)]
marker_combinations <- unique(Areas$marker_combination)
main_markers <- unique(Areas[, main_marker]) 

###################### Boxplots by Marker ##########################
n_categories <- length(unique(Samples[comparison != "NA", comparison]))
# Boxplots should be made only when we have 2 groups to compare
if(n_categories == 2){
	sample_colors <- Samples[, color]
	names(sample_colors) <- Samples[, color]
	area_boxplots <- lapply(main_markers, function(mm){
		boxplots <- list_boxplotter(Areas[main_marker == mm], "sample_name",
			"marker_combination", "comparison", "percentage", "Marker area %", "color", sample_colors)
		boxplots <- lapply(boxplots, function(x){x$Plot})
	})
	names(area_boxplots) <- main_markers
###################### PDF output ##########################
	for (mm in main_markers){
		boxplot_output_folder <- file.path(output_folder, "Boxplots", mm)  
		dir.create(boxplot_output_folder, recursive = T, showWarnings = F)
		multi_pdf_plotter(plots = area_boxplots[[mm]], filename = file.path(boxplot_output_folder,
			paste0(mm, "_boxplots.pdf")))
	}
}
