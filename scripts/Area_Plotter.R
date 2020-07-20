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
main_markers <- unique(Areas[, main_marker]) 

for(comparison_column in comparison_columns){
	###################### Boxplots by Marker ##########################
	n_categories <- sum(unique(Samples[[comparison_column]]) != "NA")
	# Boxplots should be made only when we have 2 groups to compare
	if(n_categories == 2){
		sample_colors <- Samples[, color]
		names(sample_colors) <- Samples[, color]
		area_boxplots <- lapply(main_markers, function(mm){
			boxplots <- list_boxplotter(Areas[Areas[[comparison_column]] != "NA" & main_marker == mm], "sample_name",
				"marker_combination", comparison_column, "percentage", "Marker area %", "color", sample_colors)
			boxplots <- lapply(boxplots, function(x){x$Plot})
		})
		names(area_boxplots) <- main_markers
	###################### PDF output ##########################
		for (mm in main_markers){
			boxplot_output_folder <- file.path(output_folder, "Boxplots", mm)  
			dir.create(boxplot_output_folder, recursive = T, showWarnings = F)
			multi_pdf_plotter(plots = area_boxplots[[mm]], filename = file.path(boxplot_output_folder,
				paste0(mm, "_boxplots-", comparison_column, ".pdf")))
		}
	}
}	
