####### SETUP #######
rm(list = ls())
library(data.table)
library(uwot)
library(EBImage)
source("/opt/Plot_Functions.R")

arguments <- commandArgs(trailingOnly = TRUE)                                                                              
print(arguments)                                                                                                           
cell_file_name <- arguments[[1]]                                                                                             
sample_metadata_file_name <- arguments[[2]]
category_column <- arguments[[3]]   
cell_population_metadata_file_name <- arguments[[4]]
output_folder <- arguments[[5]]   
cell_mask_file_list <- unlist(arguments[6:length(arguments)])

######## Cell data #######
Cells <- fread(cell_file_name)
Samples <- fread(sample_metadata_file_name)
Cells <- merge(suppressWarnings(Cells[, -category_column, with = F]), Samples,
	by.x = "Metadata_sample_name", by.y = "sample_name")
Samples <- unique(Cells$Metadata_sample_name)


###################### Read cell type colors from the cell type metadata  ##########################
population_color <- fread(cell_population_metadata_file_name)
color_list <- population_color$color
color_list[length(color_list) + 1] <- "#888888"
names(color_list) <- c(population_color$cell_type, "UNASSIGNED")
rm(population_color)

###################### Barplots ##########################
category_barplot <- Barplotter(Cells, "Metadata_sample_name", category_column, cell_type, color_list, "Cell type cells / total cells %")
sample_barplot <- Barplotter(Cells, "Metadata_sample_name", "Metadata_sample_name", cell_type, color_list, "Cell type cells / total cells %")

###################### Boxplots by Population ##########################
n_categories <- length(unique(Cells[[category_column]]))

if(n_categories == 2){
	Cells[, n_cells := .N, by = c("Metadata_sample_name", "cell_type")]
	Cells[, total := .N, by = "Metadata_sample_name"]
	Cells[, percentage := n_cells / total * 100]
	population_boxplots <- list_boxplotter(Cells, "Metadata_sample_name", "cell_type", "category", "percentage",
		"Cell type cells / total cells %")
	population_boxplots <- lapply(population_boxplots, function(x){x$Plot})
}

###################### Cell Overlays ##########################
cell_mask_file_names <- sapply(Samples, function(sample_name){
 cell_mask_file_name <- grep(sample_name, cell_mask_file_list, value = T)
})
names(cell_mask_file_names) <- Samples

cell_overlays <- lapply(Samples, function(sample){
	cell_mask <- load_image(cell_mask_file_names[[sample]])
	cell_list <- lapply(names(color_list), function(type){
		Cells[Metadata_sample_name == sample & cell_type == type, ObjectNumber]
	})
	cell_overlayer(0, cell_mask, cell_list, color_list)	
})
names(cell_overlays) <- Samples

################ Output ######################

# Barlots
barplot_output_folder <- file.path(output_folder, "Barplots")  
pdf_plotter(category_barplot, filename =  paste0(barplot_output_folder,	"/barplot_by_category.pdf"))
pdf_plotter(sample_barplot, filename =  paste0(barplot_output_folder,	"/barplot_by_sample.pdf"))

# Boxplots by population
if(n_categories == 2){
	boxplot_output_folder <- file.path(output_folder, "Boxplots")  
	dir.create(boxplot_output_folder, recursive = T, showWarnings = F)
	lapply(names(color_list), function(population){
		pdf_plotter(population_boxplots[[population]], filename =  paste0(boxplot_output_folder,
			"/boxplot_", population, ".pdf"))
	})
}

# Cell Overlays
overlay_output_folder <- file.path(output_folder, "Overlays")  
dir.create(overlay_output_folder, recursive = T, showWarnings = F)
lapply(Samples, function(sample){
	writeImage(cell_overlays[[sample]], paste0(overlay_output_folder, "/overlay-", sample, ".tiff"),
		bits.per.sample = 8, compression = "LZW")
})

