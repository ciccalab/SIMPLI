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
plot_me <- arguments[[3]]
comparison_column <- arguments[[4]]   
cell_population_metadata_file_name <- arguments[[5]]
output_folder <- arguments[[6]]   
cell_mask_file_list <- unlist(arguments[7:length(arguments)])

######## Cell data #######
Cells <- fread(cell_file_name)
Samples <- fread(sample_metadata_file_name)
Samples[[comparison_column]] <- ifelse(is.na(Samples[[comparison_column]]), "NA", Samples[[comparison_column]])
suppressWarnings(Cells[, color := NULL])
Cells <- merge(suppressWarnings(Cells[, -comparison_column, with = F]), Samples,
	by.x = "Metadata_sample_name", by.y = "sample_name")
sample_names <- unique(Cells$Metadata_sample_name)


###################### Read cell type colors from the cell type metadata  ##########################
population_color <- fread(cell_population_metadata_file_name)
color_list <- population_color$color
color_list[length(color_list) + 1] <- "#888888"
names(color_list) <- c(population_color$cell_type, "UNASSIGNED")
rm(population_color)

###################### Barplots ##########################
if (plot_me){
	sample_barplot <- Barplotter(Cells, "Metadata_sample_name", "Metadata_sample_name", "cell_type",
		color_list, "Cell type cells / total cells %")
}
category_barplot <- Barplotter(Cells[Cells[[comparison_column]] != "NA"], "Metadata_sample_name",
	comparison_column, "cell_type",	color_list, "Cell type cells / total cells %")

###################### Boxplots by Population ##########################
n_categories <- sum(unique(Samples[[comparison_column]]) != "NA")
sample_colors <- Samples[, color]
names(sample_colors) <- Samples[, color]

# Boxplots should be made only when we have 2 groups to compare
if(n_categories == 2){
	plot_dataset <- copy(Cells)		
	plot_dataset <- plot_dataset[plot_dataset[[comparison_column]] != "NA"]
	plot_dataset[, n_cells := .N, by = c("Metadata_sample_name", "cell_type")]
	plot_dataset[, total := .N, by = "Metadata_sample_name"]
	plot_dataset[, percentage := n_cells / total * 100]
	population_boxplots <- list_boxplotter(plot_dataset, "Metadata_sample_name", "cell_type", comparison_column, "percentage",
		"Cell type cells / total cells %", "color", sample_colors)
	population_boxplots <- lapply(population_boxplots, function(x){x$Plot})
}

###################### Cell Overlays ##########################
if(plot_me){
	cell_mask_file_names <- sapply(sample_names, function(sample_name){
	cell_mask_file_name <- grep(sample_name, cell_mask_file_list, value = T)
	})
	names(cell_mask_file_names) <- sample_names

	cell_overlays <- lapply(sample_names, function(sample){
		cell_mask <- load_image(cell_mask_file_names[[sample]])
		cell_list <- lapply(names(color_list), function(type){
			Cells[Metadata_sample_name == sample & cell_type == type, ObjectNumber]
		})
		cell_overlayer(0, cell_mask, cell_list, color_list)	
	})
	names(cell_overlays) <- sample_names
}

################ Output ######################

# Barlots
barplot_output_folder <- file.path(output_folder, "Barplots")  
dir.create(barplot_output_folder, recursive = T, showWarnings = F)
pdf_plotter(category_barplot, filename =  paste0(barplot_output_folder,	"/barplot_by_", comparison_column, ".pdf"))
if(plot_me){
	# Barlots
	pdf_plotter(sample_barplot, filename =  paste0(barplot_output_folder,	"/barplot_by_sample.pdf"))

	# Cell Overlays
	overlay_output_folder <- file.path(output_folder, "Overlays")  
	dir.create(overlay_output_folder, recursive = T, showWarnings = F)
	lapply(sample_names, function(sample){
		writeImage(cell_overlays[[sample]], paste0(overlay_output_folder, "/overlay-", sample, ".tiff"),
			bits.per.sample = 8, compression = "LZW")
	})
}

# Boxplots by population
if(n_categories == 2){
	boxplot_output_folder <- file.path(output_folder, "Boxplots", comparison_column)  
	dir.create(boxplot_output_folder, recursive = T, showWarnings = F)
	lapply(names(color_list), function(population){
		pdf_plotter(population_boxplots[[population]], filename =  paste0(boxplot_output_folder,
			"/boxplot_", population, ".pdf"))
	})
}
