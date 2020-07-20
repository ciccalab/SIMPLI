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
comparison_columns <- strsplit(arguments[[3]], ",")[[1]]    
cell_cell_type_metadata_file_name <- arguments[[4]]
output_folder <- arguments[[5]]   
cell_mask_file_list <- unlist(arguments[6:length(arguments)])

######## Cell data #######
Cells <- fread(cell_file_name)
Samples <- fread(sample_metadata_file_name)
print(comparison_columns)
Samples[, c(comparison_columns) := lapply(.SD, function(x){ifelse(is.na(x), "NA", x)}), .SDcols = comparison_columns]
suppressWarnings(Cells[, color := NULL])
suppressWarnings(Cells[, comparison_columns := NULL])
Cells <- merge(Cells, Samples, by.x = "Metadata_sample_name", by.y = "sample_name")

sample_names <- Samples[, sample_name]

###################### Read cell type colors from the cell type metadata  ##########################
cell_type_color <- fread(cell_cell_type_metadata_file_name)
color_list <- cell_type_color$color
color_list[length(color_list) + 1] <- "#888888"
names(color_list) <- c(cell_type_color$cell_type, "UNASSIGNED")
rm(cell_type_color)

###################### Cell Overlays ##########################
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

legend_table <- data.table(cell_type = names(color_list), color = color_list)
legend_table[, X := (.I-1)%%round(sqrt(.N))*round(sqrt(.N)) + 1]
legend_table[, text_X := X + round(sqrt(.N)) - 0.5]
legend_table[, Y := -(0:round(sqrt(.N))), by = X]
legend_plot <- ggplot(data = legend_table) +
	geom_tile(mapping = aes(x = X, y = Y, fill = cell_type), width = 0.5, height = 0.5, color = "black") +
	geom_text(mapping = aes(x = text_X, y = Y, label = cell_type), hjust = "right", size = 3) +
	scale_fill_manual(values = color_list, guide = "none") +
	coord_equal() +
	theme(axis.line=element_blank(), axis.text = element_blank(), axis.title = element_blank(),
		axis.ticks = element_blank(), panel.background = element_blank(), panel.border = element_blank(),
		panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())


###################### Barplots ##########################
sample_barplot <- Barplotter(Cells, "Metadata_sample_name", "Metadata_sample_name", "cell_type",
	color_list, "Cell type cells / total cells %")

comparison_barplots <- NULL
if(length(sample_names) > 1){
	comparison_barplots <- lapply(comparison_columns, function(comparison_column){
		Barplotter(Cells[Cells[[comparison_column]] != "NA"], "Metadata_sample_name",
			comparison_column, "cell_type",	color_list, "Cell type cells / total cells %")
	})
	names(comparison_barplots) <- comparison_columns
}

###################### Boxplots by Population ##########################
comparison_boxplots <- lapply(comparison_columns, function(comparison_column){
	n_categories <- sum(unique(Samples[[comparison_column]]) != "NA")
	sample_colors <- Samples[, color]
	names(sample_colors) <- Samples[, color]
	cell_type_boxplots <- list()
	if(n_categories == 2){
	# Boxplots should be made only when we have 2 groups to compare
		plot_dataset <- copy(Cells)		
		plot_dataset <- plot_dataset[plot_dataset[[comparison_column]] != "NA"]
		plot_dataset[, n_cells := .N, by = c("Metadata_sample_name", "cell_type")]
		plot_dataset[, total := .N, by = "Metadata_sample_name"]
		plot_dataset[, percentage := n_cells / total * 100]
		cell_type_boxplots <- list_boxplotter(plot_dataset, "Metadata_sample_name", "cell_type", comparison_column, "percentage",
			"Cell type cells / total cells %", "color", sample_colors)
		lapply(cell_type_boxplots, function(x){x$Plot})
	}
})
names(comparison_boxplots) <- comparison_columns

################ Output ######################
# Cell Overlays
overlay_output_folder <- file.path(output_folder, "Overlays")  
dir.create(overlay_output_folder, recursive = T, showWarnings = F)
lapply(sample_names, function(sample){
	writeImage(cell_overlays[[sample]], paste0(overlay_output_folder, "/overlay-", sample, ".tiff"),
		bits.per.sample = 8, compression = "LZW")
})
pdf_plotter(file.path(overlay_output_folder, "overlay_legend.pdf"), legend_plot)

# Barlots
barplot_output_folder <- file.path(output_folder, "Barplots")  
dir.create(barplot_output_folder, recursive = T, showWarnings = F)
multi_pdf_plotter(c(list(sample_barplot), comparison_barplots), filename = paste0(barplot_output_folder, "/barplots.pdf"), n_col = 2, n_row = 2)
# Boxplots by comparison
if(any(lengths(comparison_boxplots) > 0)){	
	boxplot_output_folder <- file.path(output_folder, "Boxplots")  
	dir.create(boxplot_output_folder, recursive = T, showWarnings = F)
	for(comparison_column in names(comparison_boxplots[lengths(comparison_boxplots) > 0]))
	{
		multi_pdf_plotter(comparison_boxplots[[comparison_column]],
			filename = paste0(boxplot_output_folder, "/boxplots-", comparison_column, ".pdf"))
	}
}

