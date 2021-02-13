####### SETUP #######
rm(list = ls())
library(data.table)
library(uwot)
library(EBImage)
source("/opt/Plot_Functions.R")

arguments <- commandArgs(trailingOnly = TRUE)                                                                              
cell_file_name <- arguments[[1]]                                                                                             
sample_metadata_file_name <- arguments[[2]]
cell_cell_type_metadata_file_name <- arguments[[3]]
output_folder <- arguments[[4]]   
cell_mask_file_list <- unlist(arguments[5:length(arguments)])

######## Cell data #######
Cells <- fread(cell_file_name)
Samples <- fread(sample_metadata_file_name)
Samples[is.na(comparison), comparison := "NA" ]
suppressWarnings(Cells[, color := NULL])
suppressWarnings(Cells[, comparison := NULL])
Cells <- merge(Cells, Samples, by.x = "Metadata_sample_name", by.y = "sample_name")
sample_names <- Samples[, sample_name]

Comparison_Cells <- Cells[comparison != "NA",]
comparison_sample_names <- Samples[comparison != "NA", sample_name]

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
if(length(sample_names) > 0){
	unassigned_sample_barplot <- Barplotter(Cells, "Metadata_sample_name", "Metadata_sample_name", "cell_type",
		color_list, "Cell type cells / total cells %")
	assigned_sample_barplot <- Barplotter(Cells[cell_type != "UNASSIGNED"], "Metadata_sample_name", "Metadata_sample_name",
		"cell_type", color_list, "Cell type cells / total cells %")
	if(length(comparison_sample_names) > 1){
		unassigned_comparison_barplot <- Barplotter(Comparison_Cells, "Metadata_sample_name", "comparison", "cell_type", color_list,
			"Cell type cells / total cells %")
		assigned_comparison_barplot <- Barplotter(Comparison_Cells[cell_type != "UNASSIGNED"], "Metadata_sample_name", "comparison",
			"cell_type", color_list, "Cell type cells / total cells %")
	}
}

###################### Boxplots by Population ##########################
n_categories <- length(unique(Samples[comparison != "NA", comparison]))
sample_colors <- Samples[, color]
names(sample_colors) <- Samples[, color]
if(n_categories == 2){
# Boxplots should be made only when we have 2 groups to compare
	unassigned_plot_data <- copy(Comparison_Cells)		
	unassigned_plot_data[, n_cells := .N, by = c("Metadata_sample_name", "cell_type")]
	unassigned_plot_data[, total := .N, by = "Metadata_sample_name"]
	unassigned_plot_data[, percentage := n_cells / total * 100]
	unassigned_cell_type_boxplots <- list_boxplotter(unassigned_plot_data, "Metadata_sample_name", "cell_type", "comparison", "percentage",
		"Cell type cells / total cells %", "color", sample_colors)
	unassigned_cell_type_boxplots <- lapply(unassigned_cell_type_boxplots, function(x){x$Plot})
	assigned_plot_data <- copy(Comparison_Cells)[cell_type != "UNASSIGNED", ]
	assigned_plot_data[, n_cells := .N, by = c("Metadata_sample_name", "cell_type")]
	assigned_plot_data[, total := .N, by = "Metadata_sample_name"]
	assigned_plot_data[, percentage := n_cells / total * 100]
	assigned_cell_type_boxplots <- list_boxplotter(assigned_plot_data, "Metadata_sample_name", "cell_type", "comparison", "percentage",
		"Cell type cells / total cells %", "color", sample_colors)
	assigned_cell_type_boxplots <- lapply(assigned_cell_type_boxplots, function(x){x$Plot})
}

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
if(length(sample_names) > 1){
	multi_pdf_plotter(list(unassigned_sample_barplot, unassigned_comparison_barplot), filename = paste0(barplot_output_folder, "/barplots.pdf"),
		n_col = 1, n_row = 2)
	multi_pdf_plotter(list(assigned_sample_barplot, assigned_comparison_barplot), filename = paste0(barplot_output_folder, "/assigned_only_barplots.pdf"),
		n_col = 1, n_row = 2)
} else if (length(sample_names > 0)){
	pdf_plotter(filename = paste0(barplot_output_folder, "/barplots.pdf"), plot = unassigned_sample_barplot)
	pdf_plotter(filename = paste0(barplot_output_folder, "/assigned_ony_barplots.pdf"), plot = assigned_sample_barplot)
}
# Boxplots by comparison
if(n_categories == 2){	
	boxplot_output_folder <- file.path(output_folder, "Boxplots")  
	dir.create(boxplot_output_folder, recursive = T, showWarnings = F)
	multi_pdf_plotter(unassigned_cell_type_boxplots, filename = paste0(boxplot_output_folder, "/boxplots.pdf"))
	multi_pdf_plotter(assigned_cell_type_boxplots, filename = paste0(boxplot_output_folder, "/assigned_only_boxplots.pdf"))
}
