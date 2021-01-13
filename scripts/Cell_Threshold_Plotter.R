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
threshold_metadata_file_name <- arguments[[3]]
high_color <- arguments[[4]]                                                                                               
mid_color <- arguments[[5]]                                                                                                
low_color <- arguments[[6]]                                                                                                
output_folder <- arguments[[7]]   
cell_mask_file_list <- unlist(arguments[8:length(arguments)])

### Load the threshold metadata
threshold_metadata <- fread(threshold_metadata_file_name)
populations <- unique(threshold_metadata$cell_type)

Population_Colors <- lapply(populations, function(pop){
  color_list <- threshold_metadata[cell_type == pop, color]
  color_list[length(color_list) + 1] <- "#888888"
  names(color_list) <- c(threshold_metadata[cell_type == pop,  population_name], "UNASSIGNED")
  color_list
})
names(Population_Colors) <- populations

Markers <- lapply(populations, function(pop){
  unique(unlist(sapply(threshold_metadata[cell_type == pop, plotting_markers], strsplit, "@")))
})
names(Markers) <- populations

######## Cell data #######
Cells <- fread(cell_file_name)
Samples <- fread(sample_metadata_file_name)
Samples[is.na(comparison), comparison := "NA" ]
suppressWarnings(Cells[, color := NULL])
suppressWarnings(Cells[, comparison := NULL])
Cells <- merge(Cells, Samples, by.x = "Metadata_sample_name", by.y = "sample_name")
Cells <- Cells[comparison != "NA",]

sample_names <- Samples[comparison != "NA", sample_name]

#################################### Cell overlays
cell_mask_file_names <- sapply(sample_names, function(sample_name){
  cell_mask_file_name <- grep(sample_name, cell_mask_file_list, value = T)
})
names(cell_mask_file_names) <- sample_names

cell_overlays <- lapply(populations, function(pop){
  pics <- lapply(sample_names, function(sample){
    cell_mask <- load_image(cell_mask_file_names[[sample]])
    cell_list <- lapply(names(Population_Colors[[pop]]), function(type){
      Cells[Metadata_sample_name == sample & eval(parse(text = paste0(pop, "_Thresholded"))) == type, ObjectNumber]
    })
    cell_overlayer(0, cell_mask, cell_list, Population_Colors[[pop]])	
  })
  names(pics) <- sample_names
  return(pics)
})
names(cell_overlays) <- populations

legend_tables <- lapply(populations, function(pop){
  legend_table <- data.table(cell_type = names(Population_Colors[[pop]]), color = Population_Colors[[pop]])
  legend_table[, X := (.I-1)%%round(sqrt(.N))*round(sqrt(.N)) + 1]
  legend_table[, text_X := X + round(sqrt(.N)) - 0.5]
  legend_table[, Y := -(0:round(sqrt(.N))), by = X]
  legend_plot <- ggplot(data = legend_table) +
    geom_tile(mapping = aes(x = X, y = Y, fill = cell_type), width = 0.5, height = 0.5, color = "black") +
    geom_text(mapping = aes(x = text_X, y = Y, label = cell_type), hjust = "right", size = 3) +
    scale_fill_manual(values = Population_Colors[[pop]], guide = "none") +
    coord_equal() +
    theme(axis.line=element_blank(), axis.text = element_blank(), axis.title = element_blank(),
          axis.ticks = element_blank(), panel.background = element_blank(), panel.border = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.background = element_blank())
})
names(legend_tables) <- populations

#################################### Heatmaps
heatmaps <- lapply(populations, function(pop){
  heatmapper(Cells[!is.na(eval(parse(text = paste0(pop, "_Thresholded")))), c(paste0(pop, "_Thresholded"),
  Markers[[pop]]), with = F], res_column = paste0(pop, "_Thresholded"), cols = Markers[[pop]],
  high_color, mid_color, low_color)
})
names(heatmaps) <- populations

#################################### Boxplots
n_categories <- length(unique(Samples[comparison != "NA", comparison]))

sample_colors <- Samples[, color]
names(sample_colors) <- Samples[, color]

if(n_categories == 2){
  comparison_boxplots <- lapply(populations, function(pop){
    plot_dataset <- copy(Cells)[!is.na(eval(parse(text = paste0(pop, "_Thresholded"))))]		
    plot_dataset[, total := .N, by = "Metadata_sample_name"]
    plot_dataset[, n_cells := .N, by = c("Metadata_sample_name", paste0(pop, "_Thresholded"))]
    plot_dataset[, frac := n_cells / total * 100]
    cluster_boxplots <- list_boxplotter(plot_dataset, "Metadata_sample_name",
      paste0(pop, "_Thresholded"), "comparison", "frac",
      "Cells in subpopulation / total cells %", "color", sample_colors)
    lapply(cluster_boxplots, function(x){x$Plot})
  })
  names(comparison_boxplots) <- populations
}

#################################### Barplots
if(length(sample_names) > 0){
  sample_barplots <- lapply(populations, function(pop){
    plot_dataset <- copy(Cells)[!is.na(eval(parse(text = paste0(pop, "_Thresholded"))))]		
    Barplotter(plot_dataset, "Metadata_sample_name", "Metadata_sample_name",
      paste0(pop, "_Thresholded"), Population_Colors[[pop]],
      paste0("Subpopulation type cells / ", pop, " total %"))
  })
  names(sample_barplots) <- populations
  if(length(sample_names) > 1){
    comparison_barplots <- lapply(populations, function(pop){
      plot_dataset <- copy(Cells)[!is.na(eval(parse(text = paste0(pop, "_Thresholded"))))]		
      Barplotter(plot_dataset, "Metadata_sample_name", "comparison", paste0(pop, "_Thresholded"),
        Population_Colors[[pop]],  paste0("Subpopulation type cells / ", pop, " total %"))
    })
    names(comparison_barplots) <- populations
  }
}

################## Density Plots
density_plots <- lapply(populations, function(pop){
  subpopulations <- names(Population_Colors[[pop]])
  subpopulations <- subpopulations[subpopulations != "UNASSIGNED"]
  subpopulation_plots <- lapply(subpopulations, function(sub_pop){
    plot_dataset <- copy(Cells)[!is.na(eval(parse(text = paste0(pop, "_Thresholded"))))]
    threshold_markers <- all.vars(parse(text = threshold_metadata[cell_type == pop &
      population_name == sub_pop, threshold_expression]))
    marker_densities <- lapply(threshold_markers, function(marker){
      histogram_density(plot_dataset, marker, paste(pop, sub_pop, marker, sep = " "),
        marker, paste0(pop, " cells"), geom = "density", color = "black",
        fill = Population_Colors[[pop]][[sub_pop]], alpha = 0.3, log_x = T, log_y = T)
    })
    names(marker_densities) <- threshold_markers
    return(marker_densities)
  })
  names(subpopulation_plots) <- subpopulations
  return(subpopulation_plots)
})
names(density_plots) <- populations

##################### OUTPUT ########################
# Overlays
for(pop in populations){
  overlay_output_folder <- file.path(output_folder, "Overlays", pop) 
  dir.create(overlay_output_folder, recursive = T, showWarnings = F)
  for(s_name in sample_names){
     writeImage(cell_overlays[[pop]][[s_name]], paste0(overlay_output_folder, "/", pop, "-overlay-",
       s_name, ".tiff"), bits.per.sample = 8, compression = "LZW")
  }
  pdf_plotter(file.path(overlay_output_folder, "overlay_legend.pdf"), legend_tables[[pop]])
}

# Heatmaps
heatmap_output_folder <- file.path(output_folder, "Heatmaps") 
dir.create(heatmap_output_folder, recursive = T, showWarnings = F)
for(pop in populations){
  pdf_plotter(paste0(heatmap_output_folder, "/", pop, "-heatmap.pdf"), heatmaps[[pop]])
}

# Barlots
barplot_output_folder <- file.path(output_folder, "Barplots")  
dir.create(barplot_output_folder, recursive = T, showWarnings = F)
if(length(sample_names) > 1){
  for (pop in populations){
    multi_pdf_plotter(list(sample_barplots[[pop]], comparison_barplots[[pop]]),
      filename = paste0(barplot_output_folder, "/", pop, "-barplots.pdf"), n_col = 1, n_row = 2)
  }
} else if (length(sample_names > 0)){
  for (pop in populations){
    pdf_plotter(filename = paste0(barplot_output_folder, "/", pop, "-barplots.pdf"),
      plot = sample_barplots[[pop]])
  }
}

# Boxplots by comparison
if(n_categories == 2){
  boxplot_output_folder <- file.path(output_folder, "Boxplots")  
  dir.create(boxplot_output_folder, recursive = T, showWarnings = F)
  for (pop in populations){
    multi_pdf_plotter(comparison_boxplots[[pop]], filename = paste0(boxplot_output_folder, "/",
      pop, "-boxplots.pdf"))
  }
}

# Density Plots
for(pop in populations){
  for(sub_pop in names(density_plots[[pop]])){
    density_folder <- file.path(output_folder, "Density_Plots")  
    dir.create(density_folder, recursive = T, showWarnings = F)
    for (marker in names(density_plots[[pop]][[sub_pop]])){
      pdf_plotter(plot = density_plots[[pop]][[sub_pop]][[marker]],
			filename = paste0(density_folder, "/", sub_pop, "-", marker, "-density.pdf"))
    }
  }
}

