####### SETUP #######
rm(list = ls())
library(data.table)
library(ggplot2)
library(ggrepel)
library(EBImage)
library(uwot)

arguments <- commandArgs(trailingOnly = TRUE)                                                                              
print(arguments)                                                                                                           
cell_file_name <- arguments[[1]]                                                                                             
sample_metadata_file_name <- arguments[[2]]
category_column <- arguments[[3]]   
cell_population <- arguments[[4]]                                                                                          
markers <- strsplit(arguments[[5]], "@")[[1]]                                                                              
resolutions <- as.numeric(strsplit(arguments[[6]], "@")[[1]])                                                              
high_color <- arguments[[7]]                                                                                               
low_color <- arguments[[8]]                                                                                                
output_folder <- arguments[[9]]   

resolutions <- paste0("res_", resolutions, "_ids")

######## Cell data #######
Cells <- fread(cell_file_name)
Cells <- Cells[cell_type == cell_population, ]
Samples <- fread(sample_metadata_file_name)
Cells <- merge(suppressWarnings(Cells[, -category_column, with = F]), Samples,
	by.x = "Metadata_sample_name", by.y = "sample_name")

######## Heatmap plotter #######
cluster_heatmapper <- function(plot_dataset, res_column, cols){
  plot_data <- copy(plot_dataset)
  setnames(plot_data, res_column, "cluster")
  plot_data <- plot_data[, lapply(.SD, mean), by = cluster, .SDcols = cols]
  plot_data[, cluster := as.character(cluster)]  
  plot_data <- melt(plot_data, measure.vars = cols, id.vars = "cluster")
  setnames(plot_data, c("variable", "value"), c("gene", "expression"))
  plot_data[, scaled_expression := (expression - min(expression)) / (max(expression) - min(expression)), by = gene]
  x_labels <- as.character(sort(as.numeric(unique(plot_data$cluster))))
  ggplot(data = plot_data) +
    geom_tile(mapping = aes(x = cluster, y = gene, fill = scaled_expression)) +
    geom_text(mapping = aes(x = cluster, y = gene, label = formatC(expression, digits = 2, format = "g")),
              color = "black", size = 5) + 
    scale_x_discrete(limits = x_labels, labels = x_labels, name = element_blank()) +
    scale_fill_gradient2(low = "#0000FF", high = "#FF0000", mid = "#FFFFFF", midpoint = 0.5, name = element_blank()) +
    scale_y_discrete(name = element_blank(), limits = rev(markers)) +
    theme_classic() + theme(text = element_text(size = 20), plot.margin = margin(t = 0, r = 0, b = 0,l = 0,
    unit = "pt"), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 15))
}

heatmaps <- lapply(resolutions, function(res){
  cluster_heatmapper(Cells[, c(res, markers), with = F], res_column = res, cols = markers)
})
names(heatmaps) <- resolutions

################## UMAPS #####################x
set.seed(666)
UMAPS <- lapply(resolutions, function(res){
  umap_coords <- umap(Cells, n_neighbors = 40, min_dist = 0.9, learning_rate = 0.5, init = "random",
                      metric = list("euclidean" = markers, "categorical" = res), n_sgd_threads = 1, n_threads = 1)
  umap_table <- cbind(umap_coords, Cells[, c(markers, res), with = F])
  setnames(umap_table, c("V1", "V2"), c("umap_x", "umap_y"))
  return(umap_table)
})
names(UMAPS) <- resolutions

###################### UMAP Plots by Cluster ##########################
label_dot_plotter <- function(data, x_name, y_name, marker_name)
{
  plot_data <- copy(data)
  data[[marker_name]] <- as.character(data[[marker_name]])	
  my_plot <- ggplot(data) +
    geom_point(mapping = aes_string(x = x_name, y = y_name, color = marker_name)) +
    labs(x = x_name, y = y_name) +
    theme_bw(base_size = 12, base_family = "sans") +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(),
          legend.position = "right", panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.75),
          strip.background = element_rect(colour = NA, fill = NA))
  return(my_plot)
}

UMAP_cluster_plots <- lapply(resolutions, function(res){
  label_dot_plotter(UMAPS[[res]][, c("umap_x", "umap_y", res), with = F], "umap_x", "umap_y", res)
})
names(UMAP_cluster_plots) <- resolutions

###################### UMAP Plots by Marker ##########################
color_dot_plotter <- function(data, x_name, y_name, marker, highcol, lowcol)
{
  plot_data <- copy(data)
  setnames(plot_data, c(x_name, y_name, marker), c("x_name", "y_name", "marker"))
  plot_data[, marker := marker / max(marker)]
  my_plot <- ggplot(plot_data) +
    geom_point(mapping = aes(x = x_name, y = y_name, color = marker)) +
    scale_color_gradient(low = lowcol, high = highcol) +
    theme_bw(base_size = 12, base_family = "sans") +
    labs(x = x_name, y = y_name) +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(),
          legend.position = "right", panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.75),
          strip.background = element_rect(colour = NA, fill = NA))
  return(my_plot)
}

umap_marker_plots <- lapply(resolutions, function(res){
  plots <- lapply(markers, function(marker){
    color_dot_plotter(UMAPS[[res]][, c("umap_x", "umap_y", res, marker), with = F], "umap_x", "umap_y", marker, high_color, low_color)  
  })
  names(plots) <- markers
  return(plots)
})
names(umap_marker_plots) <- resolutions


###################### Boxplots by Marker ##########################

####### Boxplot functions #######
threshold_e = 0.01
round_format_2 <- function(number, threshold = threshold_e)
{
  if (is.nan(number)) return("Nan")
  if (number == 0) return("0")
  if (number >= threshold) return(as.character(round(number, 2)))
  return(paste0(round(as.numeric(sub("e.*", "", formatC(number, format = "e"))), 2), "×10^", sub("^[^-]*", "",
    formatC(number, format = "e"))))
}

boxplotter <- function(data, yaxis_variable_name, yaxis_title, group_variable_name, sample_column_name, bp_title,
  axis_stroke = 0.75) 
{
  local_data <- as.data.frame(copy(data))
  list_return <- list()
  # Make sure no hyphen is in the name of the group levels:
  local_data[, group_variable_name] <- as.character(local_data[, group_variable_name])
  local_data[, group_variable_name] <- gsub(pattern = "-", replacement = " ", local_data[, group_variable_name])
  # Coerce group variable into alphabetically ordered factor:
  local_data[, group_variable_name] <- factor(local_data[, group_variable_name],
    levels = as.character(sort(unique(local_data[, group_variable_name]))))
  # Number of samples in each group:
  uss <- unique(local_data[,c(sample_column_name, group_variable_name)])
  sum_samples_g1 <- sum(uss[,group_variable_name] == levels(uss[,group_variable_name])[1])
  sum_samples_g2 <- sum(uss[,group_variable_name] == levels(uss[,group_variable_name])[2])
  # Two-tailed Wilcoxon test between the two groups:
  g1 <- local_data[local_data[,group_variable_name] == levels(local_data[,group_variable_name])[1], yaxis_variable_name]
  g2 <- local_data[local_data[,group_variable_name] == levels(local_data[,group_variable_name])[2], yaxis_variable_name]
  wilcox_groups <- wilcox.test(g1, g2)
  list_return[["pval"]] <- wilcox_groups$p.value
  # Creating BP:
  list_return[["bp"]] <- ggplot(data = local_data, aes_string(x = group_variable_name, y = yaxis_variable_name,
		label = sample_column_name)) +
    geom_boxplot(fill = NA, width = 0.35, outlier.shape = NA) +
    geom_point(fill = "black") +
    geom_text_repel(data = subset(local_data, get(group_variable_name) == levels(local_data[, group_variable_name])[1]),
      nudge_x = 0.3, direction = "y", hjust = 0, size = 2) +
    geom_text_repel(data = subset(local_data, get(group_variable_name) == levels(local_data[,group_variable_name])[2]),
      nudge_x = 1.7, direction = "y", hjust = 1, size = 2) +
    scale_x_discrete(labels = c(paste0(toupper(abbreviate(levels(local_data[,group_variable_name])[1])),
      " (", sum_samples_g1, ")"), paste0(toupper(abbreviate(levels(local_data[, group_variable_name])[2])),
      " (", sum_samples_g2, ")")), name = element_blank()) +  
    scale_y_continuous(name = yaxis_title, limits = c(0, max(local_data[,yaxis_variable_name])),
      breaks = seq(from = 0, to = round(max(local_data[,yaxis_variable_name])),
      length.out = 10), labels = round(seq(from = 0, to = round(max(local_data[,yaxis_variable_name])), length.out = 10))) +
    ggtitle(bp_title) + labs(subtitle = paste0("p = ", round_format_2(list_return[["pval"]]))) +
    theme_bw(base_size = 8, base_family = "sans") +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.position = "none",
          panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black", size = axis_stroke), strip.background = element_rect(colour = NA,
          fill = NA), plot.subtitle = element_text(size = 8, hjust = 0.5, vjust = 1))
  return(list_return)
}

cell_cluster_boxplotter <- function(data, sample_column, cell_type_column, group_variable_name, axis_stroke = 0.75, correct = TRUE)
{
  cell_types <- sort(unique(data[[cell_type_column]]))
  plot_dataset <- copy(data)[, .SD, .SDcols = c(sample_column, group_variable_name, cell_type_column)]
  setnames(plot_dataset, c(sample_column, group_variable_name, cell_type_column),
    c("sample_column", "x_column", "cluster"))
  plot_dataset[, N_cells := .N, by = c("sample_column", "cluster")]
  plot_dataset[, Total := .N, by = "sample_column"]
  plot_dataset[, Frac := N_cells / Total * 100]
  plot_dataset <- unique(plot_dataset[, .(sample_column, x_column, Frac), keyby = cluster])
  plots <- lapply(cell_types, function(cell_type)
  {
    boxplotter(plot_dataset[cluster == cell_type], "Frac", "Cells %", group_variable_name = "x_column",
      sample_column_name = "sample_column", bp_title = cell_type, axis_stroke = axis_stroke)
  })
  if (correct)
  {
    FDRs <-  p.adjust(sapply(plots, `[[`, 1), method = "BH")
    plots <- lapply(seq_along(plots), function(n){
      pval <- plots[[n]][[1]]
      plt <- plots[[n]][[2]]
      FDR <- FDRs[[n]]
      plt <- plt + labs(subtitle =  paste0("p = ", round_format_2(pval), "\nFDR = ", round_format_2(FDR))) +
        theme(plot.subtitle = element_text(size = 8, hjust = 0.5, vjust = 1))
      return(list(pvalue = pval, FDR_BH = FDR, Plot = plt))        
    })
  }
  names(plots) <- cell_types
  return(plots)
}

################ Make the boxplots
n_categories <- length(unique(Samples[[category_column]]))

if(n_categories == 2){
	cluster_boxplots <- lapply(resolutions, function(res){
		cluster_boxplots <- cell_cluster_boxplotter(Cells, "Metadata_sample_name", res, "category")
		cluster_boxplots <- lapply(cluster_boxplots, function(x){x$Plot})
		return(cluster_boxplots)
	})
	names(cluster_boxplots) <- resolutions
}

################## PDF Output ##################
pdf_w <- 8.27 / 2 # About 1/6th of an A4 (inches)
pdf_h <- 11.69 / 3
pdf_plotter <- function(filename, plot, width = pdf_w, height = pdf_h)
{
	pdf(filename, width = width, height = height, useDingbats = F)
	print(plot)
	dev.off()
}

lapply(resolutions, function(res){
	output_folder <- file.path(output_folder)
	# UMAPs by marker
	umap_marker_output_folder <- file.path(output_folder, "UMAPs_by_Marker", res)
	dir.create(umap_marker_output_folder, recursive = T, showWarnings = F)
	lapply(markers, function(marker){
		pdf_plotter(umap_marker_plots[[res]][[marker]], filename = paste0(umap_marker_output_folder,
			"/UMAP_", marker,"_", res, ".pdf"))
	})
	# UMAPs by cluster
	umap_cluster_output_folder <- file.path(output_folder, "UMAPs_by_Cluster")
	dir.create(umap_cluster_output_folder, recursive = T, showWarnings = F)
	pdf_plotter(UMAP_cluster_plots[[res]], filename = paste0(umap_cluster_output_folder, "/UMAP_cluster_", res, ".pdf"))
	# Heatmaps by cluster
	heatmap_output_folder <- file.path(output_folder, "Heatmaps")
	dir.create(heatmap_output_folder, recursive = T, showWarnings = F)
	pdf_plotter(heatmaps[[res]], filename = paste0(heatmap_output_folder, "/Heatmap_", res, ".pdf"))
	# Boxplots by cluster
	if(n_categories == 2){
		boxplot_output_folder <- file.path(output_folder, "Boxplots", res)  
		dir.create(boxplot_output_folder, recursive = T, showWarnings = F)
		clusters <- as.character(sort(unique(Cells[[res]])))
		lapply(clusters, function(cluster){
			pdf_plotter(cluster_boxplots[[res]][[cluster]], filename =  paste0(boxplot_output_folder,
				"/boxplot_", cluster,"_", res, ".pdf"))
		})
	}	
})
