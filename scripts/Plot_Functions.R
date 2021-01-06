library(data.table)
library(ggplot2)
library(ggrepel)
library(gridExtra)
library(EBImage)

######## Heatmaps #######
heatmapper <- function(plot_dataset, res_column, cols, high_color, mid_color, low_color)
{
  plot_data <- copy(plot_dataset)
  setnames(plot_data, res_column, "cluster")
  plot_data <- plot_data[, lapply(.SD, mean), by = cluster, .SDcols = cols]
  plot_data[, cluster := as.character(cluster)]  
  plot_data <- melt(plot_data, measure.vars = cols, id.vars = "cluster")
  setnames(plot_data, c("variable", "value"), c("gene", "expression"))
  plot_data[, scaled_expression := (expression - min(expression)) / (max(expression) - min(expression)), by = gene]
  x_labels <- as.character(sort(as.numeric(unique(plot_data$cluster))))
  if(any(is.na(as.numeric(unique(plot_data$cluster))))){
    x_labels <- sort(unique(plot_data$cluster))
  }
  ggplot(data = plot_data) +
    geom_tile(mapping = aes(x = cluster, y = gene, fill = scaled_expression)) +
    geom_text(mapping = aes(x = cluster, y = gene, label = formatC(expression, digits = 2, format = "g")),
              color = "black", size = 2) + 
    scale_x_discrete(limits = x_labels, labels = x_labels, name = element_blank()) +
    scale_fill_gradient2(low = low_color, high = high_color, mid = mid_color, midpoint = 0.5, name = element_blank()) +
    scale_y_discrete(name = element_blank(), limits = rev(cols)) + coord_flip() +
    theme_classic() + theme(text = element_text(size = 8), plot.margin = margin(t = 0, r = 0, b = 0,l = 0,
    unit = "pt"), plot.title = element_text(hjust = 0.5), axis.text.x = element_text(size = 8, angle = 90))
}

####### Scatterplots #######
label_dot_plotter <- function(data, x_name, y_name, label_name, plot_title = NULL)
{
  plot_data <- copy(data)
  data[[label_name]] <- as.character(data[[label_name]])	
  my_plot <- ggplot(data) +
    geom_point(mapping = aes_string(x = x_name, y = y_name, color = label_name), size = 0.25) +
    labs(x = x_name, y = y_name) +
    theme_bw(base_size = 8, base_family = "sans") +
	labs(title = plot_title) +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(),
          legend.position = "right", panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.75),
          strip.background = element_rect(colour = NA, fill = NA))
  return(my_plot)
}

gradient_dot_plotter <- function(data, x_name, y_name, marker, highcol, lowcol)
{
  plot_data <- copy(data)
  setnames(plot_data, c(x_name, y_name, marker), c("x_name", "y_name", "marker"))
  plot_data[, marker := marker / max(marker)]
  my_plot <- ggplot(plot_data) +
    geom_point(mapping = aes(x = x_name, y = y_name, color = marker), size = 0.25) +
    scale_color_gradient(low = lowcol, high = highcol) +
    theme_bw(base_size = 8, base_family = "sans") +
    labs(x = x_name, y = y_name, title = marker) +
    theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(),
          legend.position = "right", panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black", size = 0.75),
          strip.background = element_rect(colour = NA, fill = NA))
  return(my_plot)
}

######### Barplots  ############
Barplotter <- function(plot_dataset, x_column, annotation_column, color_column, annotation_colors, y_axis_title)
{
	plot_data <- copy(plot_dataset)
	if(x_column == annotation_column){ # Plot by sample
		setnames(plot_data, c(x_column, color_column), c("annotation_column", "color_column"))
		annotations <- sort(unique(plot_data$annotation_column))
		plot_title <- "Samples"
	} else{ # Plot by category
		setnames(plot_data, c(x_column, annotation_column, color_column), c("x_column", "annotation_column", "color_column"))
		annotations <- sort(unique(plot_data$annotation_column))
		annotations <- sapply(annotations, function(annotation){
			paste0(annotation, " (", length(unique(plot_data[annotation_column == annotation, x_column])), ")")
		})	
		plot_title <- annotation_column
	}
	names(annotations) <- sort(unique(plot_data$annotation_column))
	plot_data[, count := .N, by = c("color_column", "annotation_column")]	
	plot_data[, total := .N, by = annotation_column]
	plot_data[, percent := count / total * 100]
	plot_data <- unique(plot_data[, .(annotation_column, color_column, percent)])
	plot_data <- plot_data[, percent_text := paste0(round(percent), " %")]
	plot_data[, color_column := factor(color_column, levels = names(annotation_colors))]
	plot_data[order(color_column, decreasing = T), percent_y := cumsum(percent) - percent * 0.5, by = annotation_column]
	plot_data <- plot_data[, .SD, keyby = color_column]

	ggplot(data = plot_data) +
		geom_bar(aes(x = annotation_column, y = percent, fill = color_column), stat = "identity", color = "#000000") +
		geom_text(aes(x = annotation_column, y = percent_y, label = percent_text), size = 2) +
	    scale_x_discrete(labels = annotations, breaks = sort(unique(plot_data$annotation_column)), name = element_blank()) +
		scale_y_continuous(name = y_axis_title) +
		scale_fill_manual(values = annotation_colors, name = element_blank(), labels = names(annotation_colors)) +
		labs(title = plot_title) +
		theme_bw(base_size = 8, base_family = "sans") +
		theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.position = "right",
			panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
			axis.line = element_line(colour = "black", size = 0.75), strip.background = element_rect(colour = NA, fill = NA))
}

####### Boxplots #######
threshold_e = 0.01
round_format_n <- function(numbers, n = 2, threshold = threshold_e)
{
	sapply(numbers, function(number){
		if (is.nan(number)) return("Nan")
		if (number == 0) return("0")
		if (number >= threshold) return(as.character(round(number, n)))
		return(paste0(round(as.numeric(sub("e.*", "", formatC(number, format = "e"))), n), "×10^", sub("^[^-]*", "",
			formatC(number, format = "e"))))
	})
}

boxplotter <- function(data, y_axis_variable_name, y_axis_title, group_variable_name, sample_column_name, bp_title,
	color_variable_name, color_list, axis_stroke = 0.75) 
{
	local_data <- as.data.frame(copy(data))
	list_return <- list()
	# Coerce group variable into alphabetically ordered factor:
	local_data[, group_variable_name] <- factor(local_data[, group_variable_name],
	levels = as.character(sort(unique(local_data[, group_variable_name]))))
	# Number of samples in each group:
	uss <- unique(local_data[,c(sample_column_name, group_variable_name)])
	sum_samples_g1 <- sum(uss[,group_variable_name] == levels(uss[,group_variable_name])[1])
	sum_samples_g2 <- sum(uss[,group_variable_name] == levels(uss[,group_variable_name])[2])
	# Two-tailed Wilcoxon test between the two groups:
	g1 <- local_data[local_data[,group_variable_name] == levels(local_data[,group_variable_name])[1], y_axis_variable_name]
	g2 <- local_data[local_data[,group_variable_name] == levels(local_data[,group_variable_name])[2], y_axis_variable_name]
	wilcox_groups <- try(wilcox.test(g1, g2))
	list_return[["pval"]] <- ifelse(class(wilcox_groups) != "try-error", wilcox_groups$p.value, 1)
	# Creating Boxplot:
	y_ticks <- seq(from = 0, to = max(local_data[,y_axis_variable_name]), length.out = 10)
	list_return[["bp"]] <- ggplot(data = local_data, aes_string(x = group_variable_name, y = y_axis_variable_name,
		label = sample_column_name)) +
	geom_boxplot(fill = NA, width = 0.35, outlier.shape = NA) +
	geom_point(mapping = aes_string(color = color_variable_name), size = 1) +
	geom_text_repel(data = subset(local_data, get(group_variable_name) == levels(local_data[, group_variable_name])[1]),
		nudge_x = 0.3, direction = "y", hjust = 0, size = 2, segment.size = 0.5) +
	geom_text_repel(data = subset(local_data, get(group_variable_name) == levels(local_data[, group_variable_name])[2]),
		nudge_x = 1.7, direction = "y", hjust = 1, size = 2, segment.size = 0.5) +
	scale_x_discrete(labels = c(paste0(toupper(abbreviate(levels(local_data[,group_variable_name])[1])),
		" (", sum_samples_g1, ")"), paste0(toupper(abbreviate(levels(local_data[, group_variable_name])[2])),
		" (", sum_samples_g2, ")")), name = element_blank()) +  
	scale_y_continuous(name = y_axis_title, limits = c(0, max(local_data[,y_axis_variable_name])), breaks = y_ticks,
		labels = round_format_n(y_ticks, 1)) +
	scale_color_manual(values = color_list) +		   
	ggtitle(bp_title) + labs(subtitle = paste0("p = ", round_format_n(list_return[["pval"]]))) +
	theme_bw(base_size = 8, base_family = "sans") +
	theme(plot.title = element_text(hjust = 0.5), legend.title = element_blank(), legend.position = "none",
			panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
			axis.line = element_line(colour = "black", size = axis_stroke), strip.background = element_rect(colour = NA,
			fill = NA), plot.subtitle = element_text(size = 8, hjust = 0.5, vjust = 1))
	return(list_return)
}

list_boxplotter <- function(data, sample_column, type_column, group_variable_name, y_column, y_axis_title, color_variable_name,
	color_list, axis_stroke = 0.75, correct = TRUE)
{
  types <- sort(unique(data[[type_column]]))
  plot_dataset <- copy(data)[, .SD, .SDcols = c(sample_column, group_variable_name, type_column, y_column, color_variable_name)]
  setnames(plot_dataset, c(sample_column, group_variable_name, type_column, y_column, color_variable_name),
    c("sample_column", "x_column", "type_column", "y_column", "color_column"))
  plot_dataset <- unique(plot_dataset[, .(sample_column, x_column, y_column, color_column, type_column)])
  plots <- lapply(types, function(type)
  {
    boxplotter(plot_dataset[type_column == type], "y_column", y_axis_title, group_variable_name = "x_column",
      sample_column_name = "sample_column", "color_column", color_list, bp_title = type, axis_stroke = axis_stroke)
  })
  if (correct)
  {
    FDRs <-  p.adjust(sapply(plots, `[[`, 1), method = "BH")
    plots <- lapply(seq_along(plots), function(n){
      pval <- plots[[n]][["pval"]]
      plt <- plots[[n]][["bp"]]
      FDR <- FDRs[[n]]
      plt <- plt + labs(subtitle =  paste0("p = ", round_format_n(pval), "\nFDR = ", round_format_n(FDR))) +
        theme(plot.subtitle = element_text(size = 8, hjust = 0.5, vjust = 1))
      return(list(pvalue = pval, FDR_BH = FDR, Plot = plt))        
    })
  }
  names(plots) <- types
  return(plots)
}

####### Image Processing #######
load_image <- function(filename)
{
	Img <- readImage(filename, "tiff", all = FALSE)
	if (length(dim(Img)) == 3) {return(Img[,,1])}
	return(Img)
}

color_imager <- function(pic_list, color_list)
{
	pics <- lapply(seq_along(pic_list), function(n){
		color_ratio <- col2rgb(color_list[[n]]) / 255
		rgbImage(pic_list[[n]] * color_ratio[[1]], pic_list[[n]] * color_ratio[[2]], pic_list[[n]] * color_ratio[[3]])
	})
	Reduce(`+`, pics)
 }

cell_overlayer <- function(picture, mask, cell_list, color_list)
{
	pics <- lapply(seq_along(cell_list), function(n){
		cell_mask <- mask * (2^16 - 1)
		pic <- Image(dim = c(dim(cell_mask),3), colormode = "Color")
		cell_mask[which(!(cell_mask %in% cell_list[[n]]), arr.ind = T)] <- 0
		paintObjects(cell_mask, pic, col = color_list[[n]], thick = T)
	})
	pic <- Reduce(`+`, pics)
	return(pic + ifelse(any(picture) > 0, rgbImage(picture, picture, picture), 0))
}

############# Density Plot
marker_histogram_density <- function(data, marker_col, plot_title, x_label, y_label, geom = "density",
  color = "black", fill = "grey", alpha = 0.3, log_x = T, log_y = T)
{
  x_n_breaks <- seq(0, 1, 0.05)
  x_log_breaks <- c(0.00001, 0.0001, 0.001, 0.01, 0.1, 0.5, 1)
  plot_dataset <- copy(data)
  setnames(plot_dataset, marker_col, "marker_col")
  plot_dataset[, marker_col := marker_col + (1* 10^-5)]
  plt <- ggplot(plot_dataset, aes(x = marker_col)) +
    theme_classic() + theme(plot.margin = margin(t = 0, r = 0, b = 0,l = 0, unit = "pt"),
    axis.text = element_text(size = 15), axis.title = element_text(size = 15),
    legend.text = element_text(size = 15))
  if(log_x){
    plt <- plt + scale_x_continuous(name = paste("Log ", x_label), trans = "log10",
      breaks = x_log_breaks, labels = x_log_breaks)
  }else{
    plt <- plt + scale_x_continuous(name = x_label, breaks = x_n_breaks, labels = x_n_breaks)
  }
  if(log_y){
    plt <- plt +
      stat_bin(aes(y = ..count.. + 1), binwidth = NULL, bins = 100, fill = fill, alpha = alpha,
        color = color, geom = geom) +
      scale_y_continuous(name = paste("Log ", y_label), trans = "log10")
  }else{
    plt <- plt + stat_bin(aes(y = ..count..), binwidth = NULL, bins = 100, fill = fill,
        alpha = alpha, color = color, geom = geom) +
      scale_y_continuous(name = y_label)
  }
  if(log_x | log_y){
    plt <- plt + annotation_logticks(sides = ifelse(log_x & log_y, "lb", ifelse(log_x, "b", "l")))
  }
  plt <- plt + labs(title = plot_title)
  return(plt)
}

################## PDF Output ##################
single_pdf_w <- 8.27 / 2 # About 1/6th of an A4 (inches)
single_pdf_h <- 11.69 / 3
multi_pdf_w <- 8.27 # A4 Portrait
multi_pdf_h <- 11.69

pdf_plotter <- function(filename, plot, width = single_pdf_w, height = single_pdf_h)
{
	pdf(filename, width = width, height = height, useDingbats = F)
	print(plot)
	dev.off()
}

multi_pdf_plotter <- function(filename, plots, n_col = 2, n_row = 3, width = multi_pdf_w, height = multi_pdf_h)
{
	pdf(NULL)
	arranged_grobs <- marrangeGrob(grobs = plots, nrow = n_row, ncol = n_col, top = NULL)
	dev.off()
	pdf(filename, width = width, height = height, useDingbats = F)
	print(arranged_grobs)
    dev.off()
}
