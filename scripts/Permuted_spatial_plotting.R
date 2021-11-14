###### Setup ######
library(data.table)
library(ggplot2)
source("/opt/Plot_Functions.R")

arguments <- commandArgs(trailingOnly = TRUE)
distance_file_name <- arguments[[1]]
shuffled_distances_file_name <- arguments[[2]]
heterotypic_metadata_file_name <- arguments[[3]]
sample_metadata_file_name <- arguments[[4]]
output_folder <- arguments[[5]]

######## Distances #######
distances <- fread(distance_file_name)
samples <- fread(sample_metadata_file_name)
metadata <- fread(heterotypic_metadata_file_name)

samples[is.na(comparison), comparison := "NA" ]
suppressWarnings(distances[, color := NULL])
suppressWarnings(distances[, comparison := NULL])
distances <- merge(distances, samples, by.x = "Metadata_sample_name", by.y = "sample_name")
distances <- distances[comparison != "NA",]
distances[, pair := paste0(spatial_analysis_cell_type1, "-", spatial_analysis_cell_type2)]

metadata[, pair := paste0(cell_type1, "-", cell_type2)]
pairs <- metadata$pair

pair_colors <- metadata[, color]
names(pair_colors) <- pairs

comps <- sort(unique(samples[!is.na(comparison), comparison]))
n_categories <- length(comps)
######## Shuffled Distances #######
shuffled <- fread(shuffled_distances_file_name)

suppressWarnings(shuffled[, color := NULL])
suppressWarnings(shuffled[, comparison := NULL])
shuffled <- merge(shuffled, samples, by.x = "Metadata_sample_name", by.y = "sample_name")
shuffled <- shuffled[comparison != "NA",]
shuffled[, pair := paste0(spatial_analysis_cell_type1, "-", spatial_analysis_cell_type2)]

################# Median Distances ########################
median_distances_ALL <- unique(distances[, .(median_distance =  median(distance)), by = pair])
median_distances_COMPARISON <- unique(distances[, .(median_distance =  median(distance)), by = c("pair", "comparison")])

median_shuffled_ALL <- unique(shuffled[, .(median_distance =  median(distance)), by = c("pair", "permutation")])
median_shuffled_COMPARISON <- unique(shuffled[, .(median_distance =  median(distance)), by = c("pair", "comparison", "permutation")])

if (n_categories == 2){
	median_distances_DIFFERENCE <- dcast(median_distances_COMPARISON, pair ~ comparison, value.var = "median_distance")
	median_distances_DIFFERENCE <- median_distances_DIFFERENCE[, .(pair, difference = get(comps[[1]]) - get(comps[[2]]))]
	median_shuffled_DIFFERENCE <- dcast(median_shuffled_COMPARISON, pair + permutation ~ comparison, value.var = "median_distance")
	median_shuffled_DIFFERENCE <- median_shuffled_DIFFERENCE[, .(pair, difference = get(comps[[1]]) - get(comps[[2]])), by = permutation]
}

######################### Permutation tests #########################
for (my_comp in comps){
	for(my_pair in pairs){
		median_distances_COMPARISON[pair == my_pair & comparison == my_comp, pval := 
			mean(abs(median_shuffled_COMPARISON[pair == my_pair & comparison == my_comp, median_distance]) > abs(median_distance))]
	}
}
median_distances_COMPARISON[, fdr := p.adjust(pval, "BH"), by = comparison]

for(my_pair in pairs){
	median_distances_ALL[pair == my_pair,
		pval := mean(abs(median_shuffled_ALL[pair == my_pair, median_distance]) > abs(median_distance))]
}
median_distances_ALL[, fdr := p.adjust(pval, "BH")]

if (length(comps) == 2){
	for(my_pair in pairs){
		median_distances_DIFFERENCE[pair == my_pair,
			pval := mean(abs(median_shuffled_DIFFERENCE[pair == my_pair, difference]) > abs(difference))]
	}
	median_distances_DIFFERENCE[, fdr := p.adjust(pval, "BH")]
}

fwrite(x = median_distances_ALL, file = paste0(output_folder,"/median_distances_ALL.csv"))
fwrite(x = median_shuffled_ALL, file = paste0(output_folder,"/median_shuffled_ALL.csv"))
fwrite(x = median_distances_COMPARISON, file = paste0(output_folder,"/median_distances_COMPARISON.csv"))
fwrite(x = median_shuffled_COMPARISON, file = paste0(output_folder,"/median_shuffled_COMPARISON.csv"))
if (n_categories == 2){
	fwrite(x = median_distances_DIFFERENCE, file = paste0(output_folder,"/median_distances_DIFFERENCE.csv"))
	fwrite(x = median_shuffled_DIFFERENCE, file = paste0(output_folder,"/median_shuffled_DIFFERENCE.csv"))
}

################### Permutation Plots ########################
for (my_pair in pairs){
    out_dir <- paste0(output_folder, "/Permutations/", my_pair)
    dir.create(out_dir, recursive = T, showWarnings = F)
	observed_plot_data <- median_distances_ALL[pair == my_pair, .(pair, median_distance, fdr)]
	random_plot_data <- median_shuffled_ALL[pair == my_pair, .(pair, permutation, median_distance)]
	my_plot <- permutation_density(random_plot_data, observed_plot_data, "median_distance", my_pair,
		"Median Distance", "Permutations", pair_colors[[my_pair]])
    pdf_plotter(filename = paste0(out_dir, "/", my_pair, "-all-heterotypic_permutations.pdf"), plot = my_plot)
	if(n_categories >= 1){
		my_plots <- lapply(comps, function(my_comp){
			observed_plot_data <- median_distances_COMPARISON[comparison == my_comp & pair == my_pair,
				.(comparison, pair, median_distance, fdr)]
			random_plot_data <- median_shuffled_COMPARISON[comparison == my_comp & pair == my_pair,
				.(comparison, pair, permutation, median_distance)]
			permutation_density(random_plot_data, observed_plot_data, "median_distance", paste0(my_comp, " ", my_pair),
				"Median Distance", "Permutations", pair_colors[[my_pair]])
		})
		multi_pdf_plotter(filename = paste0(out_dir, "/", my_pair, "-category-heterotypic_permutations.pdf"),
			plots = my_plots, n_row = 2, n_col = 1)
		if(n_categories == 2){
			observed_plot_data <- median_distances_DIFFERENCE[pair == my_pair, .(pair, difference, fdr)]
			random_plot_data <- median_shuffled_DIFFERENCE[pair == my_pair, .(pair, permutation, difference)]
			my_plot <- permutation_density(random_plot_data, observed_plot_data, "difference", my_pair,
				"Median Distance Difference", "Permutations", pair_colors[[my_pair]])
			pdf_plotter(filename = paste0(out_dir, "/", my_pair, "-difference-heterotypic_permutations.pdf"), plot = my_plot)
		}
	}
}	
