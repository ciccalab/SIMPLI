####### SETUP #######
rm(list = ls())
library(dplyr)
library(data.table)

arguments <- commandArgs(trailingOnly = TRUE)
print(arguments)
in_file_names <- arguments[1:(length(arguments) -1)]
out_file_name <- arguments[[length(arguments)]]

clustered_cell_tables <- lapply(in_file_names, fread)
clustered_cells <- rbindlist(clustered_cell_tables, use.names = T, fill = T)
clustered_cells[, c("nGene", "nUMI") := NULL]

fwrite(clustered_cells, out_file_name)

