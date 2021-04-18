####### SETUP #######
rm(list = ls())
library(data.table)

args <- commandArgs(trailingOnly = TRUE)
cells_file_name <- args[[1]]
threshold_metadata_file_name <- args[[2]]
out_file_name <- args[[3]]

############ Load the cell data ##########
cells <- fread(cells_file_name)

########### Load the marker data  ##################
thresholds <- fread(threshold_metadata_file_name)

############## Threshold the cells
thresholded_cells <- copy(cells)
thresholded_cells[, index := .I]
for (pop in unique(thresholds$cell_type)){
  tresh_col <- paste0(pop, "_Thresholded")
  thresholded_cells[, c(tresh_col) := ifelse(cell_type == pop, "UNASSIGNED", "NA")]
  for (instr in thresholds[cell_type == pop, threshold_expression]){
    print(thresholded_cells[ eval(parse(text = eval(tresh_col))) != "NA", .N])
    thresholded_cells[eval(parse(text = eval(tresh_col))) != "NA",  c(tresh_col) :=
      ifelse(eval(parse(text = instr)),
        thresholds[cell_type == pop & threshold_expression == instr, population_name],
        eval(parse(text = eval(tresh_col)))), by = index]
  }
}
thresholded_cells[, index := NULL]

########### Output
fwrite(file = out_file_name, x = thresholded_cells)
