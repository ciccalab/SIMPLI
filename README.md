# SIMPLI: Single-cell Identification from Multiplexed Images

SIMPLI: Single-cell Identification from Multiplexed Images is a image processing pipeline for the analysis of highly multiplexed histological imaging data. It performs measurements of areas positive for user defined combination of markers and single-cell data analysis. The single cell data analysis includes: cell segmentation, expression based cell type identification and unsupervised clustering.

<img src="assets/SIMPLI.png" width="800" height="600">

## Requirements
- [nextflow](https://www.nextflow.io/)
- [singularity](https://sylabs.io/docs/)

## Setup
1. Install singularity.
2. Install nextlow.

## Run

### Nextflow profiles
- local
- cluster: for execution in a SGE environment

### Input
- Raw IMC aquisition data. (.mcd or .txt formats accepted)
- Channel metadata table, .csv file with an header and the following columns:
  - channel_metal = Metal the channel is associated to. Used to select which channels to extract
      from the IMC acquisition.
  - channel_label = Marker associated to the meta. Used to give meaningful names to the single
      tiff files. 
- Output format flag: 'single' for single tiffs or 'ome' for a ome.tiff".

- Raw IMC acquisition metadata, .csv file with an header with the following columns:
  - sample_name = Unique sample identifier.
  - roi_name = Name of the ROI to extract the acquisition for.
  - raw_path = Path to the raw IMC acquisition file.
  - color = Color used to represent the sample (Used for boxplots, accepts R color strings and #FFFFFF hexadecimal)
  - ... = All columns after the 4th are considered categorical variables for grouping the samples.
  
 - CellProfiler3 pipelines:
   - Image preprocessing
   - Cell segmentation

- Area measurements metadata, .csv file with an header with the two following columns:
  - marker = marker or combination of markers to measure (markers can be combined with the logical operators: "&", "|", "!" and "()").
  - main_marker = marker or combination of markers to measure.
  The ratio of marker/main_marker area will be reported as a percentage.
  
- Cell type and expression threshold metadata, .csv file with an header with the following columns:
  - cell_type = name of the cell type to identify.
  - marker = name of the marker to use as threshold (column name of the cellprofiler exported single cell data).
  - threshold = value to use as threshold for a cell to be of this cell type.
  
- Cell clustering metadata, .csv file with an header with the following columns: 
  - cell_type = cell type to cluster (from the cell type and expression threshold metadata).
  - markers = markers to include in the clusterin (column names of the cellprofiler exported single cell data).
  - resolutions = resolutions values at which to cut the nearest neighbour graph.
  Values in the markers and resolution columns should be separated by '@'.

### Output
- Raw tiff images in either:
  - Single tiff files: one for each of the selected channels.
  - One .ome.tiff file: the order of channels in the .ome.tiff file is the same as the order they
    are reported in the metadata file.
- Single raw tiff .csv metadata file.

- Normalized tiff images in either:
  - Single tiff files: one for each of the selected channels.
  - One .ome.tiff file: the order of channels in the .ome.tiff file is the same as the order they
    are reported in the metadata file.
- Normalized tiff .csv metadata file.

- Preprocessed tiff images in either:
  - Single tiff files: one for each of the selected channels.
- Preprocessed tiff .csv metadata file.

- Area measurements, and area ratios between user defined markers or their combinations.

- Cell masks in uint16 tiff format.

- Annotated cell measurements in .csv format.
  Cell measuremetnts annotated with the clusters each cell belongs to at the different levels of
  resolution selected by the user.

- Plots at three levels:
  - Areas:
    - Boxplots: The % of positive area for each user-defined combination of marker is compared between the two groups (2 groups of samples only).
  - Cell type:
    - Overlay of all cells coloured by cell type.
    - Barplots: with the percentage of each cell type in the category type or sample.
    - Boxplots: The fraction of cells in each main cell type is compared between the two groups (2 groups of samples only).  
  - Cell clusters:  
    - UMAPs: UMAPs coloured by cluster, and by marker.
    - Heatmaps: Heatmaps of mean marker intensities for every cluster.
    - Boxplots: For each cluster the fraction of cells in its main cell type is compared between the two groups (2 groups of samples only).  
