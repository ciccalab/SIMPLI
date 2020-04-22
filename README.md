# IMC_Pipeline

## Requirements
- nextflow
- python3, pip3, venv
- imctools

## Setup
1. Setup the python3 virtual environment using the pip3_requirements.txt file
2. Install nextlow

## Run

### Nextflow profiles
- local
- cluster: for execution in a SGE environment

### Input
- Raw IMC aquisition data. (.mcd or .txt formats accepted)
- Channel metadata table, .csv file with an header and the following columns:
  - channel_metal = Metal the channel is associated to. Used to select which channels to extract from the IMC acquisition
  - channel_label = Marker associated to the meta. Used to give meaningful names to the single tiff files. 

- Raw IMC acquisition metadata, .csv with an header and the following columns:
  - sample_name = Unique sample identifier
  - roi_name = Name of the ROI to extract the acquisition for
  - raw_path = Path to the raw IMC acquisition file
  - tiff_path = Path where to output the single tiff images

### Output
- Single raw tiff images (one for each of the selected channels) 
- Single raw tiff .csv metadata file.



