# IMC_Pipeline

## Requirements
- [nextflow](https://www.nextflow.io/)
- python3 (you may need `export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/path/to/python/lib/`, pip3 (updated to latest), venv
- [imctools](https://github.com/BodenmillerGroup/imctools)
- [singularity](https://sylabs.io/docs/)

## Setup
1. Setup the python3 virtual environment using the pip3_requirements.txt file
```
python3 -m venv venv
source venv/bin/activate
pip3 install -r pip3_requirements.txt
deactivate
```
2. Install nextlow

## Run

### Nextflow profiles
- local
- cluster: for execution in a SGE environment

### Input
- Raw IMC aquisition data. (.mcd or .txt formats accepted)
- Channel metadata table, .csv file with an header and the following columns:
  - channel_metal = Metal the channel is associated to. Used to select which channels to extract
      from the IMC acquisition
  - channel_label = Marker associated to the meta. Used to give meaningful names to the single
      tiff files. 
- Output format flag: 'single' for single tiffs or 'ome' for a ome.tiff"

- Raw IMC acquisition metadata, .csv with an header with the following columns:
  - sample_name = Unique sample identifier
  - roi_name = Name of the ROI to extract the acquisition for
  - raw_path = Path to the raw IMC acquisition file
  - tiff_path = Path where to output the single tiff images
  
- CellProfiler3 preprocessing pipeline

- Area measurements metadata, .csv with an header with the two following columns:
  - marker = marker or combination of markers to measure (markers can be combined with the logical operators: "&", "|", "!" and "()")
  - main_marker = marker or combination of markers to measure
  The ratio of marker/main_marker area will be reported as a percentage

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
  - One .ome.tiff file: the order of channels in the .ome.tiff file is the same as the order they
    are reported in the metadata file.
- Preprocessed tiff .csv metadata file.

- Area measurements, and area ratios between user defined markers or their combinations.


