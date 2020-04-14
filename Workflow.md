# Workflow

## Raw file conversion
### Input
- Raw Metadata table (.csv file with input/output paths and sample metadata)
- .txt or .mcd data files
### Output
- Single tiff raw data 	
- Single tiff raw metadata table (.csv file with paths to the single tiff files)	
### Tools
- Python3 script using the [imctools](https://github.com/BodenmillerGroup/imctools) python3 library.

## Tiff file preprocessing
### Input
- Single tiff raw data 	
- Single tiff metadata table (.csv file with paths to the single tiff files)
- Cellprofiler3 pipeline
### Output
- Cleaned single tiff files.
- Tissue mask tiff files.
- Processed tiff and tissue mask metadata table (.csv file with paths to the single tiff files)	
### Tools
- CellProfiler3 to run a custom pipeline provided by the user.

## Pixel Analysis
### Input
- Cleaned single tiff files.
- Mask tiff files.
- Processed tiff and mask metadata table (.csv file with paths to the single tiff files)
- .txt file with the measurements to perform.
### Output
- Pixel analysis measurements
### Tools
- R script (EBImage package)

## Cell Segmentation
### Input
- Cleaned single tiff files.
- Mask tiff files.
- Processed tiff and mask metadata table (.csv file with paths to the single tiff files)
- Cellprofiler3 pipeline
### Output
- Cell mask tiff files
- Single Cell data table
- Cell mask metadata table (.csv file with paths to the single tiff files)	
### Tools
- CellProfiler3 to run a custom pipeline provided by the user.

## Cell population identification
### Input
- Single Cell data table
- .txt file with population names and thresholds
- Mask tiff files (optional).
- Cleaned single tiff files (optional).
### Output
- Annotated Single Cell data table	
### Tools
- R script (EBImage package, optional)

## Cell clustering identification
### Input
- Annotated Single Cell data table	
- .txt file with population names and markers to use for clustering
### Output
- Annotated Single Cell data table with clusters
### Tools
- R script (Seurat package v 2.3.4)
