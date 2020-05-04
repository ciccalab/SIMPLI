script_folder = "$baseDir/scripts"

// Paths to activate a python3 virtual environment with imctools
python_environment_folder = "$baseDir/venv"
python_path = "/share/apps/python-3.7.2-shared/bin"
py_library_path = "/share/apps/python-3.7.2-shared/lib"

// Input parameters: the default values will be moved to a test profile nextflow configuration
params.output_folder ="$baseDir/example_output"
params.data_folder = "$baseDir/example_data"
params.channel_metadata = "$params.data_folder/channel_metadata.csv" 
params.raw_metadata_file = "$baseDir/example_data/sample_metadata.csv" 
params.tiff_type = "single" 

/* Reads the raw metadata file line by line to extract the sample metadata for the raw IMC acquisition files.
   It expects an header line and it extracts the following fields into the sample_metadata channel:
   - sample_name
   - roi_name
   - raw_path -> Converted to a file type 
   - tiff_path
*/

Channel
    .fromPath(params.raw_metadata_file)
    .splitCsv(header:true)
    .map{row -> tuple(row.sample_name, row.roi_name, file(row.raw_path), row.tiff_path)}
    .set{sample_metadata_raw}


/* For each aquisition spefied in the $raw_metadata_file:
    - Activates a python3 virtual environment with the imctools modules
    - Creates the $tiff_path directory if it doesn't already exist
    - Extracts the raw tiff files into the working directory
    - Generates the raw tiff .csv metadata file in the working directory
    - Copies the output files into $tiff_path 
  For ome.tiff output the channels in the .ome are in the same order they are
  written in the metadata.
*/

process convert_raw_data_to_tiffs {

    publishDir "$tiff_path/raw", mode:'copy', overwrite: true

    input:
        set sample_name, roi_name, raw_path, tiff_path from sample_metadata_raw

    output:
        file "$sample_name*raw*tiff" into raw_tiff_images
        file "${sample_name}_raw_tiff_metadata.csv" into raw_tiff_metadata_by_sample
    script:
    """
    export PATH="$python_path:$PATH"
    export LD_LIBRARY_PATH="$py_library_path:$LD_LIBRARY_PATH"
    source $python_environment_folder/bin/activate
    mkdir -p $tiff_path
    mkdir -p $sample_name
    python3 $script_folder/tiff_extracter.py \\
        $sample_name \\
        '$roi_name' \\
        $raw_path \\
        $params.tiff_type \\
        . \\
        $params.channel_metadata \\
        ${sample_name}_raw_tiff_metadata.csv > $sample_name/extract_log.txt 2>&1
    """
}

/* Collects all the raw_tiff_metadata_by_sample metadata files:
    - Concatenates them into a single $raw_tiff_metadata file
    - Removes extra header lines from the middle of the file, as each of the
        original files starts with an header line.
*/

process collect_raw_tiff_metadata {

    publishDir "$params.output_folder", mode:'copy', overwrite: true
    
    input:
        file metadata_list from raw_tiff_metadata_by_sample.collect()
    
    output:
        file "raw_tiff_metadata.csv" into raw_tiff_metadata

    script:
    """
    cat $metadata_list > raw_tiff_metadata.csv
    sed -i '1!{/sample_name,roi_name,metal,label,raw_tiff_file_name/d;}' raw_tiff_metadata.csv
    """
}

Channel
    .fromPath(params.raw_metadata_file)
    .splitCsv(header:true)
    .map{row -> tuple(row.sample_name, row.roi_name, row.tiff_path)}
    .set{sample_metadata_normalized}

process normalize_tiff {
    container = 'library://michelebortol/default/imcpipeline-rbioconductor:test'
    containerOptions = "--bind $script_folder:/opt"
    publishDir "$tiff_path/normalized", mode:'copy', overwrite: true
    
    input:
        set sample_name, roi_name, tiff_path from sample_metadata_normalized
        file raw_metadata from raw_tiff_metadata

    output:
        file "$sample_name*normalized*tiff" into normalized_tiff_images
        file "${sample_name}_normalized_tiff_metadata.csv" into normalized_tiff_metadata_by_sample
    script:
    """
    mkdir -p $sample_name
    Rscript /opt/tiff_normalizer.R \\
        $sample_name \\
        $raw_metadata\\
        $params.tiff_type \\
        . \\
        ${sample_name}_normalized_tiff_metadata.csv > $sample_name/normalization_log.txt 2>&1
    """
}

process collect_normalized_tiff_metadata {

    publishDir "$params.output_folder", mode:'copy', overwrite: true
    
    input:
        file metadata_list from normalized_tiff_metadata_by_sample.collect()
    
    output:
        file "normalized_tiff_metadata.csv" into normalized_tiff_metadata

    script:
    """
    cat $metadata_list > normalized_tiff_metadata.csv
    sed -i '1!{/sample_name,metal,label,normalized_file/d;}' normalized_tiff_metadata.csv
    """
}
