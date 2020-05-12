script_folder = "$baseDir/scripts"
cp3_pipeline_folder = "$baseDir/cell_profiler_pipelines"

// Paths to activate a python3 virtual environment with imctools
python_environment_folder = "$baseDir/venv"
python_path = "/share/apps/python-3.7.2-shared/bin"
py_library_path = "/share/apps/python-3.7.2-shared/lib"

// Input parameters: the default values will be moved to a test profile nextflow configuration
params.output_folder ="$baseDir/example_output"
params.data_folder = "$baseDir/example_data"
params.channel_metadata = "$params.data_folder/channel_metadata.csv" 
params.raw_metadata_file = "$baseDir/example_data/sample_metadata.csv" 
params.tiff_type = "ome" 
params.cp3_preprocessing_cppipe = "ome_preprocessing_example.cppipe"

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
    - Copies the output files into "$tiff_path/raw" 
  For ome.tiff output the channels in the .ome are in the same order they are
  written in the metadata.
*/

process convert_raw_data_to_tiffs {

    label 'big_memory'
    publishDir "$tiff_path/raw", mode:'copy', overwrite: true

    input:
        set sample_name, roi_name, raw_path, tiff_path from sample_metadata_raw

    output:
        file "$sample_name*raw*tiff" into raw_tiff_images
        file "${sample_name}-raw_tiff_metadata.csv" into raw_tiff_metadata_by_sample
    script:
    """
    export PATH="$python_path:$PATH"
    export LD_LIBRARY_PATH="$py_library_path:$LD_LIBRARY_PATH"
    source $python_environment_folder/bin/activate
    python3 $script_folder/tiff_extracter.py \\
        $sample_name \\
        '$roi_name' \\
        $raw_path \\
        $params.tiff_type \\
        . \\
        $params.channel_metadata \\
        ${sample_name}-raw_tiff_metadata.csv > extract_log.txt 2>&1
    """
}

raw_tiff_metadata_by_sample.into{raw_tiff_metadata_to_collect; raw_tiff_metadata_to_normalize}

/* Collects all the raw_tiff_metadata_by_sample metadata files:
    - Concatenates them into raw_tiff_metadata.csv 
    - Removes extra header lines from the middle of the file, as each of the
        original files starts with an header line.
    - Copies the metadata to "$params.output_folder/raw_tiff_metadata.csv" 
*/

process collect_raw_tiff_metadata {

    label 'small_memory'
    publishDir "$params.output_folder", mode:'copy', overwrite: true
    
    input:
        file metadata_list from raw_tiff_metadata_to_collect.collect()
    
    output:
        file "raw_tiff_metadata.csv" into raw_tiff_metadata

    script:
    """
    cat $metadata_list > raw_tiff_metadata.csv
    sed -i '1!{/sample_name,roi_name,metal,label,raw_tiff_file_name/d;}' raw_tiff_metadata.csv
    """
}

raw_tiff_metadata_to_normalize
    .map{file ->
            def key = file.name.toString().tokenize('-').get(0)
            return tuple(key, file)
         }
    .groupTuple()
    .set{normalize_tiff_metadata}

/* For each sample:
    - Performs 99th percentile normalization and scaling
    - Copies the results to "$params.output_folder/$sample_name/normalized"
  For ome.tiff output the channels in the .ome are in the same order they are
  written in the metadata.
*/

process normalize_tiff {
    
    label 'big_memory'
    container = 'library://michelebortol/default/imcpipeline-rbioconductor:test'
    containerOptions = "--bind $script_folder:/opt"
    
    publishDir "$params.output_folder/$sample_name/normalized", mode:'copy', overwrite: true
    
    input:
        set sample_name, file(raw_metadata_file) from normalize_tiff_metadata

    output:
        file "$sample_name*normalized*tiff" into normalized_tiff_images
        file "${sample_name}-normalized_tiff_metadata.csv" into normalized_tiff_metadata_by_sample
        file "${sample_name}-cp3_normalized_tiff_metadata.csv" into cp3_normalized_tiff_metadata_by_sample
    script:
    """
    mkdir -p $sample_name
    Rscript /opt/tiff_normalizer.R \\
        $sample_name \\
        $raw_metadata_file \\
        $params.tiff_type \\
        ./ \\
        ${sample_name}-normalized_tiff_metadata.csv \\
        ${sample_name}-cp3_normalized_tiff_metadata.csv > normalization_log.txt 2>&1
    """
}

/* Collects all the normalized_tiff_metadata_by_sample metadata files:
    - Concatenates them into normalized_tiff_metadata.csv
    - Removes extra header lines from the middle of the file, as each of the
        original files starts with an header line.
    - Copies the metadata to "$params.output_folder/normalized_tiff_metadata.csv" 
*/

process collect_normalized_tiff_metadata {

    label 'small_memory'
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

cp3_normalized_tiff_metadata_by_sample
    .map { file ->
            def key = file.name.toString().tokenize('-').get(0)
            return tuple(key, file)
            }
    .groupTuple()
    .set{ cp3_preprocessing_metadata  }

process cell_profiler_image_preprocessing {

    label 'big_memory'
    container = 'library://michelebortol/default/cellprofiler3_imcplugins:example'
    containerOptions = "--bind $cp3_pipeline_folder:/mnt"

    publishDir "$params.output_folder/$sample_name/preprocessed", mode:'copy', overwrite: true
                                                                                                
    input:
        set sample_name, file(cp3_normalized_metadata) from cp3_preprocessing_metadata 

    output:
        file "*-Preprocessed.tiff" into preprocessed_tiff_files
        file "${sample_name}-preprocessed_metadata.csv" into preprocessed_tiff_metadata_by_sample
    script:
    """
    cellprofiler \\
        --run-headless \\
        --data-file $cp3_normalized_metadata \\
        --pipeline /mnt/$params.cp3_preprocessing_cppipe \\
        --plugins-directory /opt/CellProfiler/plugins/ \\
        --image-directory /data \\
        --output-directory ./ \\
        --log-level DEBUG \\
        --temporary-directory ./tmp > cp3_preprocessing_log.txt 2>&1

    echo "sample_name,label,preprocessed_file_name" > "${sample_name}-preprocessed_metadata.csv"
    find "\$(pwd)" -name "*-Preprocessed.tiff" > filename.csv
    sed "s@.*-\\(.*\\)-Preprocessed.*tiff@\\1@" filename.csv > label.csv
    paste -d , label.csv filename.csv >> "${sample_name}-preprocessed_metadata.csv"
    sed -i -e "1!{s/^/${sample_name},/}" "${sample_name}-preprocessed_metadata.csv"
    rm filename.csv label.csv
    """
}

/* Collects all the preprocessed_tiff_metadata_by_sample metadata files:
    - Concatenates them into preprocessed_tiff_metadata.csv
    - Removes extra header lines from the middle of the file, as each of the
        original files starts with an header line.
    - Copies the metadata to "$params.output_folder/preprocessed_tiff_metadata.csv" 
*/

process collect_cp3_preprocessed_metadata {

    label 'small_memory'
    publishDir "$params.output_folder", mode:'copy', overwrite: true
    
    input:
        file metadata_list from preprocessed_tiff_metadata_by_sample.collect()
    
    output:
        file "preprocessed_tiff_metadata.csv" into preprocessed_tiff_metadata

    script:
    """
    cat $metadata_list > preprocessed_tiff_metadata.csv
    sed -i '1!{/sample_name,label,preprocessed_file_name/d;}' preprocessed_tiff_metadata.csv
    """
}
