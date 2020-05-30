script_folder = "$baseDir/scripts"
cp3_pipeline_folder = "$baseDir/cell_profiler_pipelines"

// Input parameters: the default values will be moved to a test profile nextflow configuration
params.output_folder ="$baseDir/example_output"
params.data_folder = "$baseDir/example_data"
params.channel_metadata = "$params.data_folder/channel_metadata.csv" 
params.raw_metadata_file = "$baseDir/example_data/sample_metadata.csv" 
params.tiff_type = "single" 
params.cp3_preprocessing_cppipe = "single_preprocessing_example.cppipe"
params.cp3_segmentation_cppipe = "example_segmentation_pipeline.cppipe" 
params.area_measurements_metadata = "$baseDir/example_data/marker_area_metadata.csv"
params.cell_threshold_metadata = "$baseDir/example_data/cell_threshold_metadata.csv"
params.cell_clustering_metadata = "$baseDir/example_data/cell_clustering_metadata.csv"

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
    container = 'library://michelebortol/default/imctools1.0.7:test'
    containerOptions = "--bind $script_folder:/opt"

    input:
        set sample_name, roi_name, file(raw_path), tiff_path from sample_metadata_raw

    output:
        file "$sample_name*raw*tiff" into raw_tiff_images
        file "${sample_name}-raw_tiff_metadata.csv" into raw_tiff_metadata_by_sample
    script:
    """
    python3.8 /opt/tiff_extracter.py \\
        $sample_name \\
        '$roi_name' \\
        $raw_path \\
        $params.tiff_type \\
        ./ \\
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
    .set{ cp3_preprocessing_metadata }

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

/* Processess all the preprocessed_tiff_metadata_by_sample metadata files:
    - Concatenates them into preprocessed_tiff_metadata.csv
    - Removes extra header lines from the middle of the file, as each of the
        original files starts with an header line.
    - Copies the metadata to "$params.output_folder/preprocessed_tiff_metadata.csv" 
*/

process process_cp3_preprocessed_metadata {

    label 'mid_memory'
    container = 'library://michelebortol/default/imcpipeline-rbioconductor:test'
    containerOptions = "--bind $script_folder:/opt"
    
    input:
        file metadata_list from preprocessed_tiff_metadata_by_sample.collect()
    
    output:
        file "preprocessed_tiff_metadata.csv" into preprocessed_tiff_metadata
        file "*-cp3-preprocessed_metadata.csv" into cp3_preprocessed_tiff_metadata_by_sample

    script:
    """
    Rscript /opt/CP3_metadata_maker.R \\
        $params.tiff_type \\
        ./ \\
        $metadata_list
    cat $metadata_list > preprocessed_tiff_metadata.csv
    sed -i '1!{/sample_name,label,preprocessed_file_name/d;}' preprocessed_tiff_metadata.csv
    """
}

/* Measures the user specified areas and their ratios:
    - Measures the areas specified in: $params.area_measurements_metadata
    - Copies the results to "$params.output_folder/area_measurements.csv" 
*/

measurement_metadata = Channel.fromPath(params.area_measurements_metadata)

process pixel_area_measurements {

    label 'big_memory'
    container = 'library://michelebortol/default/imcpipeline-rbioconductor:test'
    containerOptions = "--bind $script_folder:/opt"

    publishDir "$params.output_folder", mode:'copy', overwrite: true
    
    input:
        file tiff_metadata from preprocessed_tiff_metadata 
        file area_metadata from measurement_metadata
    
    output:
        file "area_measurements.csv" into area_measurements

    script:
    """
    Rscript /opt/pixel_measurements.R \\
        $area_metadata \\
        $tiff_metadata \\
        area_measurements.csv > pixel_area_measurements.log 2>&1
    """
}

/* Perform cell segmentation and size, shape, position and intensity measurements
    - For each sample produce:
        - A table with single cell measurements.
        - A uint16 cell mask
    - Copies the results to "$params.output_folder/$ample_name/" 
*/

cp3_preprocessed_tiff_metadata_by_sample
    .flatten()
    .map { file ->
            def key = file.name.toString().tokenize('-').get(0)
            return tuple(key, file)
            }
    .groupTuple()
    .set{ cp3_segmentation_metadata }

process cell_segmentation {

    label 'big_memory'
    container = 'library://michelebortol/default/cellprofiler3_imcplugins:example'
    containerOptions = "--bind $cp3_pipeline_folder:/mnt"

    publishDir "$params.output_folder/$sample_name/cell_data", mode:'copy', overwrite: true
                                                                                                
    input:
        set sample_name, file(cp3_preprocessed_metadata) from cp3_segmentation_metadata 

    output:
        file "${sample_name}-Cells.csv" into cell_data_csv_by_sample
        file "${sample_name}-Cell_Mask.tiff" into cell_mask_tiff
script:
    """
    cellprofiler \\
        --run-headless \\
        --data-file $cp3_preprocessed_metadata \\
        --pipeline /mnt/$params.cp3_segmentation_cppipe \\
        --plugins-directory /opt/CellProfiler/plugins/ \\
        --image-directory /data \\
        --output-directory ./ \\
        --log-level DEBUG \\
        --temporary-directory ./tmp > cp3_segmentation_log.txt 2>&1
    """
}

/* Collects all the cell_data_csv single cell data files:
    - Concatenates them into unannotated_cells.csv
    - Removes extra header lines from the middle of the file, as each of the
        original files starts with an header line.
*/

process collect_single_cell_data {

    label 'small_memory'
    
    input:
        file cell_data_list from cell_data_csv_by_sample.collect()
    
    output:
        file "unannotated_cells.csv" into unannotated_cell_data

    script:
    """
    cat $cell_data_list > unannotated_cells.csv
    sed -i '1!{/ImageNumber,ObjectNumber,.*/d;}' unannotated_cells.csv
    """
}

/* Assign each cell to a main cell population
    - Outputs annotated_cells.csv
*/

cell_threshold_metadata = Channel.fromPath(params.cell_threshold_metadata)
process cell_population_identification {
    label 'mid_memory'
    container = 'library://michelebortol/default/imcpipeline-rbioconductor:test'
    containerOptions = "--bind $script_folder:/opt"
    publishDir "$params.output_folder", mode:'copy', overwrite: true
    
    input:
        file unannotated_cell_data_file from unannotated_cell_data 
        file cell_threshold_metadata_file from cell_threshold_metadata
    
    output:
        file "annotated_cells.csv" into annotated_cell_data

    script:
    """
    Rscript /opt/cell_selecter.R \\
        $unannotated_cell_data_file \\
        $cell_threshold_metadata_file \\
        annotated_cells.csv
    """
}

/* Reads the raw metadata file line by line to extract the sample metadata for the raw IMC acquisition files.
   It expects an header line and it extracts the following fields into the sample_metadata channel:
   - sample_name
   - roi_name
   - raw_path -> Converted to a file type 
   - tiff_path
*/

Channel
    .fromPath(params.cell_clustering_metadata)
    .splitCsv(header:true)
    .map{row -> tuple(row.population_name, row.markers, row.resolutions)}
    .combine(annotated_cell_data)
    .set{clustering_metadata}

/* For each population in the params.cell_clustering_metadata file:
    - Performs clustering with Seurat with the given markers for the given resolutions 
    - Copies the output files into params.output_folder/CellClusters
*/

process cluster_cells {

    label 'big_memory'
    publishDir "$params.output_folder/CellClusters", mode:'copy', overwrite: true
    container = 'library://michelebortol/default/imcpipeline-rbioconductor:test'
    containerOptions = "--bind $script_folder:/opt"

    input:
        set population_name, markers, resolutions, file(annotated_cell_file) from clustering_metadata

    output:
        file "*_clusters.csv" into clusters_csf_files
        file "*_clusters.RData" into clusters_rdata_files
    
    script:
    """
    Rscript /opt/Seurat_Runner.R \\
        $annotated_cell_file \\
        $population_name \\
        $markers \\
        $resolutions \\
        $population_name \\
        . > clustering_log.txt 2>&1
    """
}
