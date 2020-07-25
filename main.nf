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
params.cell_analysis_metadata = "$baseDir/example_data/cell_analysis_metadata.csv"

params.high_color = "'#FF0000'"
params.mid_color = "'#FFFFFF'"
params.low_color = "'#0000FF'"

image_folder = "$params.output_folder/Images"

/* Gets the Siungulairty key to verify the containers used by the pipeline */

process get_singularity_key {

    label 'small_memory'

    output:
        val true into singularity_key_got
                            
    script:
    """
    singularity keys pull 25892DAEC49C3AAA9691A5CF8661311A5FB2DD90
    """                             
}

/* Reads the raw metadata file line by line to extract the sample metadata for the raw IMC acquisition files.
   It expects an header line and it extracts the following fields into the sample_metadata channel:
   - sample_name
   - roi_name
   - raw_path -> Converted to a file type 
*/

Channel
    .fromPath(params.raw_metadata_file)
    .splitCsv(header:true)
    .map{row -> tuple(row.sample_name, row.roi_name, file(row.raw_path))}
    .set{sample_metadata_raw}

/* For each aquisition specified in the $raw_metadata_file:
    - Extracts the raw tiff files into the working directory
    - Generates the raw tiff .csv metadata file in the working directory
    - Copies the output files into "$image_folder/Raw/$sample_name" 
  For ome.tiff output the channels in the .ome are in the same order they are
  written in the metadata.
*/

process convert_raw_data_to_tiffs {

    label 'big_memory'
    publishDir "$image_folder/Raw/$sample_name", mode:'copy', overwrite: true
    container = 'library://michelebortol/default/simpli_imctools:reupload'
    containerOptions = "--bind $script_folder:/opt"

    input:
        val flag from singularity_key_got
        set sample_name, roi_name, file(raw_path) from sample_metadata_raw

    output:
        file "$sample_name*raw*tiff" into raw_tiff_images
        file "${sample_name}-raw_tiff_metadata.csv" into raw_tiff_metadata_by_sample

    script:
    """
    python3.8 /opt/Tiff_extracter.py \\
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
    - Copies the metadata to "$image_folder/Raw" 
*/

process collect_raw_tiff_metadata {

    label 'small_memory'
    publishDir "$image_folder/Raw/", mode:'copy', overwrite: true
    
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
    - Copies the results to "$image_folder/Normalized/$sample_name"
  For ome.tiff output the channels in the .ome are in the same order they are
  written in the metadata.
*/

process normalize_tiff {
    
    label 'big_memory'
    container = 'library://michelebortol/default/simpli_rbioconductor:ggrepel'
    containerOptions = "--bind $script_folder:/opt"
    
    publishDir "$image_folder/Normalized/$sample_name", mode:'copy', overwrite: true
    
    input:
        set sample_name, file(raw_metadata_file) from normalize_tiff_metadata

    output:
        file "$sample_name*normalized*tiff" into normalized_tiff_images
        file "${sample_name}-normalized_tiff_metadata.csv" into normalized_tiff_metadata_by_sample
        file "${sample_name}-cp3_normalized_tiff_metadata.csv" into cp3_normalized_tiff_metadata_by_sample
    script:
    """
    Rscript /opt/Tiff_normalizer.R \\
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
    - Copies the metadata to "$image_folder/Normalized/" 
*/

process collect_normalized_tiff_metadata {

    label 'small_memory'
    publishDir "$image_folder/Normalized", mode:'copy', overwrite: true
    
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

/* Collects all the cp3_normalized_tiff_metadata_by_sample metadata files and
    pairs them with the corresponding sample name
*/

cp3_normalized_tiff_metadata_by_sample
    .map { file ->
            def key = file.name.toString().tokenize('-').get(0)
            return tuple(key, file)
            }
    .groupTuple()
    .set{ cp3_preprocessing_metadata }

/* Process each sample with CellProfiler3 using:
    - CellProfiler3 pipeline file: $params.cp3_preprocessing_cppipe
    - Image metadata emitted by: cp3_normalized_tiff_metadata_by_sample
    It expects the pipeline to produce output images as single tiffs with this name pattern:
        "SAMPLE-MARKER-Preprocessed.tiff"
    The output images and a metadata file are exported to: "$image_folder/Preprocessed/$sample_name"  
*/

process cell_profiler_image_preprocessing {

    label 'big_memory'
    container = 'library://michelebortol/default/simpli_cp3:reupload'
    containerOptions = "--bind $cp3_pipeline_folder:/mnt"

    publishDir "$image_folder/Preprocessed/$sample_name", mode:'copy', overwrite: true
                                                                                                
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
    - Copies the metadata to "$image_folder/Preprocessed/preprocessed_tiff_metadata.csv" 
*/

process process_cp3_preprocessed_metadata {

    label 'mid_memory'
    container = 'library://michelebortol/default/simpli_rbioconductor:ggrepel'
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
    container = 'library://michelebortol/default/simpli_rbioconductor:ggrepel'
    containerOptions = "--bind $script_folder:/opt"

    publishDir "$params.output_folder", mode:'copy', overwrite: true
    
    input:
        file tiff_metadata from preprocessed_tiff_metadata 
        file area_metadata from measurement_metadata
    
    output:
        file "area_measurements.csv" into area_measurements

    script:
    """
    Rscript /opt/Area_measurer.R \\
        $area_metadata \\
        $tiff_metadata \\
        area_measurements.csv > pixel_area_measurements.log 2>&1
    """
}

/* Reads the $params.raw_metadata_file to extract the names of the columns indicating the categories
    to compare (all columns after the 4th)
*/

Channel
    .fromPath(params.raw_metadata_file)
    .splitCsv(header:false)
    .first()
    .map{row -> row.drop(4).join(",")}
    .into{categories_area; categories_type; categories_cluster}

/* For each category to compare if there are 2 types of samples:
    For each main marker:
        compare all its combinationations with other markers between the two types of samples,
        and make one boxplot each. Output the boxplots in a single multipage pdf file:
            $params.output_folder/Plots/Area_Plots/Boxplots/MAIN_MARKER/MAIN_MARKER_boxplots-CATEGORY.pdf
*/

process area_visualization {

    label 'mid_memory'
    publishDir "$params.output_folder/Plots/Area_Plots/", mode:'copy', overwrite: true
    container = 'library://michelebortol/default/simpli_rbioconductor:ggrepel'
    containerOptions = "--bind $script_folder:/opt"

    input:
        val category_columns from categories_area
        file area_file from area_measurements 
        file sample_metadata_file from file(params.raw_metadata_file)

    output:
        file "*/**/*.pdf" optional true into area_plots
    
    script:
    """
    Rscript /opt/Area_Plotter.R \\
        $area_file \\
        $sample_metadata_file \\
        $category_columns \\
        . > area_plotting_log.txt 2>&1
    """
}

/* Collects all the cp3_preprocessed_tiff_metadata_by_sample metadata files and
    pairs them with the corresponding sample name
*/

cp3_preprocessed_tiff_metadata_by_sample
    .flatten()
    .map { file ->
            def key = file.name.toString().tokenize('-').get(0)
            return tuple(key, file)
            }
    .groupTuple()
    .set{ cp3_segmentation_metadata }

/* Perform cell segmentation and size, shape, position and intensity measurements
    - For each sample the $params.cp3_segmentation_cppipei should produce:
        - A table with single cell measurements: ${sample_name}-Cells.csv
        - A uint16 cell mask: ${sample_name}-Cell_Mask.tiff
    - Copies the results to "$params.output_folder/Segmentation/$sample_name" 
*/

process cell_segmentation {

    label 'big_memory'
    container = 'library://michelebortol/default/simpli_cp3:reupload'
    containerOptions = "--bind $cp3_pipeline_folder:/mnt"

    publishDir "$params.output_folder/Segmentation/$sample_name", mode:'copy', overwrite: true
                                                                                                
    input:
        set sample_name, file(cp3_preprocessed_metadata) from cp3_segmentation_metadata 

    output:
        file "${sample_name}-Cells.csv" into cell_data_csv_by_sample
        file "${sample_name}-Cell_Mask.tiff" into cell_mask_tiffs
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

/* Assign each cell to a main cell type 
    - Outputs: $params.output_folder/annotated_cells.csv
*/

process cell_type_identification {
    label 'mid_memory'
    container = 'library://michelebortol/default/simpli_rbioconductor:ggrepel'
    containerOptions = "--bind $script_folder:/opt"
    publishDir "$params.output_folder", mode:'copy', overwrite: true
    
    input:
        file unannotated_cell_data_file from unannotated_cell_data 
        file cell_threshold_metadata_file from file(params.cell_analysis_metadata)
    
    output:               
        file "annotated_cells.csv" into annotated_cell_data

    script:
    """
    Rscript /opt/Cell_Type_Selecter.R \\
        $unannotated_cell_data_file \\
        $cell_threshold_metadata_file \\
        annotated_cells.csv
    """
}

annotated_cell_data.into{cell_data_to_plot; cell_data_to_cluster}


/* Visualization of cell annotations as main cell types:
   - For each category to compare if there are 2 types of samples:
        - 1 Boxplot for each main cell type, collected in: $params.output_folder/Plots/Cell_Type_Plots/Boxplots/boxplots-CATEGORY.pdf
   - For each category and once more by sample:
        - Barplot with the percentage of each cell type in the category type or sample, collected into
          a single file in: $params.output_folder/Plots/Cell_Type_Plots/Barplots/barplots.pdf
   - For each sample:
        - Overlay of all cells coloured by cell type: $params.output_folder/Plots/Cell_Type_Plots/Overlays/overlay-SAMPLE_NAME.tiff 
   -PDF color legend: $params.output_folder/Plots/Cell_Type_Plots/Overlays/overlay_legend.pdf
   Colors are specified in: $params.cell_analysis_metadata
*/

process cell_type_visualization {

    label 'big_memory'
    publishDir "$params.output_folder/Plots/Cell_Type_Plots", mode:'copy', overwrite: true
    container = 'library://michelebortol/default/simpli_rbioconductor:ggrepel'
    containerOptions = "--bind $script_folder:/opt"

    input:
        val category_columns from categories_type
        file(annotated_cell_file) from cell_data_to_plot
        file sample_metadata_file from file(params.raw_metadata_file)
        file cell_metadata_file from file(params.cell_analysis_metadata)
        file cell_mask_list from cell_mask_tiffs.collect()

    output:
        file "*/*.pdf" optional true into cell_type_plots
        file "*/*.tiff" optional true into cell_type_overlays
    
    script:
    """
    Rscript /opt/Cell_Type_Plotter.R \\
        $annotated_cell_file \\
        $sample_metadata_file \\
        $category_columns \\
        $cell_metadata_file \\
        . \\
        $cell_mask_list > cell_type_plotting.txt 2>&1
    """
}

/* Reads the clustering metadata file line by line to extract:
   It expects an header line and it extracts the following fields into the sample_metadata channel:
   - cell_type
   - markers: "@" separated list of valid measurement names from $params.output_folder/annotated_cells.csv
   - resolutions: "@" separated list of floats
   If any of these fields is "NA" the line is skipped.
*/

Channel
    .fromPath(params.cell_analysis_metadata)
    .splitCsv(header:true)
    .map{row -> tuple(row.cell_type, row.clustering_markers, row.clustering_resolutions)}
    .filter{!it.contains("NA")}
    .combine(cell_data_to_cluster)
    .set{clustering_metadata}

/* For each cell_type (line) not containing "NA"in the cell_type, markers, or resolutions, fields
   in $params.cell_analysis_metadata file:
    - Performs clustering with Seurat with the given markers for the given resolutions 
    - Copies the output files into params.output_folder/Cell_Clusters
*/

process cluster_cells {

    label 'big_memory'
    publishDir "$params.output_folder/Cell_Clusters", mode:'copy', overwrite: true
    container = 'library://michelebortol/default/simpli_rbioconductor:ggrepel'
    containerOptions = "--bind $script_folder:/opt"

    input:
        set cell_type, markers, resolutions, file(annotated_cell_file) from clustering_metadata

    output:
        file "$cell_type/*_clusters.csv" into cluster_csv_files
        file "$cell_type/*_clusters.RData" into cluster_rdata_files
    
    script:
    """
    Rscript /opt/Seurat_Runner.R \\
        $annotated_cell_file \\
        $cell_type \\
        $markers \\
        $resolutions \\
        $cell_type \\
        $cell_type > clustering_log.txt 2>&1
    """
}

/* Collects all the cluster_csv_files single clustered cell data files:
    - Concatenates them into $params.output/Cell_Clusters/clustered_cells.csv
    - Removes extra header lines from the middle of the file, as each of the
        original files starts with an header line.
*/

process collect_clustering_data {

    label 'small_memory'
    publishDir "$params.output_folder/Cell_Clusters", mode:'copy', overwrite: true
    
    input:
        file cluster_list from cluster_csv_files.collect()
    
    output:
        file "clustered_cells.csv" into clustered_cell_data

    script:
    """
    cat $cluster_list > clustered_cells.csv
    sed -i '1!{/.*Metadata_sample_name.*/d;}' clustered_cells.csv
    """
}

/* Visualization of cell annotations as main cell types, for each cell type and clustering resolution:
   - For each category to compare:
        - An heatmap with marker expression by cluster, and a boxplot for each cluster (if there are two groups in the category): 
        $params.output_folder/Plots/Cell_Cluster_Plots/CELL_TYPE/Cluster_Comparisons/CATEGORY-RESOLUTION-plots.pdf
   - UMAPS by cluster, marker, and sample:
        $params.output_folder/Plots/Cell_Cluster_Plots/CELL_TYPE/UMAPs/UMAPs-RESOLUTION.pdf
*/

Channel
    .fromPath(params.cell_analysis_metadata)
    .splitCsv(header:true)
    .map{row -> tuple(row.cell_type, row.clustering_markers, row.clustering_resolutions)}
    .filter{!it.contains("NA")}
    .combine(clustered_cell_data)
    .set{cell_visualization_metadata}

process cell_cluster_visualization {

    label 'mid_memory'
    publishDir "$params.output_folder/Plots/Cell_Cluster_Plots", mode:'copy', overwrite: true
    container = 'library://michelebortol/default/simpli_rbioconductor:ggrepel'
    containerOptions = "--bind $script_folder:/opt"

    input:
        val category_columns from categories_cluster
        set cell_type, markers, resolutions, file(clustered_cell_file) from cell_visualization_metadata
        file sample_metadata_file from file(params.raw_metadata_file)

    output:
        file "$cell_type/**/*.pdf" optional true into cell_plots
    
    script:
    """
    Rscript /opt/Cell_Cluster_Plotter.R \\
        $clustered_cell_file \\
        $sample_metadata_file \\
        $category_columns \\
        $cell_type \\
        $markers \\
        $resolutions \\
        $params.high_color \\
        $params.mid_color \\
        $params.low_color \\
        $cell_type \\
        $cell_type > cell_plotting_log.txt 2>&1
    """
}
