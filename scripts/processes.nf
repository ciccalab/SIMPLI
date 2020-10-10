script_folder = "$baseDir/scripts"
image_folder = "$params.output_folder/Images"

if (params.cp3_preprocessing_cppipe){
    cp3_preprocessing_pipeline_folder = file(params.cp3_preprocessing_cppipe).getParent()
    cp3_preprocessing_pipeline = file(params.cp3_preprocessing_cppipe).getName()
}
if (params.cp3_segmentation_cppipe){
    cp3_segmentation_pipeline_folder = file(params.cp3_segmentation_cppipe).getParent()
    cp3_segmentation_pipeline = file(params.cp3_segmentation_cppipe).getName() 
}
/* Gets the Singularity key to verify the containers used by the pipeline */

process get_singularity_key {

    label 'small_memory'

    output:
        val(true, emit: singularity_key_got)
                            
    script:
    """
    singularity keys pull 25892DAEC49C3AAA9691A5CF8661311A5FB2DD90
    """                             
}

process cp3_format_convert {
    
    label 'mid_memory'
    container = 'library://michelebortol/default/simpli_rbioconductor:ggrepel'
    containerOptions = "--bind $script_folder:/opt"
    
    input:
        val(singularity_key_got)
        val(output_suffix)
        path(metadata_files_to_convert)

    output:
        path("*$output_suffix", emit: cp3_metadata)
    
    script:
    """
    Rscript /opt/Convert_to_cp3_metadata.R \\
        $params.tiff_type \\
        ./ \\
        $output_suffix \\
        $metadata_files_to_convert > conversion_log.txt 2>&1
    """
}

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
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    input:
        val(singularity_key_got)
        tuple val(sample_name), val(roi_name), path(raw_path)

    output:
        path("$sample_name*raw*tiff", emit: raw_tiff_images)
        path("${sample_name}-raw_tiff_metadata.csv", emit: raw_tiff_metadata_by_sample)

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
        path(metadata_list)

    output:
        path("raw_tiff_metadata.csv", emit: converted_tiff_metadata)

    script:
    """
    cat $metadata_list > raw_tiff_metadata.csv
    sed -i '1!{/sample_name,roi_name,metal,label,file_name/d;}' raw_tiff_metadata.csv
    """
}

/* For each sample:
    - Performs 99th percentile normalization and scaling
    - Copies the results to "$image_folder/Normalized/$sample_name"
  For ome.tiff output the channels in the .ome are in the same order they are
  written in the metadata.
*/

process normalize_tiffs {
    
    label 'big_memory'
    container = 'library://michelebortol/default/simpli_rbioconductor:ggrepel'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"
    
    publishDir "$image_folder/Normalized/$sample_name", mode:'copy', overwrite: true
    
    input:
        val(singularity_key_got)
        val(sample_name)
        path(converted_metadata_file)

    output:
        path("$sample_name*normalized*tiff", emit: normalized_tiff_images)
        path("${sample_name}-normalized_tiff_metadata.csv", emit: normalized_tiff_metadata_by_sample)
        path("${sample_name}-cp3_normalized_tiff_metadata.csv", emit: cp3_normalized_tiff_metadata_by_sample)
    
    script:
    """
    Rscript /opt/Tiff_normalizer.R \\
        $sample_name \\
        $converted_metadata_file \\
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
        path(metadata_list)
    
    output:
        path("normalized_tiff_metadata.csv", emit: normalized_tiff_metadata)

    script:
    """
    cat $metadata_list > normalized_tiff_metadata.csv
    sed -i '1!{/sample_name,metal,label,Frame,URL/d;}' normalized_tiff_metadata.csv
    """
}

/* Process each sample with CellProfiler3 using:
    - CellProfiler3 pipeline file: $params.cp3_preprocessing_cppipe
    - Image metadata emitted by: cp3_normalized_tiff_metadata_by_sample
    It expects the pipeline to produce output images as single tiffs with this name pattern:
        "SAMPLE-MARKER-Preprocessed.tiff"
    The output images and a metadata file are exported to: "$image_folder/Preprocessed/$sample_name"  
*/

process image_preprocessing {

    label 'big_memory'
    container = 'library://michelebortol/default/simpli_cp3:cp3_fix_dependencies'
    containerOptions = "--bind $cp3_preprocessing_pipeline_folder:/mnt,$workflow.launchDir/:/data"

    publishDir "$image_folder/Preprocessed/$sample_name", mode:'copy', overwrite: true
                                                                                                
    input:
       val(singularity_key_got)
       tuple val(sample_name), path(cp3_normalized_metadata) 

    output:
        path("*-Preprocessed.tiff", emit: preprocessed_tiff_files)
        path("${sample_name}-preprocessed_metadata.csv", emit: preprocessed_tiff_metadata_by_sample)

    script:
    """
    cellprofiler \\
        --run-headless \\
        --data-file $cp3_normalized_metadata \\
        --pipeline /mnt/$cp3_preprocessing_pipeline \\
        --plugins-directory /opt/CellProfiler/plugins/ \\
        --image-directory /data \\
        --output-directory ./ \\
        --log-level DEBUG \\
        --temporary-directory ./tmp > cp3_preprocessing_log.txt 2>&1

    echo "sample_name,label,file_name" > "${sample_name}-preprocessed_metadata.csv"
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

process process_preprocessed_metadata {

    label 'mid_memory'
    container = 'library://michelebortol/default/simpli_rbioconductor:ggrepel'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    publishDir "$image_folder/Preprocessed/", mode:'copy', overwrite: true
    
    input:
        path(metadata_list)
    
    output:
        path("preprocessed_tiff_metadata.csv", emit: preprocessed_tiff_metadata)
        path("*-cp3-preprocessed_metadata.csv", emit: cp3_preprocessed_tiff_metadata_by_sample)

    script:
    """
    Rscript /opt/CP3_metadata_maker.R \\
        $params.tiff_type \\
        ./ \\
        $metadata_list
    cat $metadata_list > preprocessed_tiff_metadata.csv
    sed -i '1!{/sample_name,label,file_name/d;}' preprocessed_tiff_metadata.csv
    """
}

/* Measures the user specified areas and their ratios:
    - Measures the areas specified in: $params.area_measurements_metadata
    - Copies the results to "$params.output_folder/area_measurements.csv" 
*/

process measure_positive_areas {

    label 'big_memory'
    container = 'library://michelebortol/default/simpli_rbioconductor:ggrepel'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    publishDir "$params.output_folder", mode:'copy', overwrite: true
    
    input:
        val(singularity_key_got)
        path(tiff_metadata) 
        path(area_metadata)
    
    output:
        path("area_measurements.csv", emit: area_measurements)

    script:
    """
    Rscript /opt/Area_measurer.R \\
        $area_metadata \\
        $tiff_metadata \\
        area_measurements.csv > pixel_area_measurements.log 2>&1
    """
}

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
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    input:
        val(singularity_key_got)
        val(category_columns)
        path(area_file) 
        path(sample_metadata_file)

    output:
        path("*/**/*.pdf", emit: area_plots) optional true
    
    script:
    """
    Rscript /opt/Area_Plotter.R \\
        $area_file \\
        $sample_metadata_file \\
        $category_columns \\
        . > area_plotting_log.txt 2>&1
    """
}

/* Perform cell segmentation and size, shape, position and intensity measurements
    - For each sample the $params.cp3_segmentation_cppipei should produce:
        - A table with single cell measurements: ${sample_name}-Cells.csv
        - A uint16 cell mask: ${sample_name}-Cell_Mask.tiff
    - Copies the results to "$params.output_folder/Segmentation/$sample_name" 
*/

process cell_segmentation {

    label 'big_memory'
    container = 'library://michelebortol/default/simpli_cp3:cp3_fix_dependencies'
    containerOptions = "--bind $cp3_segmentation_pipeline_folder:/mnt,$workflow.launchDir/:/data"

    publishDir"$params.output_folder/Segmentation/$sample_name", mode:'copy', overwrite: true
                                                                                                
    input:
        val(singularity_key_got)
        tuple val(sample_name), path(cp3_preprocessed_metadata) 

    output:
        path("${sample_name}-Cells.csv", emit: cell_data_csv_by_sample)
        path("${sample_name}-Cell_Mask.tiff", emit: cell_mask_tiffs)

script:
    """
    cellprofiler \\
        --run-headless \\
        --data-file $cp3_preprocessed_metadata \\
        --pipeline /mnt/$cp3_segmentation_pipeline \\
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
    publishDir"$params.output_folder/Segmentation", mode:'copy', overwrite: true
    
    input:
        path(cell_data_list)
        path(cell_mask_list)
    
    output:
        path("unannotated_cells.csv", emit: unannotated_cell_data)
        path("cell_mask_metadata.csv", emit: cell_mask_metadata)

    script:
    """
    cat $cell_data_list > unannotated_cells.csv
    sed -i '1!{/ImageNumber,ObjectNumber,.*/d;}' unannotated_cells.csv

    echo "sample_name,label,file_name" > cell_mask_metadata.csv
    readlink -e $cell_mask_list > filename.csv
    sed "s@.*/\\(.*\\)-\\Cell_Mask\\.tiff@\\1@" filename.csv > sample.csv
    sed "s@.*-\\(Cell_Mask\\).tiff@\\1@" filename.csv > label.csv
    paste -d , sample.csv label.csv filename.csv >> cell_mask_metadata.csv
    rm filename.csv label.csv
    """
}

/* Assign each cell to a main cell type 
    - Outputs: $params.output_folder/annotated_cells.csv
*/

process cell_type_identification_expression {
    label 'mid_memory'
    container = 'library://michelebortol/default/simpli_rbioconductor:ggrepel'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    publishDir "$params.output_folder", mode:'copy', overwrite: true
    
    input:
        val(singularity_key_got)
        path(unannotated_cell_data_file) 
        path(cell_threshold_metadata_file)
    
    output:               
        path("annotated_cells.csv", emit: annotated_cell_data)

    script:
    """
    Rscript /opt/Cell_type_selecter_expression.R \\
        $unannotated_cell_data_file \\
        $cell_threshold_metadata_file \\
        annotated_cells.csv
    """
}

process cell_type_identification_mask {
    label 'big_memory'
    container = 'library://michelebortol/default/simpli_rbioconductor:ggrepel'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    publishDir "$params.output_folder", mode:'copy', overwrite: true
    
    input:
        val(singularity_key_got)
        path(unannotated_cell_data_file) 
        path(cell_threshold_metadata_file)
        path(image_metadata_file)
        path(mask_metadata_file)
    
    output:               
        path("annotated_cells.csv", emit: annotated_cell_data)

    script:
    """
    Rscript /opt/Cell_type_selecter_mask.R \\
        $unannotated_cell_data_file \\
        $cell_threshold_metadata_file \\
        $image_metadata_file \\
        $mask_metadata_file \\
        annotated_cells.csv
    """
}

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
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    input:
        val(singularity_key_got)
        val(category_columns)
        path(annotated_cell_file)
        path(sample_metadata_file)
        path(cell_metadata_file)
        path(cell_mask_list) 

    output:
        path("*/*.pdf", emit: cell_type_plots)  optional true 
        path("*/*.tiff", emit: cell_type_overlays) optional true

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

/* For each cell_type (line) not containing "NA" in the cell_type, markers, or resolutions, fields
   in $params.cell_analysis_metadata file:
    - Performs clustering with Seurat with the given markers for the given resolutions 
    - Copies the output files into params.output_folder/Cell_Clusters
*/

process cell_clustering {

    label 'huge_memory'
    publishDir "$params.output_folder/Cell_Clusters", mode:'copy', overwrite: true
    container = 'library://michelebortol/default/simpli_rbioconductor:ggrepel'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    input:
        val(singularity_key_got)
        path(annotated_cell_file)
        tuple val(cell_type), val(markers), val(resolutions)
        val(category_column)
        path(sample_metadata_file)

    output:
        path("$cell_type-$category_column/*-clusters.csv", emit: cluster_csv_files)
        path("$cell_type-$category_column/*-clusters.RData", emit: cluster_rdata_files)
    
    script:
    """
    Rscript /opt/Seurat_Runner.R \\
        $annotated_cell_file \\
        $sample_metadata_file \\
        $category_column \\
        $cell_type \\
        $markers \\
        $resolutions \\
        "$cell_type-$category_column" \\
        "$cell_type-$category_column" > clustering_log.txt 2>&1
    """
}

/* Collects all the cluster_csv_files single clustered cell data files:
    - Concatenates them into $params.output/Cell_Clusters/clustered_cells.csv
    - Removes extra header lines from the middle of the file, as each of the
        original files starts with an header line.
*/

process collect_clustering_data {

    label 'mid_memory'
    publishDir "$params.output_folder/", mode:'copy', overwrite: true
    container = 'library://michelebortol/default/simpli_rbioconductor:ggrepel'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"
    
    input:
        path(cluster_list)
    
    output:
        path("clustered_cells.csv", emit: clustered_cell_data)

    script:
    """
    Rscript /opt/Clustered_Cell_Collecter.R $cluster_list \\
        clustered_cells.csv > clustered_cells_collecting_log.txt 2>&1
    """
}

/* Visualization of cell annotations as main cell types, for each cell type and clustering resolution:
   - For each category to compare:
        - An heatmap with marker expression by cluster, and a boxplot for each cluster (if there are two groups in the category): 
        $params.output_folder/Plots/Cell_Cluster_Plots/CELL_TYPE/Cluster_Comparisons/CATEGORY-RESOLUTION-plots.pdf
   - UMAPS by cluster, marker, and sample:
        $params.output_folder/Plots/Cell_Cluster_Plots/CELL_TYPE/UMAPs/UMAPs-RESOLUTION.pdf
*/

process cell_cluster_visualization {

    label 'big_memory'
    publishDir "$params.output_folder/Plots/Cell_Cluster_Plots", mode:'copy', overwrite: true
    container = 'library://michelebortol/default/simpli_rbioconductor:ggrepel'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    input:
        val(singularity_key_got)
        val(category_column)
        tuple val(cell_type), val(markers), val(resolutions)
        path(clustered_cell_file)
        path(sample_metadata_file)

    output:
        path("$cell_type-$category_column/**/*.pdf", emit: cell_cluster_plots) optional true
         
    script:
    """
    Rscript /opt/Cell_Cluster_Plotter.R \\
        $clustered_cell_file \\
        $sample_metadata_file \\
        $category_column \\
        $cell_type \\
        $markers \\
        $resolutions \\
        $params.high_color \\
        $params.mid_color \\
        $params.low_color \\
        $cell_type-$category_column > cell_plotting_log.txt 2>&1
    """
}
