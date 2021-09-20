script_folder = "$baseDir/scripts"
image_folder = "$params.output_folder/Images"

if (params.cp4_preprocessing_cppipe){
    cp4_preprocessing_pipeline_folder = file(params.cp4_preprocessing_cppipe).getParent()
    cp4_preprocessing_pipeline = file(params.cp4_preprocessing_cppipe).getName()
}
if (params.cp4_segmentation_cppipe){
    cp4_segmentation_pipeline_folder = file(params.cp4_segmentation_cppipe).getParent()
    cp4_segmentation_pipeline = file(params.cp4_segmentation_cppipe).getName() 
}

process cp4_format_convert {
    
    label 'mid_memory'
    container = 'library://michelebortol/default/simpli_rbioconductor:dbscan'
    containerOptions = "--bind $script_folder:/opt"
    
    input:
        val(output_suffix)
        path(metadata_files_to_convert)

    output:
        path("*$output_suffix", emit: cp4_metadata)
    
    script:
    """
    Rscript --vanilla /opt/Convert_to_cp4_metadata.R \\
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
        
        tuple val(sample_name), val(roi_name), path(raw_path)
		path(channel_metadata_file)

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
        $channel_metadata_file \\
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
    sed -i '1!{/sample_name,roi_name,marker,label,file_name/d;}' raw_tiff_metadata.csv
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
    container = 'library://michelebortol/default/simpli_rbioconductor:dbscan'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"
    
    publishDir "$image_folder/Normalized/$sample_name", mode:'copy', overwrite: true
    
    input:
        
        val(sample_name)
        path(tiff_input_metadata_file)

    output:
        path("$sample_name*normalized*tiff", emit: normalized_tiff_images)
        path("${sample_name}-normalized_tiff_metadata.csv", emit: normalized_tiff_metadata_by_sample)
        path("${sample_name}-cp4_normalized_tiff_metadata.csv", emit: cp4_normalized_tiff_metadata_by_sample)
    
    script:
    """
    Rscript --vanilla /opt/Tiff_normalizer.R \\
        $sample_name \\
        $tiff_input_metadata_file \\
        $params.tiff_type \\
        ./ \\
        ${sample_name}-normalized_tiff_metadata.csv \\
        ${sample_name}-cp4_normalized_tiff_metadata.csv > normalization_log.txt 2>&1
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
    sed -i '1!{/sample_name,marker,label,Frame,URL/d;}' normalized_tiff_metadata.csv
    """
}

/* Process each sample with CellProfiler4 using:
    - CellProfiler4 pipeline file: $params.cp4_preprocessing_cppipe
    - Image metadata emitted by: cp4_normalized_tiff_metadata_by_sample
    It expects the pipeline to produce output images as single tiffs with this name pattern:
        "SAMPLE-MARKER-Preprocessed.tiff"
    The output images and a metadata file are exported to: "$image_folder/Preprocessed/$sample_name"  
*/

process image_preprocessing {

    label 'big_memory'
    container = 'library://michelebortol/default/simpli_cp4_imcplugins:plugins'
    containerOptions = "--bind $cp4_preprocessing_pipeline_folder:/mnt,$workflow.launchDir/:/data"

    publishDir "$image_folder/Preprocessed/$sample_name", mode:'copy', overwrite: true
                                                                                                
    input:
       
       tuple val(sample_name), path(cp4_normalized_metadata) 

    output:
        path("*-Preprocessed.tiff", emit: preprocessed_tiff_files)
        path("${sample_name}-preprocessed_metadata.csv", emit: preprocessed_tiff_metadata_by_sample)

    script:
    """
    cellprofiler \\
        --run-headless \\
        --data-file $cp4_normalized_metadata \\
        --pipeline /mnt/$cp4_preprocessing_pipeline \\
        --plugins-directory /opt/CellProfiler/plugins/ \\
        --image-directory /data \\
        --output-directory ./ \\
        --log-level DEBUG \\
        --temporary-directory ./tmp > cp4_preprocessing_log.txt 2>&1

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
    container = 'library://michelebortol/default/simpli_rbioconductor:dbscan'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    publishDir "$image_folder/Preprocessed/", mode:'copy', overwrite: true
    
    input:
        path(metadata_list)
    
    output:
        path("preprocessed_tiff_metadata.csv", emit: preprocessed_tiff_metadata)
        path("*-cp4-preprocessed_metadata.csv", emit: cp4_preprocessed_tiff_metadata_by_sample)

    script:
    """
    Rscript --vanilla /opt/CP4_metadata_maker.R \\
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
    container = 'library://michelebortol/default/simpli_rbioconductor:dbscan'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    publishDir "$params.output_folder", mode:'copy', overwrite: true
    
    input:
        
        path(tiff_metadata) 
        path(area_metadata)
    
    output:
        path("area_measurements.csv", emit: area_measurements)

    script:
    """
    Rscript --vanilla /opt/Area_measurer.R \\
        $area_metadata \\
        $tiff_metadata \\
        area_measurements.csv > pixel_area_measurements.log 2>&1
    """
}

/* If there are 2 types of samples:
    For each main marker:
        compare all its combinationations with other markers between the two types of samples,
        and make one boxplot each. Output the boxplots in a single multipage pdf file:
            $params.output_folder/Plots/Area_Plots/Boxplots/MAIN_MARKER/MAIN_MARKER_boxplots-CATEGORY.pdf
*/

process area_visualization {

    label 'mid_memory'
    publishDir "$params.output_folder/Plots/Area_Plots/", mode:'copy', overwrite: true
    container = 'library://michelebortol/default/simpli_rbioconductor:dbscan'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    input:
        
        path(area_file) 
        path(sample_metadata_file)

    output:
        path("*/**/*.pdf", emit: area_plots) optional true
    
    script:
    """
    Rscript --vanilla /opt/Area_Plotter.R \\
        $area_file \\
        $sample_metadata_file \\
        . > area_plotting_log.txt 2>&1
    """
}

/* Perform cell segmentation and size, shape, position and intensity measurements
    - For each sample the $params.cp4_segmentation_cppipe should produce:
        - A table with single cell measurements: ${sample_name}-Cells.csv
        - A uint16 cell mask: ${sample_name}-Cell_Mask.tiff
    - Copies the results to "$params.output_folder/Segmentation/$sample_name" 
*/

process cell_segmentation {

    label 'big_memory'
    container = 'library://michelebortol/default/simpli_cp4_imcplugins:plugins'
    containerOptions = "--bind $cp4_segmentation_pipeline_folder:/mnt,$workflow.launchDir/:/data"

    publishDir"$params.output_folder/Segmentation/$sample_name", mode:'copy', overwrite: true
                                                                                                
    input:
        
        tuple val(sample_name), path(cp4_preprocessed_metadata) 

    output:
        path("${sample_name}-Cells.csv", emit: cell_data_csv_by_sample)
        path("${sample_name}-Cell_Mask.tiff", emit: cell_mask_tiffs)

script:
    """
    cellprofiler \\
        --run-headless \\
        --data-file $cp4_preprocessed_metadata \\
        --pipeline /mnt/$cp4_segmentation_pipeline \\
        --plugins-directory /opt/CellProfiler/plugins/ \\
        --image-directory /data \\
        --output-directory ./ \\
        --log-level DEBUG \\
        --temporary-directory ./tmp > cp4_segmentation_log.txt 2>&1
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

process cell_type_identification_mask {
    label 'big_memory'
    container = 'library://michelebortol/default/simpli_rbioconductor:dbscan'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    publishDir "$params.output_folder", mode:'copy', overwrite: true
    
    input:
        
        path(unannotated_cell_data_file) 
        path(cell_threshold_metadata_file)
        path(image_metadata_file)
        path(mask_metadata_file)
    
    output:               
        path("annotated_cells.csv", emit: annotated_cell_data)

    script:
    """
    Rscript --vanilla /opt/Cell_type_selecter_mask.R \\
        $unannotated_cell_data_file \\
        $cell_threshold_metadata_file \\
        $image_metadata_file \\
        $mask_metadata_file \\
        annotated_cells.csv
    """
}

/* Visualization of cell annotations as main cell types:
   - If there are 2 types of samples:
        - 1 Boxplot for each main cell type, collected in: $params.output_folder/Plots/Cell_Type_Plots/Boxplots/boxplots-CATEGORY.pdf
   - Once by category and once more by sample:
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
    container = 'library://michelebortol/default/simpli_rbioconductor:dbscan'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    input:
        
        path(annotated_cell_file)
        path(sample_metadata_file)
        path(cell_metadata_file)
        path(cell_mask_list) 

    output:
        path("*/*.pdf", emit: cell_type_plots)
        path("*/*.tiff", emit: cell_type_overlays)

    script:
    """
    Rscript --vanilla /opt/Cell_Type_Plotter.R \\
        $annotated_cell_file \\
        $sample_metadata_file \\
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
    container = 'library://michelebortol/default/simpli_rbioconductor:dbscan'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    input:
        
        path(annotated_cell_file)
        tuple val(cell_type), val(markers), val(resolutions)
        path(sample_metadata_file)

    output:
        path("$cell_type/*-clusters.csv", emit: cluster_csv_files)
        path("$cell_type/*-clusters.RData", emit: cluster_rdata_files)
    
    script:
    """
    Rscript --vanilla /opt/Seurat_Runner.R \\
        $annotated_cell_file \\
        $sample_metadata_file \\
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

    label 'mid_memory'
    publishDir "$params.output_folder/", mode:'copy', overwrite: true
    container = 'library://michelebortol/default/simpli_rbioconductor:dbscan'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"
    
    input:
        path(cluster_list)
    
    output:
        path("clustered_cells.csv", emit: clustered_cell_data)

    script:
    """
    Rscript --vanilla /opt/Clustered_Cell_Collecter.R $cluster_list \\
        clustered_cells.csv > clustered_cells_collecting_log.txt 2>&1
    """
}

/* For each cell_type (line) not containing "NA" in the cell_type, markers, or resolutions, fields
   in $params.cell_analysis_metadata file:
    - Performs thresholding on the expressions of the given markers 
    - Copies the output files into params.output_folder/Cell_Thresholds
*/

process threshold_cells {
    label 'mid_memory'
    container = 'library://michelebortol/default/simpli_rbioconductor:dbscan'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    publishDir "$params.output_folder", mode:'copy', overwrite: true
    
    input:
        
        path(annotated_cell_data_file) 
        path(cell_threshold_metadata_file)
    
    output:               
        path("thresholded_cells.csv", emit: thresholded_cell_data)

    script:
    """
    Rscript --vanilla /opt/expression_threshold.R \\
        $annotated_cell_data_file \\
        $cell_threshold_metadata_file \\
        thresholded_cells.csv
    """
}

/* Visualization of cell annotations as main cell types, for each cell type and clustering resolution:
    - An heatmap with marker expression by cluster, and a boxplot for each cluster (if there are two groups in the category): 
        $params.output_folder/Plots/Cell_Cluster_Plots/CELL_TYPE/Cluster_Comparisons/CATEGORY-RESOLUTION-plots.pdf
   - UMAPS by cluster, marker, and sample:
        $params.output_folder/Plots/Cell_Cluster_Plots/CELL_TYPE/UMAPs/UMAPs-RESOLUTION.pdf
*/

process cell_cluster_visualization {

    label 'big_memory'
    publishDir "$params.output_folder/Plots/Cell_Cluster_Plots", mode:'copy', overwrite: true
    container = 'library://michelebortol/default/simpli_rbioconductor:dbscan'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    input:
        
        tuple val(cell_type), val(markers), val(resolutions)
        path(clustered_cell_file)
        path(sample_metadata_file)

    output:
        path("$cell_type/**/*.pdf", emit: cell_cluster_plots)
         
    script:
    """
    Rscript --vanilla /opt/Cell_Cluster_Plotter.R \\
        $clustered_cell_file \\
        $sample_metadata_file \\
        $cell_type \\
        $markers \\
        $resolutions \\
        $params.high_color \\
        $params.mid_color \\
        $params.low_color \\
        $cell_type > clustering_plotting_log.txt 2>&1
    """
}

/* Visualization of thresholded cell subpopulations:
    - Heatmaps with marker expression by main cell type and supopulation:
    - Boxplot for each subpopulation (2 categories only): 
    - Barplots for each main cell type, by samples and comparison (2 categories only): 
    - Density plots for each main cell type, subpopulation, thresholded marker:
*/

process cell_threshold_visualization {

    label 'big_memory'
    publishDir "$params.output_folder/Plots/Cell_Threshold_Plots", mode:'copy', overwrite: true
    container = 'library://michelebortol/default/simpli_rbioconductor:dbscan'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    input:
        
        path(thresholded_cell_file)
        path(sample_metadata_file)
        path(threshold_metadata_file)
        path(cell_mask_list) 

    output:
        path("**/*.pdf", emit: cell_threshold_plots) 
        path("**/*.tiff", emit: cell_threshold_overlays)
         
    script:
    """
    Rscript --vanilla /opt/Cell_Threshold_Plotter.R \\
        $thresholded_cell_file \\
        $sample_metadata_file \\
        $threshold_metadata_file \\
        $params.high_color \\
        $params.mid_color \\
        $params.low_color \\
        . \\
        $cell_mask_list > threshold_plotting_log.txt 2>&1
    """
}

process homotypic_interaction_analysis {

    label 'big_memory'
    publishDir "$params.output_folder/Homotypic_interactions", mode:'copy', overwrite: true
    container = 'library://michelebortol/default/simpli_rbioconductor:dbscan'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    input:
        
        tuple path(coordinates_file_name), val(cell_type_column), val(cell_type_to_cluster), val(reachability_distance), val(min_cells)
    output:
        path("$cell_type_to_cluster/$cell_type_to_cluster-homotypic_clusters.csv", emit: homotypic_clusters) 
    script:
    """
    Rscript --vanilla /opt/Homotypic_spatial_analysis.R \\
        $coordinates_file_name \\
        $reachability_distance \\
        $min_cells \\
        $cell_type_column \\
        $cell_type_to_cluster \\
        $cell_type_to_cluster \\
        $cell_type_to_cluster-homotypic_clusters.csv > homotypic_analysis_log.txt 2>&1
    """
}

process collect_homotypic_interactions {

    label 'small_memory'
    publishDir "$params.output_folder/Homotypic_interactions", mode:'copy', overwrite: true
    
    input:
        path(homotypic_interactions_list)

    output:
        path("homotypic_interactions.csv", emit: collected_homotypic_interactions)

    script:
    """
    cat $homotypic_interactions_list > homotypic_interactions.csv
    sed -i '1!{/CellName,Metadata_sample_name,Location_Center_X,Location_Center_Y,spatial_analysis_cell_type,cluster,isseed/d;}' homotypic_interactions.csv
    """
}

process homotypic_interaction_visualization {

    label 'big_memory'
    publishDir "$params.output_folder/Plots/Homotypic_Interaction_Plots", mode:'copy', overwrite: true
    container = 'library://michelebortol/default/simpli_rbioconductor:dbscan'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    input:
        
        path(dbscan_file_name)
        path(metadata_file_name)
        path(cell_mask_list) 

    output:
        path("**/*.pdf", emit: homotypic_interaction_plots) 
    script:
    """
    Rscript --vanilla /opt/Homotypic_spatial_plotting.R \\
        $dbscan_file_name \\
        $metadata_file_name \\
        ./ \\
        $cell_mask_list > homotypic_plotting_log.txt 2>&1
    """
}

process get_heterotypic_distances {

    label 'big_memory'
    publishDir "$params.output_folder/Heterotypic_interactions", mode:'copy', overwrite: true
    container = 'library://michelebortol/default/simpli_rbioconductor:dbscan'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    input:
        
        tuple path(coordinate_file1, stageAs: "cf1.csv"), val(cell_type_column1), val(cell_type1),
            path(coordinate_file2, stageAs: "cf2.csv"), val(cell_type_column2), val(cell_type2)
    output:
        path("$cell_type1-$cell_type2/$cell_type1-$cell_type2-distances.csv", emit: heterotypic_distances) 
    script:
    """
    Rscript --vanilla /opt/Distance_Calculator.R \\
        $coordinate_file1 \\
        $cell_type_column1 \\
        $cell_type1 \\
        $coordinate_file2 \\
        $cell_type_column2 \\
        $cell_type2 \\
        $cell_type1-$cell_type2 \\
        $cell_type1-$cell_type2-distances.csv > heterotypic_distances_log.txt 2>&1
    """
}

process collect_heterotypic_distances {
    label 'small_memory'
    publishDir "$params.output_folder/Heterotypic_interactions", mode:'copy', overwrite: true
    
    input:
        path(heterotypic_interactions_list)

    output:
        path("heterotypic_interactions.csv", emit: collected_heterotypic_interactions)

    script:
    """
    cat $heterotypic_interactions_list > heterotypic_interactions.csv
    sed -i '1!{/Metadata_sample_name,CellName1,Location_Center_X1,Location_Center_Y1,spatial_analysis_cell_type1,CellName2,Location_Center_X2,Location_Center_Y2,spatial_analysis_cell_type2,distance/d;}' heterotypic_interactions.csv
    """
}

process heterotypic_interaction_visualization {

    label 'big_memory'
    publishDir "$params.output_folder/Plots/Heterotypic_Interaction_Plots", mode:'copy', overwrite: true
    container = 'library://michelebortol/default/simpli_rbioconductor:dbscan'
    containerOptions = "--bind $script_folder:/opt,$workflow.launchDir/:/data"

    input:
        path(distance_file_name)
        path(metadata_file_name)
        path(sample_file_name)

    output:
        path("**/*.pdf", emit: heterotypic_interaction_plots) 
    script:
    """
    Rscript --vanilla /opt/Heterotypic_spatial_plotting.R \\
        $distance_file_name \\
        $metadata_file_name \\
        $sample_file_name \\
        ./ > heterotypic_plotting_log.txt 2>&1
    """
}
