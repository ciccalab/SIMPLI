nextflow.enable.dsl=2

script_folder = "$baseDir/scripts"

include {singularity_key_getter} from "$script_folder/workflows.nf"
include {convert_metadata_to_cp4} from "$script_folder/workflows.nf"

include {convert_raw_data} from "$script_folder/workflows.nf"

include {normalize_images} from "$script_folder/workflows.nf"

include {preprocess_images} from "$script_folder/workflows.nf"

include {measure_areas} from "$script_folder/workflows.nf"

include {segment_cells} from "$script_folder/workflows.nf"

include {identify_cell_types_mask} from "$script_folder/workflows.nf"

include {cluster_cells} from "$script_folder/workflows.nf"
include {threshold_expression} from "$script_folder/workflows.nf"

include {analyse_homotypic_interactions} from "$script_folder/workflows.nf"

include {calculate_heterotypic_distances} from "$script_folder/workflows.nf"

include {visualize_areas} from "$script_folder/workflows.nf"
include {visualize_cell_types} from "$script_folder/workflows.nf"
include {visualize_cell_clusters} from "$script_folder/workflows.nf"
include {visualize_cell_thresholds} from "$script_folder/workflows.nf"
include {visualize_homotypic_interactions} from "$script_folder/workflows.nf"
include {visualize_heterotypic_interactions} from "$script_folder/workflows.nf"


workflow {
    sample_names = channel.fromPath(params.sample_metadata_file).splitCsv(header: true).map{row -> row.sample_name}
    singularity_key_getter()
    if(params.raw_metadata_file && !params.skip_conversion){
        convert_raw_data(singularity_key_getter.out.singularity_key_got)
    }    
    if((params.converted_metadata_file || !params.skip_conversion) && !params.skip_normalization){
        normalization_metadata = (params.skip_conversion) ? channel.fromPath(params.converted_metadata_file) : convert_raw_data.out.converted_tiff_metadata
        normalize_images(singularity_key_getter.out.singularity_key_got, sample_names, normalization_metadata)
    }    
    if(!params.skip_preprocessing){
        if(!params.skip_normalization){
            preprocessing_metadata = normalize_images.out.cp4_normalized_tiff_metadata_by_sample
        }
        else if(params.skip_normalization && !params.skip_conversion){
            convert_metadata_to_cp4(singularity_key_getter.out.singularity_key_got, "-cp4_metadata.csv", convert_raw_data.out.converted_tiff_metadata.collect()) 
            preprocessing_metadata = convert_metadata_to_cp4.out.cp4_metadata.flatten()
        }    
        else{
            metadata_to_convert = channel.fromPath(params.normalized_metadata_file)
            convert_metadata_to_cp4(singularity_key_getter.out.singularity_key_got, "-cp4_metadata.csv", metadata_to_convert) 
            preprocessing_metadata = convert_metadata_to_cp4.out.cp4_metadata.flatten()
        }
        preprocessing_metadata = preprocessing_metadata
            .map{ file ->
                def key = file.name.toString().tokenize('-').get(0)
                return tuple(key, file)
            }
            .groupTuple()
        preprocess_images(singularity_key_getter.out.singularity_key_got, preprocessing_metadata) 
    }    
    if(!params.skip_area){
        if(!params.skip_preprocessing){
            area_measurement_metadata = preprocess_images.out.preprocessed_tiff_metadata
        }    
        else{
            area_measurement_metadata = params.preprocessed_metadata_file
        }
        measure_areas(singularity_key_getter.out.singularity_key_got, area_measurement_metadata, params.area_measurements_metadata)
    }    
    if(!params.skip_segmentation){
        if(!params.skip_preprocessing)
            segmentation_metadata = preprocess_images.out.cp4_preprocessed_tiff_metadata_by_sample
        else{
            convert_metadata_to_cp4(singularity_key_getter.out.singularity_key_got, "-cp4_metadata.csv", params.preprocessed_metadata_file) 
            segmentation_metadata = convert_metadata_to_cp4.out.cp4_metadata
        }
        segmentation_metadata = segmentation_metadata.flatten()
            .map { file ->
                def key = file.name.toString().tokenize('-').get(0)
                return tuple(key, file)
             }
            .groupTuple()
        segment_cells(singularity_key_getter.out.singularity_key_got, segmentation_metadata)
    }    
    if(!params.skip_cell_type_identification){
        if(!params.skip_segmentation){
            unannotated_cells = segment_cells.out.unannotated_cell_data
            cell_mask_metadata = segment_cells.out.cell_mask_metadata
        }    
        else{
            unannotated_cells = params.single_cell_data_file
            cell_mask_metadata = params.single_cell_masks_metadata
        }
        if(!params.skip_preprocessing)
            image_metadata = preprocess_images.out.preprocessed_tiff_metadata
        else{
            image_metadata = params.preprocessed_metadata_file
        }
        identify_cell_types_mask(singularity_key_getter.out.singularity_key_got, unannotated_cells, params.cell_masking_metadata,
            image_metadata, cell_mask_metadata)
        annotated_cell_data = identify_cell_types_mask.out.annotated_cell_data
    }    
    
    if(!params.skip_cell_thresholding){
        if(!params.skip_cell_type_identification)
            annotated_cells = annotated_cell_data 
        else{
            annotated_cells = params.annotated_cell_data_file
        }
        threshold_expression(singularity_key_getter.out.singularity_key_got, annotated_cells, params.cell_thresholding_metadata)
    }    

    if(!params.skip_cell_clustering){
        if(!params.skip_cell_type_identification)
            annotated_cells = annotated_cell_data 
        else{
            annotated_cells = params.annotated_cell_data_file
        }
        clustering_metadata = channel.fromPath(params.cell_clustering_metadata)
            .splitCsv(header:true)
            .map{row -> tuple(row.cell_type, row.clustering_markers, row.clustering_resolutions)}
            .filter{!it.contains("NA")}
        cluster_cells(singularity_key_getter.out.singularity_key_got, annotated_cells, clustering_metadata, params.sample_metadata_file)
    }

    coord_selecter_map = ["identification": identify_cell_types_mask.out.annotated_cell_data,
        "clustering": cluster_cells.out.clustered_cell_data,
        "thresholding": threshold_expression.out.thresholded_cell_data].withDefault{ key -> key}
    
    if(!params.skip_homotypic_interactions){
        homotypic_metadata = channel.fromPath(params.homotypic_interactions_metadata)
            .splitCsv(header:true)
            .map{row -> tuple(row.cell_type_column, row.cell_type_to_cluster, row.reachability_distance, row.min_cells)}
            .filter{!it.contains("NA")}
     homotypic_interactions_input = coord_selecter_map[params.homotypic_interactions_input]
        analyse_homotypic_interactions(singularity_key_getter.out.singularity_key_got,
            homotypic_interactions_input, homotypic_metadata)
    }
    
    if(!params.skip_heterotypic_interactions){
        heterotypic_metadata = channel.fromPath(params.heterotypic_interactions_metadata)
            .splitCsv(header:true)
            .map{row -> tuple(coord_selecter_map[row.cell_file1].get(), row.cell_type_column1, row.cell_type1,
                    coord_selecter_map[row.cell_file2].get(), row.cell_type_column2, row.cell_type2)}
        calculate_heterotypic_distances(singularity_key_getter.out.singularity_key_got, heterotypic_metadata)
    }
    
    if(!params.skip_visualization){
        if(!params.skip_area && !params.skip_area_visualization){
            visualize_areas(singularity_key_getter.out.singularity_key_got, measure_areas.out.area_measurements,
                params.sample_metadata_file)
        }
        if(!params.skip_type_visualization){
            if(!params.skip_segmentation){
                cell_masks = segment_cells.out.cell_mask_tiffs.collect()
            }
            if(params.skip_segmentation){
                cell_masks = channel.fromPath(params.single_cell_masks_metadata)
                    .splitCsv(header:true)
                    .map{row -> row.file_name}
                    .collect() 
            }
            if(!params.skip_cell_type_identification){
                cell_types = annotated_cell_data 
            }
            if(params.skip_cell_type_identification){
                cell_types = params.annotated_cell_data_file
            }
            visualize_cell_types(singularity_key_getter.out.singularity_key_got, cell_types,
                params.sample_metadata_file, params.cell_masking_metadata, cell_masks)
        }
        if(!params.skip_cluster_visualization){
            if(!params.skip_cell_clustering){
                clustered_cell_file = cluster_cells.out.clustered_cell_data
            }
            else{
                clustered_cell_file = params.clustered_cell_data_file
            }
            clustering_metadata = channel.fromPath(params.cell_clustering_metadata)
                .splitCsv(header:true)
                .map{row -> tuple(row.cell_type, row.clustering_markers, row.clustering_resolutions)}
                .filter{!it.contains("NA")}
            visualize_cell_clusters(singularity_key_getter.out.singularity_key_got, clustering_metadata,
                clustered_cell_file, params.sample_metadata_file)
        }
        if(!params.skip_thresholding_visualization){
            if(!params.skip_cell_thresholding){
                thresholded_cell_file = threshold_expression.out.thresholded_cell_data
            }
            else{
                thresholded_cell_file = params.thresholded_cell_data_file
            }
            if(!params.skip_segmentation){
                cell_masks = segment_cells.out.cell_mask_tiffs.collect()
            }
            else{
                cell_masks = channel.fromPath(params.single_cell_masks_metadata)
                    .splitCsv(header:true)
                    .map{row -> row.file_name}
                    .collect() 
            }
            visualize_cell_thresholds(singularity_key_getter.out.singularity_key_got,
                thresholded_cell_file, params.sample_metadata_file, params.cell_thresholding_metadata, cell_masks)
        }
        if(!params.skip_homotypic_visualization){
            if(!params.skip_homotypic_interactions){
                homotypic_interactions_file = analyse_homotypic_interactions.out.collected_homotypic_interactions
            }
            else{
                homotypic_interactions_file = params.homotypic_interactions_file
            }
            if(!params.skip_segmentation){
                cell_masks = segment_cells.out.cell_mask_tiffs.collect()
            }
            else{
                cell_masks = channel.fromPath(params.single_cell_masks_metadata)
                    .splitCsv(header:true)
                    .map{row -> row.file_name}
                    .collect() 
            }
            visualize_homotypic_interactions(singularity_key_getter.out.singularity_key_got,
                homotypic_interactions_file, params.homotypic_interactions_metadata, cell_masks)
        }
        if(!params.skip_heterotypic_visualization){
            if(!params.skip_heterotypic_interactions){
                heterotypic_interactions_file = calculate_heterotypic_distances.out.collected_heterotypic_interactions
            }
            else{
                heterotypic_interactions_file = params.heterotypic_interactions_file
            }
        }
        visualize_heterotypic_interactions(singularity_key_getter.out.singularity_key_got,
            heterotypic_interactions_file, params.heterotypic_interactions_metadata, params.sample_metadata_file)
    }
}
