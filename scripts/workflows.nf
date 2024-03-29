nextflow.enable.dsl=2

script_folder = "$baseDir/scripts"

include {get_singularity_key} from "$script_folder/processes.nf"
include {cp4_format_convert} from "$script_folder/processes.nf"

include {convert_raw_data_to_tiffs} from "$script_folder/processes.nf"
include {collect_raw_tiff_metadata} from "$script_folder/processes.nf"

include {normalize_tiffs} from "$script_folder/processes.nf"
include {collect_normalized_tiff_metadata} from "$script_folder/processes.nf"

include {image_preprocessing} from "$script_folder/processes.nf"
include {process_preprocessed_metadata} from "$script_folder/processes.nf"

include {measure_positive_areas} from "$script_folder/processes.nf"
include {area_visualization} from "$script_folder/processes.nf"

include {cell_segmentation} from "$script_folder/processes.nf"
include {collect_single_cell_data} from "$script_folder/processes.nf"

include {cell_type_identification_mask} from "$script_folder/processes.nf"
include {cell_type_visualization} from "$script_folder/processes.nf"

include {cell_clustering} from "$script_folder/processes.nf"
include {collect_clustering_data} from "$script_folder/processes.nf"
include {cell_cluster_visualization} from "$script_folder/processes.nf"

include {threshold_cells} from "$script_folder/processes.nf"
include {cell_threshold_visualization} from "$script_folder/processes.nf"

include {homotypic_interaction_analysis} from "$script_folder/processes.nf"
include {collect_homotypic_interactions} from "$script_folder/processes.nf"
include {homotypic_interaction_visualization} from "$script_folder/processes.nf"

include{get_heterotypic_distances} from "$script_folder/processes.nf"
include{collect_heterotypic_distances} from "$script_folder/processes.nf"
include {heterotypic_interaction_visualization} from "$script_folder/processes.nf"

workflow singularity_key_getter {
    get_singularity_key()
    emit:
        singularity_key_got = get_singularity_key.out.singularity_key_got 
}

workflow convert_metadata_to_cp4{
    take:
        singularity_key_got
        output_suffix
        metadata_to_convert
    main:
        cp4_format_convert(singularity_key_got, output_suffix, metadata_to_convert)
    emit:
        cp4_metadata = cp4_format_convert.out.cp4_metadata
}

workflow convert_raw_data{
    take:
        singularity_key_got
    main:
        raw_file_metadata = channel.fromPath(params.raw_metadata_file)
            .splitCsv(header:true)
            .map{row -> tuple(row.sample_name, row.roi_name, file(row.file_name))}
        convert_raw_data_to_tiffs(singularity_key_got, raw_file_metadata)
        collect_raw_tiff_metadata(convert_raw_data_to_tiffs.out.raw_tiff_metadata_by_sample.collect())
    emit:
        tiff_images = convert_raw_data_to_tiffs.out.raw_tiff_images
        converted_tiff_metadata = collect_raw_tiff_metadata.out.converted_tiff_metadata
}

workflow normalize_images{
    take:
        singularity_key_got
        sample_names
        normalization_metadata_file
    main:
        normalize_tiffs(singularity_key_got, sample_names, normalization_metadata_file)
        collect_normalized_tiff_metadata(normalize_tiffs.out.normalized_tiff_metadata_by_sample.collect())
    emit:
        normalized_tiff_images = normalize_tiffs.out.normalized_tiff_images
        normalized_tiff_metadata = collect_normalized_tiff_metadata.out.normalized_tiff_metadata
        cp4_normalized_tiff_metadata_by_sample = normalize_tiffs.out.cp4_normalized_tiff_metadata_by_sample
}

workflow preprocess_images{
    take:
        singularity_key_got
        preprocessing_metadata_files
    main:
        image_preprocessing(singularity_key_got, preprocessing_metadata_files)
        process_preprocessed_metadata(image_preprocessing.out.preprocessed_tiff_metadata_by_sample.collect())
    emit:
        preprocessed_tiff_images = image_preprocessing.out.preprocessed_tiff_files
        preprocessed_tiff_metadata = process_preprocessed_metadata.out.preprocessed_tiff_metadata
        cp4_preprocessed_tiff_metadata_by_sample = process_preprocessed_metadata.out.cp4_preprocessed_tiff_metadata_by_sample
}

workflow measure_areas{
    take:
        singularity_key_got
        image_metadata_file
        area_metadata
    main:
        measure_positive_areas(singularity_key_got, image_metadata_file, area_metadata)
    emit:
        area_measurements = measure_positive_areas.out.area_measurements
}

workflow segment_cells{
    take:
        singularity_key_got
        segmentation_metadata_files
    main:
        cell_segmentation(singularity_key_got, segmentation_metadata_files)
        collect_single_cell_data(cell_segmentation.out.cell_data_csv_by_sample.collect(),
            cell_segmentation.out.cell_mask_tiffs.collect())
    emit:
        cell_mask_tiffs = cell_segmentation.out.cell_mask_tiffs
        unannotated_cell_data = collect_single_cell_data.out.unannotated_cell_data
        cell_mask_metadata = collect_single_cell_data.out.cell_mask_metadata
}


workflow identify_cell_types_mask{
    take:
        singularity_key_got
        unannotated_cell_data
        cell_type_metadata
        image_metadata
        mask_metadata
    main:
        cell_type_identification_mask(singularity_key_got, unannotated_cell_data, cell_type_metadata,
            image_metadata, mask_metadata)
    emit:
        annotated_cell_data = cell_type_identification_mask.out.annotated_cell_data
}

workflow cluster_cells{
    take:
        singularity_key_got
        annotated_cell_data
        cell_type_metadata
        sample_metadata
    main:
        cell_clustering(singularity_key_got, annotated_cell_data, cell_type_metadata, sample_metadata)
        collect_clustering_data(cell_clustering.out.cluster_csv_files.collect())
    emit:
        cluster_csv_files = cell_clustering.out.cluster_csv_files
        cluster_rdata_files = cell_clustering.out.cluster_rdata_files
        clustered_cell_data = collect_clustering_data.out.clustered_cell_data
}

workflow threshold_expression{
    take:
        singularity_key_got
        annotated_cell_data
        cell_thresholding_metadata
    main:
        threshold_cells(singularity_key_got, annotated_cell_data, cell_thresholding_metadata)
    emit:
        thresholded_cell_data = threshold_cells.out.thresholded_cell_data
}

workflow analyse_homotypic_interactions{
    take:
        singularity_key_got
        homotypic_metadata
    main:
        homotypic_interaction_analysis(singularity_key_got, homotypic_metadata)
        collect_homotypic_interactions(homotypic_interaction_analysis.out.homotypic_clusters.collect())
    emit:
        collected_homotypic_interactions = collect_homotypic_interactions.out.collected_homotypic_interactions
}

workflow calculate_heterotypic_distances{
    take:
        singularity_key_got
        heterotypic_metadata
    main:
        get_heterotypic_distances(singularity_key_got, heterotypic_metadata)
        collect_heterotypic_distances(get_heterotypic_distances.out.heterotypic_distances.collect())
    emit:
        collected_heterotypic_interactions = collect_heterotypic_distances.out.collected_heterotypic_interactions
}

workflow visualize_areas{
    take:
        singularity_key_got
        area_measurement_file 
        sample_metadata_file
    main:
        area_visualization(singularity_key_got, area_measurement_file, sample_metadata_file)
    emit:
        area_plots =  area_visualization.out.area_plots
}

workflow visualize_cell_types{
    take:
        singularity_key_got
        annotated_cell_file
        sample_metadata_file
        cell_metadata_file
        cell_mask_files 
    main:
        cell_type_visualization(singularity_key_got, annotated_cell_file, sample_metadata_file, cell_metadata_file, cell_mask_files)
    emit:
        cell_type_plots = cell_type_visualization.out.cell_type_plots
        cell_type_overlays = cell_type_visualization.out.cell_type_overlays
}

workflow visualize_cell_clusters{
    take:
        singularity_key_got
        cluster_visualization_metadata
        clustered_cell_file
        sample_metadata_file
    main:
        cell_cluster_visualization(singularity_key_got, cluster_visualization_metadata, clustered_cell_file, sample_metadata_file)
    emit:
        cell_cluster_plots = cell_cluster_visualization.out.cell_cluster_plots
}

workflow visualize_cell_thresholds{
    take:
        singularity_key_got
        thresholded_cell_file
        sample_metadata_file
        threshold_metadata_file
        cell_mask_list
    main:
        cell_threshold_visualization(singularity_key_got, thresholded_cell_file, sample_metadata_file, threshold_metadata_file, cell_mask_list)
    emit:
        cell_threshold_plots = cell_threshold_visualization.out.cell_threshold_plots
}

workflow visualize_homotypic_interactions{
    take:
        singularity_key_got
        dbscan_file_name
        metadata_file_name
        cell_mask_list
    main:
        homotypic_interaction_visualization(singularity_key_got, dbscan_file_name, metadata_file_name, cell_mask_list)
    emit:
        homotypic_interaction_plots = homotypic_interaction_visualization.out.homotypic_interaction_plots
}

workflow visualize_heterotypic_interactions{
    take:
        singularity_key_got
        distance_file_name
        metadata_file_name
        sample_file_name
    main:
        heterotypic_interaction_visualization(singularity_key_got, distance_file_name, metadata_file_name, sample_file_name)
    emit:
        heterotypic_interaction_plots = heterotypic_interaction_visualization.out.heterotypic_interaction_plots
}

