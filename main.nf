nextflow.enable.dsl=2

script_folder = "$baseDir/scripts"

include {get_singularity_key} from "$script_folder/processes.nf"

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

include {cell_type_identification} from "$script_folder/processes.nf"
include {cell_type_visualization} from "$script_folder/processes.nf"

include {cell_clustering} from "$script_folder/processes.nf"
include {collect_clustering_data} from "$script_folder/processes.nf"
include {cell_cluster_visualization} from "$script_folder/processes.nf"

workflow singularity_key_getter {
    get_singularity_key()
    emit:
        singularity_key_got = get_singularity_key.out.singularity_key_got 
}

workflow convert_raw_data{
    take:
        singularity_key_got
    main:
        raw_file_metadata = channel.fromPath(params.raw_metadata_file)
            .splitCsv(header:true)
            .map{row -> tuple(row.sample_name, row.roi_name, file(row.raw_path))}
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
        cp3_normalized_tiff_metadata_by_sample = normalize_tiffs.out.cp3_normalized_tiff_metadata_by_sample
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
        cp3_preprocessed_tiff_metadata_by_sample = process_preprocessed_metadata.out.cp3_preprocessed_tiff_metadata_by_sample
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
        collect_single_cell_data(cell_segmentation.out.cell_data_csv_by_sample.collect())
    emit:
        cell_mask_tiffs = cell_segmentation.out.cell_mask_tiffs
        unannotated_cell_data = collect_single_cell_data.out.unannotated_cell_data
}

workflow identify_cell_types{
    take:
        singularity_key_got
        unannotated_cell_data
        cell_type_metadata
    main:
        cell_type_identification(singularity_key_got, unannotated_cell_data, cell_type_metadata)
    emit:
        annotated_cell_data = cell_type_identification.out.annotated_cell_data
}

workflow cluster_cells{
    take:
        singularity_key_got
        annotated_cell_data
        cell_type_metadata
    main:
        cell_clustering(singularity_key_got, annotated_cell_data, cell_type_metadata)
        collect_clustering_data(cell_clustering.out.cluster_csv_files)
    emit:
        cluster_csv_files = cell_clustering.out.cluster_csv_files
        cluster_rdata_files = cell_clustering.out.cluster_rdata_files
        clustered_cell_data = collect_clustering_data.out.clustered_cell_data
}

workflow visualize_areas{
    take:
        singularity_key_got
        categories
        area_measurement_file 
        sample_metadata_file
    main:
        area_visualization(singularity_key_got, categories, area_measurement_file, sample_metadata_file)
    emit:
        area_plots =  area_visualization.out.area_plots
}

workflow visualize_cell_types{
    take:
        singularity_key_got
        categories
        annotated_cell_file
        sample_metadata_file
        cell_metadata_file
        cell_mask_files 
    main:
        cell_type_visualization(singularity_key_got, categories, annotated_cell_file, sample_metadata_file, cell_metadata_file, cell_mask_files)
    emit:
        cell_type_plots = cell_type_visualization.out.cell_type_plots
        cell_type_overlays = cell_type_visualization.out.cell_type_overlays
}

workflow visualize_cell_clusters{
    take:
        singularity_key_got
        categories
        cluster_visualization_metadata
        clustered_cell_file
        sample_metadata_file
    main:
        cell_cluster_visualization(singularity_key_got, categories, cluster_visualization_metadata, clustered_cell_file, sample_metadata_file)
    emit:
        cell_cluster_plots = cell_cluster_visualization.out.cell_cluster_plots
}

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
            preprocessing_metadata = normalize_images.out.cp3_normalized_tiff_metadata_by_sample
        }
        else if(params.skip_normalization && !params.skip_conversion){}
            // preprocessing_metadata = convert_raw_data.out.cp3_converted_tiff_metadata_by_sample *** NOT IMPLEMENTED YET! ***
        else{}
            // metadata_to_convert = channel.fromPath(params.preprocessing_image_metadata_file)
            // convert_to_cp3_format(metadata_to_convert)
            // preprocessing_metadata = convert_to_cp3_format.out.cp3_preprocessing_image_metadata *** NOT IMPLEMENTED YET! ***
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
        else{}// *** NOT IMPLEMENTED YET ***
        measure_areas(singularity_key_getter.out.singularity_key_got, area_measurement_metadata, params.area_measurements_metadata)
    }    
    if(!params.skip_segmentation){
        if(!params.skip_preprocessing)
            segmentation_metadata = preprocess_images.out.cp3_preprocessed_tiff_metadata_by_sample
                .flatten()
                .map { file ->
                    def key = file.name.toString().tokenize('-').get(0)
                    return tuple(key, file)
                    }
                .groupTuple()
        else{}// *** NOT IMPLEMENTED YET ***
        segment_cells(singularity_key_getter.out.singularity_key_got, segmentation_metadata)
    }    
    if(!params.skip_cell_type_identification){
        if(!params.skip_segmentation)
            unannotated_cells = segment_cells.out.unannotated_cell_data
        else{}// *** NOT IMPLEMENTED YET ***
        identify_cell_types(singularity_key_getter.out.singularity_key_got, unannotated_cells, params.cell_analysis_metadata)
    }    
    if(!params.skip_cell_clustering){
        if(!params.skip_cell_type_identification)
            annotated_cells = identify_cell_types.out.annotated_cell_data
        else{}// *** NOT IMPLEMENTED YET ***
        clustering_metadata = channel.fromPath(params.cell_analysis_metadata)
            .splitCsv(header:true)
            .map{row -> tuple(row.cell_type, row.clustering_markers, row.clustering_resolutions)}
            .filter{!it.contains("NA")}
        cluster_cells(singularity_key_getter.out.singularity_key_got, annotated_cells, clustering_metadata)
    }    
    if(!params.skip_visualization){
        categories = channel.fromPath(params.sample_metadata_file)
            .splitCsv(header:false)
            .first()
            .map{row -> row.drop(2).join(",")}
        if(!params.skip_area){
            visualize_areas(singularity_key_getter.out.singularity_key_got, categories, measure_areas.out.area_measurements, params.sample_metadata_file)
        }
        if(!params.skip_cell_type_identification){
            if(!params.skip_segmentation){
                cell_masks = segment_cells.out.cell_mask_tiffs.collect()
            }
            else{}// *** NOT IMPLEMENTED YET ***
            visualize_cell_types(singularity_key_getter.out.singularity_key_got, categories, identify_cell_types.out.annotated_cell_data,
                params.sample_metadata_file, params.cell_analysis_metadata, cell_masks)
        }
        if(!params.skip_cell_clustering){
            if(!params.skip_cell_type_identification){
                clustered_cell_file = cluster_cells.out.clustered_cell_data
            }
            else{}// *** NOT IMPLEMENTED YET ***
            cluster_visualization_metadata = channel.fromPath(params.cell_analysis_metadata)
                .splitCsv(header:true)
                .map{row -> tuple(row.cell_type, row.clustering_markers, row.clustering_resolutions)}
                .filter{!it.contains("NA")}
            visualize_cell_clusters(singularity_key_getter.out.singularity_key_got, categories, cluster_visualization_metadata,
                clustered_cell_file, params.sample_metadata_file)
        }
    }
}



//raw_tiff_metadata_by_sample.into{raw_tiff_metadata_to_collect; raw_tiff_metadata_to_normalize}



//raw_tiff_metadata_to_normalize
//    .map{file ->
//            def key = file.name.toString().tokenize('-').get(0)
//            return tuple(key, file)
//         }
//    .groupTuple()
//    .set{normalize_tiff_metadata}



/* Collects all the cp3_normalized_tiff_metadata_by_sample metadata files and
    pairs them with the corresponding sample name
*/

//cp3_normalized_tiff_metadata_by_sample
//    .map { file ->
//            def key = file.name.toString().tokenize('-').get(0)
//            return tuple(key, file)
//            }
//    .groupTuple()
//    .set{ cp3_preprocessing_metadata }


//cp3_preprocessing_pipeline_folder = file(params.cp3_preprocessing_cppipe).getParent()
//cp3_preprocessing_pipeline = file(params.cp3_preprocessing_cppipe).getName() 


/* Measures the user specified areas and their ratios:
    - Measures the areas specified in: $params.area_measurements_metadata
    - Copies the results to "$params.output_folder/area_measurements.csv" 
*/

//measurement_metadata = Channel.fromPath(params.area_measurements_metadata)

/* Reads the $params.raw_metadata_file to extract the names of the columns indicating the categories
    to compare (all columns after the 4th)
*/



/* Collects all the cp3_preprocessed_tiff_metadata_by_sample metadata files and
    pairs them with the corresponding sample name
*/

//cp3_preprocessed_tiff_metadata_by_sample
//    .flatten()
//    .map { file ->
//            def key = file.name.toString().tokenize('-').get(0)
//            return tuple(key, file)
//            }
//    .groupTuple()
//    .set{ cp3_segmentation_metadata }


//cp3_segmentation_pipeline_folder = file(params.cp3_segmentation_cppipe).getParent()
//cp3_segmentation_pipeline = file(params.cp3_segmentation_cppipe).getName() 




//annotated_cell_data.into{cell_data_to_plot; cell_data_to_cluster}



/* Reads the clustering metadata file line by line to extract:
   It expects an header line and it extracts the following fields into the sample_metadata channel:
   - cell_type
   - markers: "@" separated list of valid measurement names from $params.output_folder/annotated_cells.csv
   - resolutions: "@" separated list of floats
   If any of these fields is "NA" the line is skipped.
*/

//Channel
//    .fromPath(params.cell_analysis_metadata)
//    .splitCsv(header:true)
//    .map{row -> tuple(row.cell_type, row.clustering_markers, row.clustering_resolutions)}
//    .filter{!it.contains("NA")}
//    .combine(cell_data_to_cluster)
//    .set{clustering_metadata}



//Channel
//    .fromPath(params.cell_analysis_metadata)
//    .splitCsv(header:true)
//    .map{row -> tuple(row.cell_type, row.clustering_markers, row.clustering_resolutions)}
//    .filter{!it.contains("NA")}
//    .combine(clustered_cell_data)
//    .set{cell_visualization_metadata}

