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

include {cluster_cells} from "$script_folder/processes.nf"
include {collect_clustering_data} from "$script_folder/processes.nf"
include {cell_cluster_visualization} from "$script_folder/processes.nf"

workflow singularity_key_getter {
    get_singularity_key()
    emit:
        singularity_key_got = get_singularity_key.out.singularity_key_got 
}

/* Reads the raw metadata file line by line to extract the sample metadata for the raw IMC acquisition files.
   It expects an header line and it extracts the following fields into the sample_metadata channel:
   - sample_name
   - roi_name
   - raw_path -> Converted to a file type 
*/


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

workflow {
    sample_names = channel.fromPath(params.sample_metadata_file).splitCsv(header: true).map{row -> row.sample_name}
    singularity_key_getter()
    if(params.raw_metadata_file && !params.skip_conversion)
        convert_raw_data(singularity_key_getter.out.singularity_key_got) 
    if((params.converted_metadata_file || !params.skip_conversion) && !params.skip_normalization)
        normalization_metadata = (params.skip_conversion) ? channel.fromPath(params.converted_metadata_file) : convert_raw_data.out.converted_tiff_metadata
        normalize_images(singularity_key_getter.out.singularity_key_got, sample_names, normalization_metadata) 
    if(!params.skip_preprocessing)
        if(!params.skip_normalization)
            preprocessing_metadata = normalize_images.out.cp3_normalized_tiff_metadata_by_sample
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

//Channel
//    .fromPath(params.raw_metadata_file)
//    .splitCsv(header:false)
//    .first()
//    .map{row -> row.drop(4).join(",")}
//    .into{categories_area; categories_type; categories_cluster}


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

