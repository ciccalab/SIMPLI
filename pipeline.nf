script_folder = "$baseDir/scripts"
python_environment_folder = "$baseDir/venv"

python_path = "/share/apps/python-3.7.2-shared/bin"
py_library_path = "/share/apps/python-3.7.2-shared/lib"

params.output_folder ="$baseDir/example_output"
params.data_folder = "$baseDir/example_data"
params.channel_metadata = "$params.data_folder/channel_metadata.csv" 
params.raw_metadata_file = "$baseDir/example_data/sample_metadata.csv" 

Channel
    .fromPath(params.raw_metadata_file)
    .splitCsv(header:true)
    .map{row -> tuple(row.sample_name, row.roi_name, file(row.raw_path), row.tiff_path)}
    .set{sample_metadata}

process convert_raw_data_to_single_tiffs {

    publishDir "$tiff_path", mode:'copy', overwrite: true

    input:
        set sample_name, roi_name, raw_path, tiff_path from sample_metadata

    output:
        file "$sample_name*.tiff" into raw_tiff_images
        file "${sample_name}_raw_tiff_metadata.txt" into raw_tiff_metadata_by_sample
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
        . \\
        $params.channel_metadata \\
        ${sample_name}_raw_tiff_metadata.txt > $sample_name/extract_log.txt 2>&1
    """
}

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

