import sys
import os
import argparse
import csv
from imctools.io import txtparser
from imctools.io import mcdparser
import numpy

def get_arguments():	
	parser = argparse.ArgumentParser(description = "Converts raw .mcd/.txt files into single tiff images")
	parser.add_argument("sample_name", help = "name of the sample")
	parser.add_argument("roi_name", help = "name of the ROI")
	parser.add_argument("raw_path", help = "path to the raw data .mcd/.txt file")
	parser.add_argument("single_tiff_path", help = "path to a folder for single tiff output")
	parser.add_argument("channel_metadata_file", help = "path to a .csv file with the metadata")
	parser.add_argument("output_metadata_file", help = "path to a .csv file with the metadata")
	args = parser.parse_args()
	return args.sample_name, args.roi_name, args.raw_path, args.single_tiff_path, args.channel_metadata_file, \
		args.output_metadata_file

def parse_channel_metadata(file_name):
	with open(file_name) as metadata_file:
		reader = csv.DictReader(metadata_file)
		if "channel_metal" not in reader.fieldnames or "channel_label" not in reader.fieldnames:
			print(f'"channel_metal" or channel_label missing in the channel metadata header:')
			print("Channel metadata header: " + " ".join(reader.fieldnames))
			sys.exit(1)
		channel_metadata = {row["channel_metal"] : row["channel_label"] for row in reader} 
	return channel_metadata

def get_txt_aquisition(file_name):
	txt = txtparser.TxtParser(file_name)
	return txt.get_imc_acquisition()

def get_mcd_aquisition(file_name, roi_name):
	mcd = mcdparser.McdParser(file_name)
	roi_ids = {mcd.get_acquisition_description(id) : id for id in mcd.acquisition_ids}
	if roi_name not in roi_ids.keys():
		print(f'ROI not present: {roi_name}')
		print("Available ROIs: " + "\n".join(roi_ids.keys()))
		sys.exit(1)
	return mcd.get_imc_acquisition(roi_ids[roi_name])

def output_tiffs(acquisition, output_path, sample_name, roi_name, channel_metadata):
	single_metadata = []
	for metal, label in channel_metadata.items():
		tiff_file_name = os.path.abspath(os.path.join(output_path, sample_name + "-" + label + "-raw.tiff"))
		tiff_image_writer = acquisition.get_image_writer(tiff_file_name, metals = [metal], mass = None)
		tiff_image_writer.save_image(mode = "imagej", compression = 0, dtype = numpy.int16().dtype, bigtiff = False)
		single_metadata.append({"sample_name" : sample_name, "roi_name" : roi_name, "metal" : metal, "label" : label,
			"raw_tiff_file_name" : tiff_file_name})
	return single_metadata	
					
if __name__ == "__main__":
	sample_name, roi_name, raw_path, single_tiff_path, channel_metadata_file, output_metadata_file = get_arguments()
	if not os.path.isfile(raw_path):
		print(f'The raw data file path specified does not exist: {raw_path}')
		sys.exit(1)
	if not os.path.isdir(single_tiff_path):
		print(f'The single tiff output path specified does not exist: {single_tiff_path}')
		sys.exit(1)
	if not os.path.isfile(channel_metadata_file):
		print(f'The channel metadata file specified does not exist: {channel_metadata_file}')
		sys.exit(1)
	channel_metadata = parse_channel_metadata(channel_metadata_file)
	extension = raw_path.split(sep = ".")[-1]
	if extension == "txt":
		acquisition = get_txt_aquisition(raw_path)
	elif extension == "mcd":
		acquisition = get_mcd_aquisition(raw_path, roi_name)
	else:
		print(f'The raw data file extension is not .mcd or .txt: {extension}')
		sys.exit(1)
	single_tiff_metadata = output_tiffs(acquisition, single_tiff_path, sample_name, roi_name, channel_metadata)
	with open(output_metadata_file, mode = 'w') as output_metadata:
		metadata_writer = csv.DictWriter(output_metadata, fieldnames = single_tiff_metadata[0].keys())
		metadata_writer.writeheader()
		metadata_writer.writerows(single_tiff_metadata)
