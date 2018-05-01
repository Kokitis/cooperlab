from pathlib import Path
import pandas
import shutil
import argparse

LOGFILE_PATH = Path.home() / "fastq_concat_log_file.log"

parser = argparse.ArgumentParser(
	description = "Combines fastq files based on the SampleSheet.csv file used during a sequencing run.")

parser.add_argument(
	'-i', '--sheet',
	action = "store",
	help = "The filename of the sample sheet file used during the sequencing run, or a folder containing sample configuration folders with sample sheets",
	dest = "sheet"
)

parser.add_argument(
	'-o', '--output',
	action = "store",
	help = "Name of the output folder. The concatenated fastq files will be saved to this folder.",
	dest = 'output'
)

args = parser.parse_args()

def print_log(string):
	path = LOGFILE_PATH
	with path.open('a') as log_file:
		log_file.write(string + "\n")
	print(string)

def search_for_sample_sheets(folder:Path):
	pattern = "*/SampleSheet.csv"
	return list(folder.glob(pattern))

def groupby(iterable, by):
	groups = dict()
	for i in iterable:
		key = by(i)
		if key not in groups:
			groups[key] = [i]
		else:
			groups[key].append(i)
	return groups


def concatenate_files(output_file, files):
	with output_file.open('wb') as output:
		for f in files:
			with f.open('rb') as fd:
				shutil.copyfileobj(fd, output, 1024 * 1024 * 10)


def combine_files(sample_sheet_path: Path, output_folder: Path):
	"""

	Parameters
	----------
	path: Path
		path to a sample sheet.

	Returns
	-------

	"""

	dmux_folder = Path("/home/dmux")
	sample_sheet = pandas.read_csv(str(sample_sheet_path), skiprows = 9)

	for index, row in sample_sheet.iterrows():
		# Sample_ID	Sample_Name	Species	Project	NucleicAcid	Sample_Well	I7_Index_ID	index	I5_Index_ID	index2
		sample_name = row['Sample_Name']
		candidates = dmux_folder.glob("*/*/*/{}_*.fastq.gz".format(sample_name))
		candidates = sorted(candidates)
		groups = groupby(candidates, lambda s: s.stem.split('_')[-2])
		print(sample_sheet_path.parent.stem, sample_name)

		for key, values in sorted(groups.items()):
			output_filename = output_folder / sample_sheet_path.parent.stem/ "{}_{}.fastq.gz".format(sample_name, key)
			if not output_filename.parent.exists():
				output_filename.parent.mkdir()

			values = sorted(values)
			print("\tCombining '{}-{}' ".format(sample_name, key))
			for v in values:
				print("\t\t", v)
			print("\tOutput File: ", output_filename)

			concatenate_files(output_filename, values)


if __name__ == "__main__":
	data_folder = Path("/home/data")
	dmux_folder = Path("/home/dmux")
	if args.output:
		output_folder = Path(args.output)
	else:
		output_folder = Path.home() / "concatenated_fastq_files"
	LOGFILE_PATH = output_folder / "fastq_concat_log_file.txt"
	if args.sheet:
		sample_sheet = Path(args.sheet)
	else:
		sample_sheet = data_folder / "180416_NB501145_0100_AHJNJ5BGX5" / "SampleSheet.csv"

	if not output_folder.exists():
		output_folder.mkdir()
	if sample_sheet.is_dir():
		sample_sheets = search_for_sample_sheets(sample_sheet)
	else:
		sample_sheets = [sample_sheet]
	for sample_sheet in sample_sheets:
		print("Using {}...".format(sample_sheet))
		if sample_sheet.parent.stem.split('_')[0] not in ['180416', '180423']: continue
		combine_files(sample_sheet, output_folder)
