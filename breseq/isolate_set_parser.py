from pathlib import Path
import pandas

try:
	from .isolate_parser import Isolate
except:
	from isolate_parser import Isolate


class IsolateSet:
	def __init__(self, path: Path):
		""" Parses a folder containing a number of breseq output folders."""
		self.output_folder = path / "isolate_set_output"
		if not self.output_folder.exists():
			self.output_folder.mkdir()
		self.samples = list()
		for sample in path.iterdir():
			s = Isolate(sample)
			self.samples.append(s)

	def combineIsolateTables(self):

		tables = list()
		for sample in self.samples:
			df = pandas.read_table(sample.output_folder / "output_table.tsv", sep = "\t")
			tables.append(df)

		df = pandas.concat(tables)
		output_filename = self.output_folder / "isolate_set_combined_table.tsv"
		df.to_csv(output_filename, sep = "\t")
		return output_filename


if __name__ == "__main__":
	path = Path("/home/cld100/projects/moreira_por/output/P148")
	isolate_set = IsolateSet(path)
	isolate_set.combineIsolateTables()
