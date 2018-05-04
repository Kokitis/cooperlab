from pathlib import Path
from Bio import SeqIO
import pandas
try:
	from .genome_diff_parser import GenomeDiff
except:
	from genome_diff_parser import GenomeDiff
class Isolate:
	def __init__(self, path:Path):
		assert path.is_dir()
		self.path = path
		self.sample_id = path.stem
		self.output_folder = path / "sample_output"
		if not self.output_folder: self.output_folder.mkdir()
		self.output_vcf = path / "data" / "output.vcf"
		self.output_gd_basic = path / "output" / "output.gd"
		self.output_gd_evidence = path / "output" / "evidence" / "evidence.gd"
		self.output_gd_annotated = path / "output" / "evidence" / "annotated.gd"
		self.index = path / "output" / "index.html"
		self.reference = path / "data" / "reference.fasta"

		self.output_gd_annotated = GenomeDiff(self.output_gd_annotated)
		self.output_gd_evidence = GenomeDiff(self.output_gd_evidence)
		self.output_gd_basic = GenomeDiff(self.output_gd_basic)

		#self.generate_output_table()
		reference = SeqIO.to_dict(SeqIO.parse(self.reference, "fasta"))
		#self.output_gd_annotated.to_vcf(self.reference)


	def combine_output_files(self):
		""" Combines all relevant outputfiles into a single table."""
		pass

	def generate_output_table(self, path:Path=None)->Path:
		output_table = list()
		record_dict = SeqIO.to_dict(SeqIO.parse(self.reference, "fasta"))
		for index, mutation in enumerate(self.output_gd_annotated.mutations):
			seq_id = mutation.get('seq_id')
			ref_seq = record_dict[seq_id]
			vcf_record = mutation.to_vcf(ref_seq)

			row = {
				'sample': self.sample_id,
				#'gene': mutation.get('gene_name'),
				'ref': vcf_record.REF,
				'alt': vcf_record.ALT,
				'position': mutation.position,
				'sequenceId': mutation.seq_id,
				'mutationType': mutation.type
			}
			output_table.append(row)

		output_file = self.output_folder / "annotated_table.tsv"
		df = pandas.DataFrame(output_table)
		df.to_csv(str(output_file), sep = '\t', index = False)
		return output_file

if __name__ == "__main__":
	path = Path.home() / "Documents" / "projects" / "Moreira-POR" / "Breseq Output" / "P148-1"

	isolate = Isolate(path)
	isolate_output = path / "isolate_output"
	if not isolate_output.exists():
		isolate_output.mkdir()

	isolate.generate_output_table(isolate_output / "output_table.tsv")