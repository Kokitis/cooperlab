from pathlib import Path
from typing import List, Dict, Union, Iterable, Tuple
import re
from dataclasses import dataclass
import yaml
import json
from Bio import SeqIO

Row = List[str]

regex = "([\w]+)[=]([^\t]+)"
regex = re.compile(regex)

TYPE_SPECIFIC_FIELDS = {
	'SNP':  ('seq_id', 'position', 'new_seq'),
	'SUB':  ('seq_id', 'position', 'size', 'new_seq'),
	'DEL':  ('seq_id', 'position', 'size'),
	'INS':  ('seq_id', 'position', 'new_seq'),
	'MOB':  ('seq_id', 'position', 'repeat_name', 'strand', 'duplication_size'),
	'AMP':  ('seq_id', 'position', 'size', 'new_copy_number'),
	'CON':  ('seq_id', 'position', 'size', 'region'),
	'INV':  ('seq_id', 'position', 'size'),
	'RA':   ('seq_id', 'position', 'insert_position', 'ref_base', 'new_base'),
	'MC':   ('seq_id', 'start', 'end', 'start_range', 'end_range'),
	'JC':   ('side_1_seq_id',
			 'side_1_position',
			 'side_1_strand',
			 'side_2_seq_id',
			 'side_2_position',
			 'side_2_strand',
			 'overlap'),
	'UN':   ('seq_id', 'start', 'end'),
	'CURA': ('expert',),
	'FPOS': ('expert',),
	'PHYL': ('gd',),
	'TSEQ': ('seq_id', 'primer1_start', 'primer1_end', 'primer2_start', 'primer2_end'),
	'PFLP': ('seq_id', 'primer1_start', 'primer1_end', 'primer2_start', 'primer2_end'),
	'RFLP': ('seq_id', 'primer1_start', 'primer1_end', 'primer2_start', 'primer2_end', 'enzyme'),
	'PFGE': ('seq_id', 'restriction_enzyme'),
	'NOTE': ('note',)
}


@dataclass
class Metadata:
	pass


@dataclass
class VcfRecord:
	CHROM: str
	POS: int
	ID: str
	REF: str
	ALT: str
	QUAL: str
	FILTER: str
	INFO: List[str]

	def __str__(self):
		vinfo = ';'.join(self.INFO) if self.INFO else '.'
		string = "\t".join(map(str, [self.CHROM, self.POS, self.ID, self.REF, self.ALT, self.QUAL, self.FILTER, vinfo]))
		return string


class Mutation:
	def __init__(self, *fields, **named_fields):
		mtype, evidence_id, parent_id, *other_fields = fields
		self.type = mtype[1]
		self.id = evidence_id[1]
		self.parent_id = parent_id[1].split(',')
		self.fields = other_fields
		self.fields_dict = dict(self.fields)
		self.named_fields = named_fields

		self.position = self.fields_dict.get('position')
		if self.position: self.position = int(self.position)
		self.size = self.fields_dict.get('size')
		if self.size: self.size = int(self.size)

		self.seq_id = self.fields_dict.get('seq_id', '')
		self.new_seq = self.fields_dict.get('new_seq', '')

	def __str__(self):
		_t = sorted("=".join(i) for i in self.named_fields.items())
		_f = ["=".join(i) for i in self.fields]
		_pi = ",".join(self.parent_id)
		string = "Mutation('{}', '{}', '{}', {}, {})".format(self.type, self.id, _pi, _f, _t)
		return string

	def to_dict(self):
		data = dict(self.fields)
		data['type'] = self.type
		data['id'] = self.id
		data['parent_id'] = self.parent_id
		data = {**data, **self.named_fields}
		return data

	def get(self, item):
		return self.fields_dict.get(item, self.named_fields.get(item))

	def to_vcf(self, reference: str):
		if self.type == 'INS':
			reference_sequence = get_position(reference, self.position)
			alternate_sequence = reference_sequence + self.new_seq
		elif self.type == 'SNP':
			reference_sequence = get_position(reference, self.position)
			alternate_sequence = self.new_seq
		elif self.type == 'SUB':
			reference_sequence = get_position(reference, (self.position, self.position + self.size))
			alternate_sequence = self.new_seq
		elif self.type == 'DEL':
			reference_sequence = get_position(reference, (self.position, self.position + self.size))
			alternate_sequence = '.'
		elif self.type == 'INV':
			reference_sequence = get_position(reference, (self.position, self.position + self.size))
			alternate_sequence = reference_sequence[::-1]
		else:
			message = "'{}' Invalid mutation type!".format(self.type)
			raise ValueError(message)

		vcf_pass = "PASS"
		vcf_info = [
			('TP', self.type),
			('P', len(self.parent_id)),
			('GP', self.get('gene_position')),
			# ('GN', self.get('gene_name')),
			# ('LT', self.get('locus_tag')),
			('CAT', self.get('mutation_category'))
		]
		vcf_info = sorted(vcf_info)
		vcf_info = ['{}={}'.format(i, j).replace(' ', '') for i, j in vcf_info]

		record = VcfRecord(self.seq_id, self.position, '.', reference_sequence, alternate_sequence, vcf_pass, '.',
						   vcf_info)

		return record


def get_position(genome: str, position: Union[int, Iterable[int]]) -> str:
	if isinstance(position, str): position = int(position)
	if isinstance(position, int):

		result = genome[position - 1]
	elif isinstance(position, (list, tuple)) and len(position) == 2:
		position = slice(position[0] - 1, position[1] - 1)
		result = genome[position].seq
	else:
		message = "Invalid position: '{}'".format(position)
		raise ValueError(message)
	return result


class GenomeDiff:
	""" http://barricklab.org/twiki/pub/Lab/ToolsBacterialGenomeResequencing/documentation/gd_format.html
		Parameters
		----------
		path: Path
			Path to the .gd file
	"""

	def __init__(self, path: Path):
		with path.open('r', encoding='utf-8') as gd_file:
			contents = gd_file.read().split('\n')
			contents = [r.split('\t') for r in contents]

		self.mutations, self.evidence = self.parse(contents)

		self.mutation_map = {(key.seq_id, key.position):index for index, key in enumerate(self.mutations)}
		self.evidence_map = {key.id: index for index, key in enumerate(self.evidence)}
	def __getitem__(self, item:Tuple[str,int]):
		return self.mutation_map.get(item)
	def get_evidence(self, mutation: Mutation):
		evidence = (self.evidence_map[i] for i in mutation.parent_id)
		evidence = [self.evidence[i] for i in evidence]

		return evidence

	def parse(self, contents: List[Row]):
		evidence = list()
		mutations = list()
		for line in contents:
			if len(line) < 1:
				continue

			row_type = line[0]

			if len(row_type) == 3:  # Is a mutation
				result = self.parse_row(line)
				mutations.append(result)
			elif len(row_type) == 2:
				result = self.parse_row(line)
				evidence.append(result)
			else:
				pass

		return mutations, evidence

	def parse_row(self, row: Row) -> Mutation:
		mutation_type = row[0]
		field_keys = TYPE_SPECIFIC_FIELDS[mutation_type]
		field_keys = ['type', 'id', 'parent_id'] + list(field_keys)
		fields = [(f, r) for f, r in zip(field_keys, row[:len(field_keys) + 3])]

		named_fields = self.get_named_fields(row)
		mutation = Mutation(*fields, **named_fields)
		return mutation

	@staticmethod
	def get_named_fields(row: Row) -> Dict[str, str]:
		named_fields = [regex.search(i) for i in row]

		named_fields = [i for i in named_fields if i]
		named_fields = dict([i.groups() for i in named_fields])
		return named_fields

	def to_yaml(self, path: Path):
		data = dict()
		for mutation in self.mutations:
			evidence = self.get_evidence(mutation)
			data[mutation.id] = {
				'mutation': mutation.to_dict(),
				'evidence': [i.to_dict() for i in evidence]
			}
		if path.suffix == '.yaml':
			with path.open('w') as yaml_file:
				yaml.dump(data, yaml_file, default_flow_style = False)
		else:
			with path.open('w') as file1:
				file1.write(json.dumps(data, sort_keys = True, indent = 4))

	def to_vcf(self, reference: Path, path: Path = None):
		lines = list()
		reference = SeqIO.parse(reference, "fasta")
		reference = SeqIO.to_dict(reference)
		for mutation in self.mutations:
			ref_seq = reference[mutation.seq_id]
			line = mutation.to_vcf(ref_seq.seq)
			lines.append(str(line))

		if path:
			with path.open('w') as vcf_file:
				vcf_file.writelines(lines)
		print("\n".join(lines))
		return lines


if __name__ == "__main__":
	path = Path.home() / "Documents" / "output.gd"
	path = Path.home() / "Documents" / "projects" / "Moreira-POR" / "Breseq Output" / "P148-1" / "output" / "evidence" / "annotated.gd"
	gd = GenomeDiff(path)
