
from pathlib import Path
from typing import List
try:
	from ..misc_parsers import NCBIReport
except:
	from misc_parsers import NCBIReport

def remap_tree_labels(tree_path:Path, reports:List[NCBIReport]):
	reports = [(NCBIReport(i) if not isinstance(i, NCBIReport) else i) for i in reports]
	id_map = list()
	for report in reports:
		name = report.metadata.organism_name
		strain = report.metadata.infraspecific_name

		assembly_name = report.metadata.assembly_name
		accession=report.metadata.refseq_assembly_accession

		value = "{}|{}".format(name, strain).replace('(', '').replace(')', '')
		key = "{}_{}".format(accession, assembly_name)

		id_map.append((key, value))

	with tree_path.open('r') as tree_file:
		contents = tree_file.read()
		for key, value in id_map:
			print("Replacing '{}' with '{}'".format(key, value))
			contents = contents.replace(key, value)
	output_file = tree_path.with_suffix('.labeled.treefile')
	with output_file.open('w') as output:
		output.write(contents)
	print("Wrote content to ", output_file)





import pandas
if __name__ == "__main__":
	pass

