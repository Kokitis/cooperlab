from pathlib import Path
from dataclasses import dataclass
import pendulum
import re
from typing import List, Union


@dataclass
class AssemblyReportMetadata:
	assembly_name: str
	organism_name: str
	infraspecific_name: str
	taxonomic_id: str
	biosample_id: str
	bioproject_id: str
	submitter: str
	date: pendulum.Pendulum
	assembly_type: str
	release_type: str
	assembly_level: str
	genome_representation: str
	wgs_project: str
	assembly_method: str
	expected_final_version: str
	genome_coverage: str
	sequencing_technology: str
	genebank_assembly_extension: str
	refseq_assembly_accession: str
	refseq_and_genebank_asseblies_identical: bool

	@classmethod
	def from_dict(cls, **args):
		# Assembly name:  ASM222304v1
		# Organism name:  Burkholderia sp. HI4860 (b-proteobacteria)
		# Infraspecific name:  strain=HI4860
		# Taxid:          2015361
		# BioSample:      SAMN07257396
		# BioProject:     PRJNA391085
		# Submitter:      University of Michigan
		# Date:           2017-7-20
		# Assembly type:  n/a
		# Release type:   major
		# Assembly level: Contig
		# Genome representation: full
		# WGS project:    NKFO01
		# Assembly method: SPAdes v. 3.7.0
		# Expected final version: yes
		# Genome coverage: 150x
		# Sequencing technology: Illumina HiSeq 4005
		# GenBank assembly accession: GCA_002223045.1
		# RefSeq assembly accession: GCF_002223045.1
		# RefSeq assembly and GenBank assemblies identical: yes
		#
		## Assembly-Units:
		## GenBank Unit Accession	RefSeq Unit Accession	Assembly-Unit name
		## GCA_002223065.1	GCF_002223065.1	Primary Assembly
		#
		# Ordered by chromosome/plasmid; the chromosomes/plasmids are followed by
		# unlocalized scaffolds.
		# Unplaced scaffolds are listed at the end.
		# RefSeq is equal or derived from GenBank object.

		date = args.get('Date')
		if date:
			date = pendulum.parse(date)

		metadata = [
			("assembly_name", args.get('Assembly name')),
			("organism_name", args.get('Organism name')),
			("infraspecific_name", args.get('Infraspecific name')),
			("taxonomic_id", args.get('Taxid')),
			("biosample_id", args.get('BioSample')),
			("bioproject_id", args.get('BioProject')),
			("submitter", args.get('Submitter')),
			("date", date),
			("assembly_type", args.get('Assembly type')),
			("release_type", args.get('Release type')),
			("assembly_level", args.get('Assembly level')),
			("genome_representation", args.get('Genome representation')),
			("wgs_project", args.get('WGS project')),
			("assembly_method", args.get('Assembly method')),
			("expected_final_version", args.get('Expected final version')),
			("genome_coverage", args.get('Genome coverage')),
			("sequencing_technology", args.get('Sequencing technology')),
			("genebank_assembly_extension", args.get('GenBank assembly accession')),
			("refseq_assembly_accession", args.get('RefSeq assembly accession')),
			("refseq_and_genebank_asseblies_identical",
			 args.get('RefSeq assembly and GenBank assemblies identical') == 'yes')
		]
		metadata = [i[1] for i in metadata]

		return cls(*metadata)


class NCBIPackage:
	"""
	Parameters
	----------
	path: Path
		The folder containing the genome assembly and report.
	"""

	def __init__(self, path: Path):
		assert path.is_dir()
		self.genome_id = path.name
		print(repr(self.genome_id))

		print([i.name.replace(self.genome_id, '') for i in path.iterdir()])

	suffixes = [
		'annotation_hashes.txt', 'assembly_status.txt', '_assembly_report.txt', '_assembly_stats.txt',
		'_cds_from_genomic.fna.gz', '_feature_count.txt.gz', '_feature_table.txt.gz', '_genomic.fna', '_genomic.fna.gz',
		'_genomic.gbff.gz', '_genomic.gff.gz', '_protein.faa.gz', '_protein.gpff.gz', '_rna_from_genomic.fna.gz',
		'_translated_cds.faa.gz', '_wgsmaster.gbff.gz', 'md5checksums.txt'
	]

	def unpack(self, key):
		pass


class NCBIReport:
	def __init__(self, path: Path):
		self.sample_id = path.stem.replace('_assembly_report', '')
		with path.open('r') as ncbi_report:
			contents = ncbi_report.read().split('\n')
			self.metadata, self.sequences = self.parse(contents)
	def parse(self, contents: List):
		metadata = [i for i in contents if (i.startswith('#') and len(i) > 2)]

		sequences = [i for i in contents if not i.startswith('#') and len(i) > 2]

		metadata = self.parse_metadata(metadata)

		return metadata, sequences

	def parse_metadata(self, contents):
		regex = re.compile("#\s(.+)[:](.*)")

		matches = [regex.search(i) for i in contents]
		matches = [m.groups() for m in matches if m]
		matches = [(i[0].strip(), i[1].strip()) for i in matches]

		metadata = dict(matches)

		metadata = AssemblyReportMetadata.from_dict(**metadata)
		return metadata


if __name__ == "__main__":
	path = Path.home() / "Documents" / "projects" / "docusate" / "genome_assemblies" / "ncbi-genomes" / "GCF_002223045.1_ASM222304v1" / "GCF_002223045.1_ASM222304v1_assembly_report.txt"
	report = NCBIReport(path)

	print(report.metadata)