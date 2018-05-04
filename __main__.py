from phylogeny import remap_tree_labels
from pathlib import Path
if __name__ == "__main__":
	projects_folder = Path.home() / "Documents" / "projects" / "docusate"
	treename = Path.home() / "Documents" / "projects" / "docusate" / "iqtree_final_output_quick" / "nasp_full_quickbb.treefile"
	folder = projects_folder / "genome_assemblies" / "outgroups"

	folder2 = projects_folder / "genome_assemblies" / "ncbi-genomes"
	pattern = "*/*_assembly_report.txt"
	reports = list(folder.glob(pattern)) + list(folder2.glob(pattern))


	remap_tree_labels(treename, reports)