from pyensembl import EnsemblRelease
from pyensembl import Genome
from pyensembl import genome_for_reference_name
data = EnsemblRelease(77)
human = genome_for_reference_name("GRCh38")

gene_name = human.transcript_by_id("ENSG00000235205")
print(gene_name)