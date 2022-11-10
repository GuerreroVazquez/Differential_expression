from pyensembl import EnsemblRelease
from pyensembl import Genome
from pyensembl import genome_for_reference_name
# data = EnsemblRelease(77)
# human = genome_for_reference_name("GRCh38")

# gene_name = human.transcript_by_id("ENSG00000235205")
# print(gene_name)

from PyEntrezId import Conversion  # got it from https://github.com/lwgray/pyEntrezId

import sys

from Bio import Entrez

# *Always* tell NCBI who you are
Entrez.email = "your email here"


def retrieve_annotation(id_list):
    """Annotates Entrez Gene IDs using Bio.Entrez, in particular epost (to
    submit the data to NCBI) and esummary to retrieve the information.
    Returns a list of dictionaries with the annotations.

    got it from https://biopython.org/wiki/Annotate_Entrez_Gene_IDs
    """

    request = Entrez.epost("gene", id=",".join(id_list))
    try:
        result = Entrez.read(request)
    except RuntimeError as e:
        # FIXME: How generate NAs instead of causing an error with invalid IDs?
        print("An error occurred while retrieving the annotations.")
        print("The error returned was %s" % e)
        sys.exit(-1)

    webEnv = result["WebEnv"]
    queryKey = result["QueryKey"]
    data = Entrez.esummary(db="gene", webenv=webEnv, query_key=queryKey)
    annotations = Entrez.read(data)

    print("Retrieved %d annotations for %d genes" % (len(annotations), len(id_list)))

    return annotations


def get_entrez_ids_from_ensembl(ensembl_ids=None):
    if ensembl_ids is None:
        ensembl_ids = ['ENST00000407559']
    idc = Conversion('dummyemail@dummybunny.info')
    entrez_ids = []
    for EnsemblId in ensembl_ids:
        entrez_id = idc.convert_ensembl_to_entrez(EnsemblId)
        entrez_ids.append(entrez_id)
    return entrez_ids

def get_gene_name_from_entrez(entrez_ids):
    gene_names=[]
    for i in range(0,len(entrez_ids),100):
        tmp_list = entrez_ids[i:i+100]
        entrez_ids = get_entrez_ids_from_ensembl(tmp_list)
        x = retrieve_annotation(entrez_ids)
        y = x['DocumentSummarySet']
        y = dict(y)['DocumentSummary']
        y = dict(y[0])
        name = y['NomenclatureSymbol']
        gene_names.extend(name)
    return gene_names

#get_gene_name_from_entrez(['ENST00000226004', 'ENST00000407559'])

