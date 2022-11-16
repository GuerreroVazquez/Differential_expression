import csv
import re
from itertools import combinations

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from adjusting_genes import *
import ensmbl_finder


def get_data(all_experiments_de='../Results/DSeq/All_experiments_dE.csv'):
    """
    This function takes the csv in the shape
    EX1, Ex2, Ex3
    gen11, gen24, gen36
    gen12, gen25, gen37
    gen13, , gen38

    where [gen11, gen12 gen13] are genes differentiated expressed from EX1
    (i generate this file manually by copy and pasting the dif expressed genes in Calc)
    :param all_experiments_de:
    :return:
    """
    df = pd.read_csv(all_experiments_de)
    experiment_names = list(df.head())
    return df, experiment_names


def save_gene_intersections(all_experiments_de='../Results/DSeq/All_experiments_dE.csv',
                            file_name='../Results/DSeq/intersection_union.csv'):
    """
    This takes ALL the combinations of Experiments in all_experiments_de and calculates
    the union and intersections and saves them on a csv called file_name y regresa la
    cardinalidad de los conjuntos
    :return:
    """
    df, experiment_names = get_data(all_experiments_de)
    u = {}
    cardinalities = {}
    for j in range(2, len(experiment_names) + 1):
        for i in combinations(experiment_names, j):
            a = set(df[i[0]].dropna())
            inter = a
            union = a
            key_name = i[0]
            for k in range(1, j):
                b = set(df[i[k]].dropna())
                inter = inter & b
                union = inter.union(b)
                key_name += "*" + i[k]

            u[key_name.replace("*", "n")] = list(inter)
            u[key_name.replace("*", "u")] = list(union)
            cardinalities[key_name.replace("*", "n")] = len(inter)
            cardinalities[key_name.replace("*", "u")] = len(union)
            if len(inter) == 1:
                print(inter)
    # analysis = pd.DataFrame.from_dict(u)
    analysis = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in u.items()]))
    analysis.to_csv(file_name, index=False)

    return cardinalities


def print_venn(all_genes, names, output='../Results/DSeq/venn_diagram.png', ignore_version=True):
    """
    Makes a venn diagram of THREE sets
    :param all_genes: list of lists of str of all gene names
    :param names: list (str) the names of the experiments
    :param output:
    :return:
    """
    import re
    sets = []
    for set_name in names:
        gene_list = all_genes[set_name]
        set_values = set(gene_list.dropna())
        if ignore_version:
            set_values = {re.sub('\..*$', '', x) for x in set_values}
        sets.append(set_values)
    venn3(sets, names)
    plt.savefig(output)


def transforme_2_genes(data, names, ignore_version=True):
    """
    If it is not empty, removes the version if true,
    it will go to ensbml and find the gene name (display name)
    :param data: pandas dataframe with the experiments and the transcripts
    :param names: list of names of the experiments
    :return: dataframe  with the experiments and the genes
    """
    new_df = data
    for name in names[1:]:
        transcripts = data[name]
        gene_names = []
        for gene in transcripts:
            if ignore_version:
                try:
                    gene = re.sub('\..*$', '', gene)
                except:
                    gene_names.append("NF")
                    continue
            gene_names.append(ensmbl_finder.get_gene_name_from_ensmbl(gene))
        new_df[name] = gene_names
        new_df.to_csv(f"{name}_gene_names.csv")
    return new_df


def convert_2_entrez(genes_data_frame, names, ignore_version=True):
    """
    Gets the entrez number for TRANSCRIPTS
    :param genes_data_frame:
    :param names:
    :return:
    """
    new_df = genes_data_frame
    for name in names[1:]:
        transcripts = genes_data_frame[name]
        gene_names_1 = []
        for gene in transcripts:
            if ignore_version:
                try:
                    gene_names_1.append(re.sub('\..*$', '', gene))
                except:
                    continue
        gene_name = get_gene_name_from_entrez(gene_names_1)
        new_df[name] = gene_name
        new_df.to_csv(f"{name}_gene_entrez_names.csv")
    return new_df


def remove_version(genes_data_frame, experiment_names, diff_express=True,
                   output="gene_ensmbl_no_version_names"):
    """
    Creates a new file with the same genes but without the version
    :param genes_data_frame: The dataframe with the gene names
    :param experiment_names:
    :param diff_express: bool If the genes are diff express or not
    :param output: The prefix of the name the output file will have

    :return:
    """
    new_df = genes_data_frame
    for name in experiment_names:
        transcripts = genes_data_frame[name]
        gene_names_1 = []
        for gene in transcripts:
            try:
                gene_names_1.append(re.sub('\..*$', '', gene))
            except:
                gene_names_1.append("")
                continue
        new_df[name] = gene_names_1
    if diff_express:
        new_df.to_csv(f"{output}NDE.csv")
    else:
        new_df.to_csv(f"{output}DE.csv")
    return new_df


if __name__ == '__main__':
    data, names = get_data(all_experiments_de='../Results/DSeq/All_experiments_NodE.csv')
    # transforme_2_genes(data, names)
    remove_version(data, names, diff_express=False)
    data, names = get_data(all_experiments_de='../Results/DSeq/All_experiments_dE.csv')
    # transforme_2_genes(data, names)
    remove_version(data, names, diff_express=True)

if __name__ == '__main__':
    print(save_gene_intersections(all_experiments_de='../Results/DSeq/All_dE_gene_names.csv',
                                  file_name='../Results/DSeq/intersection_union_gene_name.csv'))
    data, names = get_data(all_experiments_de='../Results/DSeq/All_dE_gene_names.csv')
    print_venn(data, names, output="test.png")
    pass
