import csv
import re
from itertools import combinations

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from Analysis_code.adjusting_genes import *
import Analysis_code.ensmbl_finder
import seaborn as sns

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


def print_venn(all_genes, names, output='../Results/DSeq/venn_diagram.png',
               ignore_version=True, alias=None):
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
            set_values = set(filter(None, set_values))
        sets.append(set_values)
    if len(names) == 3:
        if not alias:
            alias = names
        venn3(sets, alias)
        plt.savefig(output)
    elif len(names) == 2:
        import venn
        if not alias:
            alias = names
        labels = venn.get_labels(sets)
        fig, ax = venn.venn2(labels, names=alias)
        fig.savefig(output)
    elif len(names) == 5:
        import venn
        if not alias:
            alias = venn.get_labels(sets)
        labels = venn.get_labels(sets)
        fig, ax = venn.venn5(labels, names=alias)
        fig.savefig(output)


def transforme_2_genes(data, names, ignore_version=True):
    """
    If it is not empty, removes the version if true,
    it will go to ensbml and find the gene name (display name)
    :param data: pandas dataframe with the experiments and the transcripts
    :param names: list of names of the experiments
    :return: dataframe  with the experiments and the genes
    """
    new_df = data
    for name in names:
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
        gene_names_1 = remove_version_from_list(transcripts)
        new_df[name] = gene_names_1
    if output:
        if diff_express:
            new_df.to_csv(f"{output}NDE.csv")
        else:
            new_df.to_csv(f"{output}DE.csv")
    return new_df


def remove_version_from_list(gene_names):
    """
    From a list of gene names, return a list with the same names but without version.
    :param gene_names:
    :return:
    """
    gene_names_1 = []
    for gene in gene_names:
        try:
            gene_names_1.append(re.sub('\..*$', '', gene))
        except:
            gene_names_1.append("")
            continue
    return gene_names_1


def get_genes_from_csv(read_count_file_name="../Results/GSE164471def_expressed.csv"):
    df, _ = get_data(all_experiments_de=read_count_file_name)
    genes = list(df["Gen"])
    genes = remove_version_from_list(genes)
    filter_read_counts(metadata_path="../Results/culo_culo.csv",
                       filered_list=genes, filter_by="Gen", output="../Results/GSE164471_DE_abundances.csv")

    pass


def filter_read_counts(metadata_path, filered_list, filter_by="Gen", output=None):
    """
    Filter from the csv in metadata_path only those that in "filter_by" have the values in
    "filtered_list" and save it in output
    :param metadata_path:
    :param filered_list:
    :param filter_by:
    :param output:
    :return:
    """
    if output is None:
        output = f"{metadata_path[:-4]}_filtered.csv"
    df_metadata = pd.read_csv(metadata_path)
    df_metadata = df_metadata[df_metadata[filter_by].isin(filered_list)]
    df_metadata.to_csv(output, index=False)


def test_filter_culo_culo():
    get_genes_from_csv("../Results/GSE164471def_expressed.csv")


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


def test_24_nov_22():
    """
    I have my genes, 17 for (1), 82 for (2) and 13 for (3).
    So let's see if they have something in common
    :return:
    """
    all_de_data = '../Results/def_expressed.csv'
    data, names = get_data(all_experiments_de=all_de_data)
    data = remove_version(genes_data_frame=data, experiment_names=names)
    print_venn(data, names, output="All_experiments_no_versions.png")
    pass
def test_24_nov_22_venn_3_4_5():
    """
    The venn diagram with all is ugly. And still not a lot of information

    :return:
    """
    all_de_data = '../Results/def_expressed.csv'
    data, names = get_data(all_experiments_de=all_de_data)
    data = remove_version(genes_data_frame=data, experiment_names=names)
    names = names[2:]
    data = data[names]
    labels = ["Young vs Old", "Young vs Middle Age", "Middle Age vs Old"]
    print_venn(data, names, output="Experiments_71_no_versions.png", labels= labels)
    pass
def test_24_nov_22_venn_MAvsOld():
    """
    The venn diagram with all is ugly. And still not a lot of information

    :return:
    """
    all_de_data = '../Results/def_expressed.csv'
    data, names = get_data(all_experiments_de=all_de_data)
    data = remove_version(genes_data_frame=data, experiment_names=names)
    names = [names[0], names[-1]]
    labels = ["GSE152558", 'GSE164471']
    print_venn(data, names, output="OldvsMA_no_versions.png", alias= labels)
    pass

def test_24_nov_get_intersection_n_stuff():
    """
    Here I just get the list of the genes that I go the venn Diagram for
    :return:
    """
    all_de_data = '../Results/def_expressed.csv'
    save_gene_intersections(all_experiments_de=all_de_data,file_name="RNASeq_gene_sets.csv")

def test_24_nov_22_venn_oldvsyoung():
    """
    GSE157585  and GSE164471 have genes in common
    :return:
    """
    all_de_data = '../Results/def_expressed.csv'
    data, names = get_data(all_experiments_de=all_de_data)
    data = remove_version(genes_data_frame=data, experiment_names=names)
    names = [names[1], names[2]]
    labels = ["GSE157585", 'GSE164471']
    print_venn(data, names, output="OldvsYoung_no_versions.png", alias= labels)
    pass

def test_24_nov_22_venn_oldvsyoung():
    """
    GSE157585  and GSE164471 have genes in common
    :return:
    """
    all_de_data = '../Results/def_expressed.csv'
    data, names = get_data(all_experiments_de=all_de_data)
    data = remove_version(genes_data_frame=data, experiment_names=names)
    names = [names[1], names[3]]
    labels = ["GSE157585", 'GSE164471']
    print_venn(data, names, output="MAnOldvsYoung_no_versions.png", alias= labels)
    pass

def test_24_nov_22_venn_MDvsoldvsyoung():
    """
    GSE152558  and GSE164471 also have genes in common, even if
    the first one is MA vs Old and the other one Young vs Old
    GSE157585nGSE164471nGSE164471_MvYoung
    :return:
    """
    all_de_data = '../Results/def_expressed.csv'
    data, names = get_data(all_experiments_de=all_de_data)
    data = remove_version(genes_data_frame=data, experiment_names=names)
    names = [names[0], names[2]]
    labels = ["GSE152558", 'GSE164471']
    print_venn(data, names, output="MAnOldvsYoung_no_versions.png", alias= labels)
    pass



def test_25_nov_22_convert_de_genes_all_2_gene_name():
    """
    What it says in the title
    :return:
    """
    all_de_data = '../Results/def_expressed.csv'
    data, names = get_data(all_experiments_de=all_de_data)

    transforme_2_genes(data=data,names=names)

def create_heatmap(df):

    sns.heatmap(df)


def test_create_heatmap_of_df_genes():
    """
    I have already the genes de for all the experiments, and I have the abundances on
    GSE164471_abundancess.csv
    :return:
    """
    glue1 = sns.load_dataset("glue")

    data, _ = get_data("../Results/GSE164471_DE_abundances.csv")
    glue = glue1.pivot("Model", "Task", "Score")
    data = data.pivot("Gen", "Condition","Abundance")
    #sns.heatmap(glue)
    sns.heatmap(data)


def test_25_nov_fitler():
    """
    Pre meeting: lets filter the txi for those that are de
    :return:
    """
    experiments = ["GSE152558", "GSE157585", "GSE164471", "GSE164471_MAvO", "GSE164471_MAvY"]
    all_de_data = '../Results/def_expressed.csv'
    data, names = get_data(all_experiments_de=all_de_data)
    for experiment in experiments:
        path = f"../Results/{experiment}_txi.csv"
        diff_expressed = data[experiment]
        filter_read_counts(metadata_path=path,filered_list=diff_expressed, filter_by="Gen", output=f"{experiment}_abundances.csv")
