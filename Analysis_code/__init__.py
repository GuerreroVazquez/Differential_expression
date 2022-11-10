import csv
from itertools import combinations

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3
from adjusting_genes import *

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
    cardinalities={}
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
            if len(inter)==1:
                print(inter)
    # analysis = pd.DataFrame.from_dict(u)
    analysis = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in u.items()]))
    analysis.to_csv(file_name, index=False)

    return cardinalities


def print_venn(all_genes, names, output='../Results/DSeq/venn_diagram.png'):
    """
    Makes a venn diagram of THREE sets
    :param all_genes:
    :param names:
    :param output:
    :return:
    """

    sets = []
    for set_name in names:
        gene_list = all_genes[set_name]
        sets.append(set(gene_list.dropna()))
    venn3(sets, names)
    plt.savefig(output)

def transforme_2_genes(data, names):
    """

    :param data: pandas dataframe with the experiments and the transcripts
    :param names: list of names of the experiments
    :return: dataframe  with the experiments and the genes
    """
    new_df = data
    for name in names:
        transcripts = data[name]
        gene_names=[]
        for i in range(0, len(transcripts), 113):
            gene_names.append( retrieve_annotation(list(transcripts[i:i+113]))  )
        new_df[name]= gene_names
    return new_df

if __name__ == '__main__':
    data, names = get_data()
    transforme_2_genes(data, names)
if __name__ == '__main__':
    # print(save_gene_intersections())
    data, names = get_data()
    print_venn(data, names)
    pass
