import csv
from itertools import combinations

import pandas as pd




def get_data(all_experiments_de='../Results/DSeq/All_experiments_dE.csv'):
    """
    This function takes the csv in the shape
    EX1, Ex2, Ex3
    gen11, gen24, gen36
    gen12, gen25, gen37
    gen13, NA, gen38

    where [gen11, gen12 gen13] are genes differentiated expressed from EX1
    (i generate this file manually by copy and pasting the dif expressed genes in Calc)
    :param all_experiments_de:
    :return:
    """
    df = pd.read_csv(all_experiments_de)
    experiment_names = list(df.head())
    return df, experiment_names


def save_gene_intersections(matplotlib_venn=None):
    df, experiment_names = get_data()
    u = {}
    for j in range(2, len(experiment_names)+1):
        for i in combinations(experiment_names, j):
            a = set(df[i[0]])
            inter = a
            union = a
            key_name = i[0]
            for k in range(1, j):
                b = set(df[i[k]])
                inter = inter & b
                union = inter.union(b)
                key_name += "*" + i[k]
            u[key_name.replace("*", "n")] = list(inter)
            u[key_name.replace("*", "u")] = list(union)

    # analysis = pd.DataFrame.from_dict(u)
    analysis = pd.DataFrame(dict([(k, pd.Series(v)) for k, v in u.items()]))
    file_name = '../Results/DSeq/intersection_union.csv'
    analysis.to_csv(file_name, index=False)

    return set(df["GSE152558"]), experiment_names
    list(set(df[0]) & set(df[1]))


def print_venn(all_genes, names, output='../Results/DSeq/venn_diagram.png'):
    import matplotlib.pyplot as plt
    from matplotlib_venn import venn3
    sets = []
    for set_name in names:
        gene_list = all_genes[set_name]
        sets.append(set(gene_list))
    plot = venn3(sets, names)
    plot.savefig(output)


if __name__ == '__main__':
    save_gene_intersections()
    # data, names = get_data()
    # print_venn(data, names)
    pass
