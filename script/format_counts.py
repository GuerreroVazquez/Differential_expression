import csv
import os
import pandas as pd


def read_sample(sample_name, colums=['target_id', 'est_counts']):
    file_name = f"/home/karen/Documents/phd/DifferentialExpression/GSE152558/{sample_name}/abundance.tsv"
    sample = pd.read_csv(file_name, sep='\t')
    filtered_sample = sample[colums].set_index(colums[0]).rename(columns={colums[1]: sample_name})

    return filtered_sample


def write_dictionary_to_csv(file_name="x.csv", dic=None):
    if dic is None:
        dic = {}
    with open(file_name, mode='w') as outfile:
        writer = csv.writer(outfile)
        mydict = {rows[0]: rows[1] for rows in dic}


def merge_files(sample_names=[], colums=['target_id', 'est_counts']):
    main_dic = pd.DataFrame()
    for sample_name in sample_names:
        sample = read_sample(sample_name)
        result = pd.concat([main_dic, sample], axis=1)
        main_dic = result
    return main_dic.rename(columns={colums[0]: 'ensgene'})


sample_names = ["SRR12021926", "SRR12021927", "SRR12021928", "SRR12021929", "SRR12021930"]
x = read_sample(sample_name=sample_names[0])
print(x[:3])

final = merge_files(sample_names=sample_names)
final.to_csv('/home/karen/Documents/phd/DifferentialExpression/GSE152558/All_samples.csv')
print(final[:5])
