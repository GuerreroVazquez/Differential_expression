import csv
import os
import pandas as pd


def read_sample(sample_name, experiment, colums=['target_id', 'est_counts']):
    file_name = f"/home/karen/Documents/phd/DifferentialExpression/{experiment}/{sample_name}/{sample_name}/abundance.tsv"
    sample = pd.read_csv(file_name, sep='\t')
    filtered_sample = sample[colums].set_index(colums[0]).rename(columns={colums[1]: sample_name})

    return filtered_sample


def write_dictionary_to_csv(file_name="x.csv", dic=None):
    if dic is None:
        dic = {}
    with open(file_name, mode='w') as outfile:
        writer = csv.writer(outfile)
        mydict = {rows[0]: rows[1] for rows in dic}


def merge_files(sample_names=[], experiment="GSE152558", colums=['target_id', 'est_counts']):
    main_dic = pd.DataFrame()
    for sample_name in sample_names:
        sample = read_sample(sample_name, experiment=experiment)
        result = pd.concat([main_dic, sample], axis=1)
        main_dic = result
    return main_dic.rename(columns={colums[0]: 'ensgene'})


experiment = "GSE164471"
sample_names = ["SRR12021926", "SRR12021927", "SRR12021928", "SRR12021929", "SRR12021930"]
sample_names = ["SRR13388732", "SRR13388737", "SRR13388742", "SRR13388747", "SRR13388752", "SRR13388733", "SRR13388738",
                "SRR13388743", "SRR13388748", "SRR13388753", "SRR13388734", "SRR13388739", "SRR13388744", "SRR13388749",
                "SRR13388754", "SRR13388735", "SRR13388740", "SRR13388745", "SRR13388750", "SRR13388755", "SRR13388736",
                "SRR13388741", "SRR13388746", "SRR13388751", "SRR13388756"]
x = read_sample(sample_name=sample_names[0], experiment=experiment)
print(x[:3])

final = merge_files(sample_names=sample_names, experiment=experiment)
final.to_csv(f'/home/karen/Documents/phd/DifferentialExpression/{experiment}/All_samples.csv')
print(final[:5])
