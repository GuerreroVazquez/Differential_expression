# This file will do the things that I couldn't do with R becasue R is ugly >:c
import csv
import os
import re

import pandas as pd

master_path = "/home/karen/Documents/phd/DifferentialExpression/"


def metadata_file(experiment):
    """
    This function will give me the path to the metadata
    :param experiment:
    :return:
    """
    count_dir = f"{master_path}{experiment}/metadata.csv"
    return count_dir


def get_existent_samples(experiment):
    files = [f for f in os.listdir(f"{master_path}{experiment}") if re.match(r'SRR', f)]
    return files


def assign_age_group(metadata_path):
    df_metadata = pd.read_csv(metadata_path)
    df_metadata['category'] = 'MiddleAge'
    df_metadata['category'].loc[df_metadata['Age'] < 35] = "Young"
    df_metadata['category'].loc[df_metadata['Age'] > 65] = "Old"
    df_metadata.to_csv(metadata_path, index=False)


experiments = ["GSE152558", "GSE157585", "GSE164471"]


def set_category_all():
    for experiment in experiments:
        metadata_path = metadata_file(experiment)
        assign_age_group(metadata_path)


def filter_metadata(metadata_path, filered_list, output=None):
    if output is None:
        output = f"{metadata_path[:-4]}_filtered.csv"
    df_metadata = pd.read_csv(metadata_path)
    df_metadata = df_metadata[df_metadata['Run'].isin(filered_list)]
    df_metadata.to_csv(output, index=False)


def filter_all(experiments):
    for experiment in experiments:
        metadata_path = metadata_file(experiment)
        filtered_list = get_existent_samples(experiment)
        filter_metadata(metadata_path=metadata_path, filered_list=filtered_list)


def get_filter_by_age(metadata_path, ages=["Young", "Old"]):
    df_metadata = pd.read_csv(metadata_path)
    filtered_df = df_metadata[df_metadata['category'].isin(ages)]
    filtered = list(filtered_df["Run"])
    return filtered


def remove_ma_on_71():
    metadata_path = f"{metadata_file(experiments[2])[:-4]}_filtered.csv"
    only_young_old = get_filter_by_age(metadata_path)
    filter_metadata(metadata_path, only_young_old, metadata_path)


def get_path_list(samples, experiment):
    abundances = []
    for sample in samples:
        abundances.append(f"{master_path}{experiment}/{sample}/{sample}/abundance.tsv")
    return abundances


def save_path_list(path_list, experiment, output=None):
    if output is None:
        output = f"{master_path}{experiment}/abundances.csv"
    n_names = ["{}\n".format(i) for i in path_list]
    with open(output, 'w') as f:
        # using csv.writer method from CSV package
        f.writelines(n_names)


metadata_path = f"{metadata_file(experiments[2])[:-4]}_filtered.csv"
existent = get_existent_samples(experiment=experiments[2])
age_range = get_filter_by_age(metadata_path)
existent_and_age_range = set(existent).intersection(set(age_range))
paths = get_path_list(list(existent_and_age_range), experiment=experiments[2])
save_path_list(path_list=paths, experiment=experiments[2])
