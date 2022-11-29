# This file will do the things that I couldn't do with R becasue R is ugly >:c
import csv
import os
import re

import pandas as pd

master_path = "/home/karen/Documents/phd/DifferentialExpression/"

experiments = ["GSE152558", "GSE157585", "GSE164471", "GSE164471_MAvO", "GSE164471_MAvY"]


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


def set_category_all():
    for experiment in experiments:
        metadata_path = metadata_file(experiment)
        assign_age_group(metadata_path)


def filter_metadata(metadata_path, filered_list, output=None):
    """
    Takes the metadata path and a list to be filtered and save it in output
    :param metadata_path:
    :param filered_list:
    :param output:
    :return:
    """
    if output is None:
        output = f"{metadata_path[:-4]}_filtered.csv"
    df_metadata = pd.read_csv(metadata_path)
    df_metadata = df_metadata[df_metadata['Run'].isin(filered_list)]
    df_metadata.to_csv(output, index=False)


def filter_all(experiments):
    """
    This function takes all the experiments and filter only fot the values that exist in the folder.

    :param experiments:
    :return:
    """
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
    """
    This function removes the middle age samples
    :return:
    """
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


def test_24_nov_22_experiment1():
    """
    I just ran this to have the same format as the other experiments,
    I had to filter MiddleAge and Old becasue this one do not have young
    :return:
    """
    metadata_path = f"{metadata_file(experiments[0])[:-4]}_filtered.csv"
    existent = get_existent_samples(experiment=experiments[0])
    age_range = get_filter_by_age(metadata_path, ages=["MiddleAge", "Old"])
    existent_and_age_range = set(existent).intersection(set(age_range))
    paths = get_path_list(list(existent_and_age_range), experiment=experiments[0])
    save_path_list(path_list=paths, experiment=experiments[0])


def test_24_nov_22_experiment3():
    """
    Experimet that ends with 71 is young, old and MA,
    so I did the process for Old and Young, now i do it for MA and old
    :return:
    """
    experiment_number = 3
    metadata_path = metadata_file(experiments[experiment_number])
    filtered_list = get_existent_samples(experiments[experiment_number])
    filter_metadata(metadata_path=metadata_path, filered_list=filtered_list)
    metadata_path = f"{metadata_file(experiments[experiment_number])[:-4]}_filtered.csv"
    existent = get_existent_samples(experiment=experiments[3])
    age_range = get_filter_by_age(metadata_path, ages=["MiddleAge", "Young"])
    existent_and_age_range = set(existent).intersection(set(age_range))
    paths = get_path_list(list(existent_and_age_range), experiment=experiments[experiment_number])
    save_path_list(path_list=paths, experiment=experiments[experiment_number])
    filter_metadata(metadata_path=metadata_path, filered_list=existent_and_age_range, output=metadata_path)


def test_24_nov_22_experiment4():
    """
    Experimet that ends with 71 is young, old and MA,
    so I did the process for Old and Young, now i do it for MA and old
    :return:
    """
    experiment_number = 4
    metadata_path = metadata_file(experiments[experiment_number])
    filtered_list = get_existent_samples(experiments[experiment_number])
    filter_metadata(metadata_path=metadata_path, filered_list=filtered_list)
    metadata_path = f"{metadata_file(experiments[experiment_number])[:-4]}_filtered.csv"
    existent = get_existent_samples(experiment=experiments[3])
    age_range = get_filter_by_age(metadata_path, ages=["MiddleAge", "Old"])
    existent_and_age_range = set(existent).intersection(set(age_range))
    paths = get_path_list(list(existent_and_age_range), experiment=experiments[experiment_number])
    save_path_list(path_list=paths, experiment=experiments[experiment_number])
    filter_metadata(metadata_path=metadata_path, filered_list=existent_and_age_range, output=metadata_path)
