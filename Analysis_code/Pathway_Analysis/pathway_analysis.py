import gseapy
import csv
import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt

import Analysis_code.analyze_de_genes as a_g
import sys
from wordcloud import WordCloud, STOPWORDS
from gseapy.plot import barplot, dotplot

orig_stdout = sys.stdout
default_pathway_file = '../../Results/Pathway_A_def_expressed_gene_names.csv'
databases_names = gseapy.get_library_name()


def get_pathways(genes_file='../../Results/def_expressed.csv'):
    """
    From the gene_files csv, gets ALL the pathways in ALL the databases of ALL the experiments
    :param genes_file: CSV file with the experiments df genes in each column.
    :return:
    """
    df, names = a_g.get_data(all_experiments_de=genes_file)
    df_results = {}
    df = a_g.remove_version(genes_data_frame=df, experiment_names=names, output=None)
    for experiment in names:
        DEGs = df[experiment].tolist()
        DEGs = [x for x in DEGs if pd.notnull(x)]
        experiment_dic = {}
        for database in databases_names:
            try:
                enr_GOBP_up = gp.enrichr(gene_list=DEGs,
                                         gene_sets=[database],
                                         organism='Human',
                                         outdir='test/enr_DEGs_GOBP_up',
                                         cutoff=0.5
                                         )
                if not enr_GOBP_up.results.empty:
                    experiment_dic[database] = enr_GOBP_up.results
            except Exception as e:
                print(f"Error", e)
        df_results[experiment] = experiment_dic
    barplot(df, title=experiment, color='r')
    return df_results


def get_pathways(genes_file='../../Results/def_expressed_gene_names.csv', sel_databases_names=None):
    """
     From the gene_files csv, gets the pathways in  the databases in databases_names
      of ALL the experiments
     :param sel_databases_names: List of string with the subset of gseapy.get_library_name()
     :param genes_file: CSV file with the experiments df genes in each column.
     :return:
     """
    if sel_databases_names is None:
        sel_databases_names = ['GO_Molecular_Function_2021',
                               'GO_Cellular_Component_2021',
                               'GO_Biological_Process_2021',
                               'Reactome_2022',
                               'KEGG_2021_Human']
    df, names = a_g.get_data(all_experiments_de=genes_file)
    df_results = pd.DataFrame()
    df = a_g.remove_version(genes_data_frame=df, experiment_names=names, output=None)

    for experiment in names:
        DEGs = df[experiment].tolist()
        DEGs = [x for x in DEGs if pd.notnull(x)]
        experiment_dic = pd.DataFrame()
        for database in sel_databases_names:
            try:
                enr_GOBP_up = gp.enrichr(gene_list=DEGs,
                                         gene_sets=[database],
                                         organism='Human',
                                         outdir=f'test/{experiment}',
                                         cutoff=0.5
                                         )
                if not enr_GOBP_up.results.empty:
                    experiment_dic = pd.concat([experiment_dic, enr_GOBP_up.results])
            except Exception as e:
                print(f"Error", e)
        experiment_dic["Experiment"] = experiment
        df_results = pd.concat([df_results, experiment_dic])
    df_results.rename(columns={'P-value': 'PValue'}, inplace=True)
    df_results.rename(columns={'Adjusted P-value': 'APValue'}, inplace=True)

    return df_results


def get_wordcloud(df, title):
    comment_words = ''
    stopwords = set(STOPWORDS)

    # iterate through the csv file
    for val in df.Term:

        # typecaste each val to string
        val = str(val)

        # split the value
        tokens = val.split()

        # Converts each token into lowercase
        for i in range(len(tokens)):
            tokens[i] = tokens[i].lower()

        comment_words += " ".join(tokens) + " "
    comment_words = comment_words.replace("go", "")
    comment_words = comment_words.replace("hsa", "")
    wordcloud = WordCloud(width=800, height=800,
                          background_color='white',
                          stopwords=stopwords,
                          min_word_length=2,
                          min_font_size=10).generate(comment_words)
    # plot the WordCloud image
    if len(wordcloud.words_) > 1:
        plt.figure(figsize=(8, 8), facecolor=None)
        plt.imshow(wordcloud)
        plt.axis("off")
        plt.tight_layout(pad=0)
        plt.title(title)

        plt.show()
        plt.savefig(f"Wordcloud_{title}")
    return wordcloud.words_


def plot_pathways(df, experiments=None):
    """
    Gets the dataframe with all the pawthays
     data from all the experiments and makes the ploting by experiment
    :param df: output of gseaspy
    :return:
    """
    if experiments is None:
        experiments = ["GSE152558", "GSE157585", "GSE164471",
                       "GSE164471_MvYoung", "GSE164471_MvOld"]
    by_experiment = []
    for experiment in experiments:
        by_experiment.append(df.query(f"Experiment == '{experiment}'"))
    for i, df in enumerate(by_experiment):
        barplot(df, title=experiments[i], cutoff=0.05, color='r')
        barplot(df, column='invRatio', title=experiments[i], cutoff=-0.033, color='b')


def test_get_pathways():
    """

    :return:
    """
    f = open('out.txt', 'w')
    sys.stdout = f
    result = get_pathways()
    sys.stdout = orig_stdout
    f.close()
    result.to_csv(default_pathway_file, index=False)
    print(result)


def filter_pathways(df=None, diff_expression=default_pathway_file, pvalue_threshold=0.1, apvalue_threshold=1,
                    ratio_treshold=0):
    """
    This function will get the csv of the pathways generated previously and get the patways that
    have a criteria of pvalue, adjusted pvalue and ratio
    :param df: dataframe The dataframe output from the gseapy
    :param apvalue_threshold: float Threshold that marks up to which value to take the adjusted P value
    :param pvalue_threshold:  float Threshold that marks up to which value to take the P value
    :param ratio_treshold: float What is the minimun ratio (genes found/ genes total in pathway)
    :param diff_expression: str the CSV file that holds the output of gseapy

    :return: dataframe with the filters applied
    """
    if df is None and diff_expression is not None:
        df = pd.read_csv(diff_expression)

    ratio = [int(x.split('/')[0]) / int(x.split('/')[1]) for x in df['Overlap']]
    df['Ratio'] = ratio
    df = df.query(f"PValue < {pvalue_threshold}")
    df = df.query(f"APValue < {apvalue_threshold}")
    df = df.query(f"Ratio > {ratio_treshold}")
    return df


def test_chicharron():
    filter_pathways()


def test_get_wordcloud():
    df = filter_pathways()
    by_experiment = []
    experiments = ["GSE152558", "GSE157585", "GSE164471",
                   "GSE164471_MvYoung", "GSE164471_MvOld"]
    experiment = experiments[3]
    df = df.query(f"Experiment == '{experiment}'")
    words = get_wordcloud(df, experiment)


def test_dec5_getting_pathways_up_down():
    """
    This is ran in order to get the same patways and wordclouds but separating on up and down
    regulated
    :return:
    """

    #up_down_file = '../../Results/up_regulated_all.csv'
    #pathways = get_pathways(genes_file=up_down_file)
    #pathways.to_csv('../../Results/Pathway_A_def_expressed_gene_names_UP_DOWN.csv', index=False)
    pathways = filter_pathways(df='../../Results/Pathway_A_def_expressed_gene_names_UP_DOWN.csv', pvalue_threshold=0.1, apvalue_threshold=0.1, ratio_treshold=0)
    experiments = ["GSE152558", "GSE157585", "GSE164471",
                   "GSE164471_MvYoung", "GSE164471_MvOld"]
    for experiment in experiments:
        df = df.query(f"Experiment == '{experiment}'")
        words = get_wordcloud(df, experiment)


def test_dec9_Tony_pathways():
    """
    This is ran in order to get the same patways and wordclouds but separating on up and down
    regulated
    :return:
    """

    up_down_file = '../../Tony.csv'
    pathways = get_pathways(genes_file=up_down_file, sel_databases_names=databases_names)
    pathways.to_csv('../../Results/Pathway_A_def_expressed_gene_names_Tony_all_db.csv', index=False)

def test_Feb_3_microarray_1_pathways():
    """
    This is ran in order to get the same patways and wordclouds but separating on up and down
    regulated
    :return:
    """

    up_down_file = '../../GSE87105.csv'
    sel_databases_names = ['GO_Molecular_Function_2021',
                           'GO_Cellular_Component_2021',
                           'GO_Biological_Process_2021',
                           'Reactome_2022',
                           'KEGG_2021_Human']
    pathways = get_pathways(genes_file=up_down_file, sel_databases_names=sel_databases_names)
    pathways.to_csv('../../Results/Pathway_A_def_expressed_gene_names_GSE87105_all_db.csv', index=False)