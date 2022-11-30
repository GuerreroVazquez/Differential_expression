import gseapy
import csv
import numpy as np
import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt

import Analysis_code.analyze_de_genes as a_g
import sys

orig_stdout = sys.stdout
default_pathway_file = '../../Results/Pathway_A_def_expressed_gene_names.csv'
databases_names = gseapy.get_library_name()



def get_pathways(genes_file='../../Results/def_expressed.csv'):

    df, names = a_g.get_data(all_experiments_de=genes_file)
    df_results = {}
    df = a_g.remove_version(genes_data_frame=df,experiment_names=names,output=None)
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
                print(f"Error",e)
        df_results[experiment] = experiment_dic
    return df_results

def get_pathways(genes_file='../../Results/def_expressed_gene_names.csv'):

    df, names = a_g.get_data(all_experiments_de=genes_file)
    df_results = pd.DataFrame()
    df = a_g.remove_version(genes_data_frame=df,experiment_names=names,output=None)

    databases_names=['GO_Molecular_Function_2021',
                     'GO_Cellular_Component_2021',
                     'GO_Biological_Process_2021',
                     'Reactome_2022',
                     'KEGG_2021_Human']

    for experiment in names:
        DEGs = df[experiment].tolist()
        DEGs = [x for x in DEGs if pd.notnull(x)]
        experiment_dic = pd.DataFrame()
        for database in databases_names:
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
                print(f"Error",e)
        experiment_dic["Experiment"]=experiment
        df_results = pd.concat([df_results, experiment_dic])
    return df_results


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


def chicharron():
    """
    This funtion will get the csv of the patways generated previously and ... we will see
    :return:
    """
    df = pd.read_csv(default_pathway_file)
    pvalue_threshold = 0.1
    apvlue_threshold = 0.1
    ratio = [int(x.split('/')[0])/int(x.split('/')[1]) for x in df['Overlap']]
    df['Ratio']= ratio
    df = df.query(f"P-value < {pvalue_threshold}")
    print(df)

def test_chicharron():
    chicharron()