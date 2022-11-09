import csv


import pandas as pd



def get_union():
    pass



def main():

    df = pd.read_csv('Results/DSeq/All_experiments_dE.csv')
    experiment_names = df.head()
    list(set(df[0]) & set(df[1]))
