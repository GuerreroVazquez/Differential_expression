import csv
import math

import pandas as pd

from Analysis_code.ensmbl_finder import get_gene_name_from_ensmbl


def test_get_gene_name_from_ensmbl():
    import csv
    with open('/home/karen/Downloads/LongHack_datasets/GSE164471_txi.csv', 'r') as csvfile:
        reader = csv.reader(csvfile)
        next(reader)
        genenames=[]# skip header row
        counter = 0
        for row in reader:
            ensmbl = row[0]
            ensmbl= ensmbl.split(".")[0]
            symbol = get_gene_name_from_ensmbl(ensmbl)
            genenames.append(ensmbl+","+symbol)
            counter += 1
            if counter>100:
                with open("/home/karen/Downloads/LongHack_datasets/genenames_71.txt", "a") as f:
                    for name in genenames:

                        f.write("%s\n" % name)
                genenames=[]
                counter=0


def test_get_logvalue():
    files = ["Age.csv", "bedding.csv","ULLS.csv"]
    for file in files:
        file = "/home/karen/Downloads/LongHack_datasets/Ready/"+file
        with open(file, "r") as data_file:
            reader = csv.reader(data_file)
            modified=pd.DataFrame()
            sum = 0
            averages=[]
            for row in reader:
                if row[0]=="GeneSymbol":
                    modified=pd.DataFrame(columns = row)
                else:
                    c=0
                    for column in row:

                        if column.replace('.', '', 1).isdigit():
                            sum += float(column)
                            c=c+1
                    if c<1:
                        print("OH NOOO")
                    average = sum/c
                    new_value=[]
                    for column in row:
                        if column.replace('.', '', 1).isdigit():
                            divition =float(column)/average
                            if divition==0:
                                new_value.append("NA")
                            else:
                                new_value.append(math.log2(divition))

                        else:
                            new_value.append(column)
                    modified.loc[len(modified)] = new_value
            modified.to_csv(file[:-3]+"_Lfc.csv")

