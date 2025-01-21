#!/storage/home/swashburn/.conda/envs/splice_project/bin/python
#This program calculate the isoform fraction of each isoform
#IF = normalized isoform counts / normalized total gene counts

#code adapted from Kate Xie

import argparse
import pandas as pd
import numpy as np
np.seterr(divide='ignore', invalid='ignore')



def main():

    #argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene_count",action='store', metavar="head_gene_tpm.tsv", type=str, required=True, help='a tab seperated normalized gene count file')
    parser.add_argument("--isoform_count",action='store', metavar='head_isoform_tpm.tsv', type=str, required=True, help='a tab seperated normalized isoform count file')

    args = parser.parse_args()


    #processing each gene
    df_gene = pd.read_csv(args.gene_count, sep='\t', header=0)
    df_isoform = pd.read_csv(args.isoform_count, sep='\t', header=0)



    df_IF = all_manipulate(df_gene,df_isoform)
    df_IF.to_csv('isoform_fraction.tsv', sep='\t', header=True, index=False, na_rep='NaN')


def all_manipulate(df_gene: pd.DataFrame, df_isoform: pd.DataFrame) -> pd.DataFrame:
    #check if the sample names in two data frame are the same
    if list(df_gene.columns[2:]) != list(df_isoform.columns[2:]):
        exit("The sample names in isoform data and gene count data are different")

    sample_size = df_isoform.shape[1]-3  #initialize sample size
    df_IF = pd.DataFrame(columns=df_isoform.columns)   #initialize the isoform fraction data frame


    for index in range(len(df_isoform)):
        row_isoform = df_isoform.iloc[[index]] #extract row from df_isoform
        row_isoform = row_isoform.values.flatten().tolist()
        transcript_id = row_isoform[0] #extract transcript id
        gene_id = row_isoform[1]  #extract gene_id


        #check if the gene id is also in gene count data frame
        if gene_id in df_gene.loc[:, 'gene_id'].values.flatten().tolist():
            row_gene = df_gene.loc[df_gene['gene_id'] == gene_id,:].values.flatten().tolist()
            row_isoform = np.array(row_isoform[2:])
            row_gene = np.array(row_gene[2:])
            IF = list(np.divide(row_isoform, row_gene))  #divide the isoform counts / gene counts
            row_IF = [transcript_id, gene_id] + IF  #concatenate with the transcript id and the gene id
        else:
            row_IF = [transcript_id, gene_id] + ['NaN']*sample_size  #concatenate list

        df_IF.loc[len(df_IF)] = row_IF



    return df_IF






if __name__ == "__main__":
    main()
