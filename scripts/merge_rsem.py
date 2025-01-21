#!/storage/home/swashburn/.conda/envs/splice_project/bin/python

#adapted from Kate Xie

import pandas as pd
import argparse
import os




dir = '/storage/home/swashburn30/splice_data/kate_data/data'
abundance_col = 'FPKM'

#check if the directory ends with back slash or not
if dir.endswith('/'):
    dir = dir[:-1]



samples = os.listdir(dir)
merged_df =pd.read_csv() #initialize the output data frame


for samples in os.listdir(dir):
    isoform = f'{dir}/{samples}/{samples}.isoforms.results'
    gene = f'{dir}/{samples}/{samples}.genes.results'
  





