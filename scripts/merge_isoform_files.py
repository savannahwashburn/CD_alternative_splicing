#!/storage/home/swashburn/.conda/envs/splice_project/bin/python

#this created a merged tsv file with the normalized isoform and gene information 
#extract the tpm values from each sample 
import pandas as pd
import os
from pathlib import Path

# Step 1: Define the directory containing your count files
directory_path = '/storage/home/swashburn30/splice_data/kate_data/data'

# Step 2: Get a list of count files in the directory
files = os.listdir(directory_path)# Update the file extension if needed
target_extension = '.isoforms.results'
isoform_files = [f for f in os.listdir(directory_path) if f.endswith(target_extension)]

# Step 3: Initialize an empty DataFrame to store the merged data
merged_df = pd.DataFrame()

for sample in isoform_files:
    #extract sample name from file
    store = Path(sample).stem
    store_name = store.split('.')
    sample_name = store_name[0]
    print(sample_name)
    #load in the gene.result file  with columns wanted 
    df = pd.read_csv(sample, sep='\t', usecols=['transcript_id', 'gene_id', 'TPM'])
    
    df = df.rename(columns={'TPM': sample_name}) # Rename the tpm column to the sample name
    
    # Merge the normalized count data into the result DataFrame
    if merged_df.empty:
        merged_df = df.copy()  # Initialize merged_df with the first DataFrame
    else:
        merged_df = merged_df.merge(df, on=['transcript_id', 'gene_id'], how='outer')
    
    print(merged_df.head())
    
    #merged_df.merge(df, on=['gene_id', 'transcript_id(s)'])
#return merged_df

#print(merged_df.head())

#save the merged dataframe
merged_df.to_csv('washburn_rsem.merged.transcript_tpm.tsv', sep="\t", index = False)


#check to make sure that the columns are the same in my gene_tpm.tsv file and Kate's gene_tpm.tsv file 
#reason - the samples are in a different order

#df = washburn_rsem.merged.transcript_tpm.tsv
#df2 = sem.merged.transcript_tpm.tsv
for column in df.columns:
    if df[column].equals(df2[column]):
            print(f"Contents in column '{column}' match between the two DataFrames.")
    else:
            print(f"Contents in column '{column}' do not match between the two DataFrames.")