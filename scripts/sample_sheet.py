import pandas as pd

file = open("/Users/swashburn30/Desktop/isoform/SRR_Acc_List.txt", "r")
data = file.read() 
#data_into_list = data.replace('\n', '""').split(",") 
data_into_list = data.split("\n") 
del data_into_list[-1] #remove last element - was empty 

data_into_list_1 = ["/storage/home/hcoda1/5/swashburn30/scratch/splice_project/kiera_fastq_files/" + x + "_1.fastq.gz" for x in data_into_list]
data_into_list_2 = ["/storage/home/hcoda1/5/swashburn30/scratch/splice_project/kiera_fastq_files/" + x + "_2.fastq.gz" for x in data_into_list]

#read in sample names 
file_sample = open("/Users/swashburn30/Desktop/isoform/sample_name.txt", "r")
data_sample = file_sample.read()
sample_list = data_sample.split("\n")
del sample_list[-1]

#strand list 
strand = ['reverse']
strand_list = strand * 119

#dictionary with sample information
data = {'sample': sample_list,
        'fastq_1': data_into_list_1,
        'fastq_2': data_into_list_2,
        'strandedness': strand_list}

df = pd.DataFrame(data)

#save the dataframe 
df.to_csv('samplesheet_1_30_24.csv', encoding='utf-8', index=False)

