#!/usr/bin/python3

import subprocess

# samples correspond to Het_1, Het_2, Imm_1, Imm_2
sra_numbers = [
    "SRR12765543", "SRR12765544", "SRR12765545", "SRR12765546", "SRR12765547", "SRR12765548", "SRR12765549", "SRR12765550", "SRR12765551", "SRR12765552", "SRR12765553", "SRR12765554", "SRR12765555", "SRR12765556", "SRR12765557", "SRR12765558", "SRR12765559", "SRR12765560", "SRR12765561", "SRR12765562", "SRR12765563", "SRR12765564", "SRR12765565", "SRR12765566", "SRR12765567", "SRR12765568", "SRR12765569", "SRR12765570", "SRR12765571", "SRR12765572", "SRR12765573", "SRR12765574", "SRR12765575", "SRR12765576", "SRR12765577", "SRR12765578", "SRR12765579", "SRR12765580", "SRR12765581", "SRR12765582", "SRR12765583", "SRR12765584", "SRR12765585", "SRR12765586", "SRR12765587", "SRR12765588", "SRR12765589", "SRR12765590", "SRR12765591", "SRR12765592", "SRR12765593", "SRR12765594", "SRR12765595", "SRR12765596", "SRR12765597", "SRR12765598", "SRR12765599", "SRR12765600", "SRR12765601", "SRR12765602", "SRR12765603", "SRR12765604", "SRR12765605", "SRR12765606", "SRR12765607", "SRR12765608", "SRR12765609", "SRR12765610", "SRR12765611", "SRR12765612", "SRR12765613", "SRR12765614", "SRR12765615", "SRR12765616", "SRR12765617", "SRR12765618", "SRR12765619", "SRR12765620", "SRR12765621", "SRR12765622", "SRR12765623", "SRR12765624", "SRR12765625", "SRR12765626", "SRR12765627", "SRR12765628", "SRR12765629", "SRR12765630", "SRR12765631", "SRR12765632", "SRR12765633", "SRR12765634", "SRR12765635", "SRR12765636", "SRR12765637", "SRR12765638", "SRR12765639", "SRR12765640", "SRR12765641", "SRR12765642", "SRR12765643", "SRR12765644", "SRR12765645", "SRR12765646", "SRR12765647", "SRR12765648", "SRR12765649", "SRR12765650", "SRR12765651", "SRR12765652"
    ]

# this will download the .sra files to ~/ncbi/public/sra/ (will create directory if not present)
# for sra_id in sra_numbers:
#     print ("Currently downloading: " + sra_id)
#     prefetch = "prefetch " + sra_id
#     print ("The command used was: " + prefetch)
#     subprocess.call(prefetch, shell=True)

# this will extract the .sra files from above into a folder named 'fastq'
for sra_id in sra_numbers:
    print ("Generating fastq for: " + sra_id)
    fastq_dump = "fasterq-dump --outdir /storage/home/swashburn30/splice_data/" + sra_id + "/fastq --split-files /storage/home/swashburn30/splice_data/" + sra_id + "/" + sra_id + ".sra"
    print ("The command used was: " + fastq_dump)
    subprocess.call(fastq_dump, shell=True)