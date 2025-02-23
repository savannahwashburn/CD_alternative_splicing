# Post-operative ileal RNA-seq 

This project examines the role of alternative splicing in recurring Crohn's disease patients. Post-operative colonscopies were collected and sequenced for bulk RNA-seq and genotype data was collected. 

Additional publicly available data of ileum and colon bulk RNA-seq was obtained through GEO accession number GSE168952.

## Methods

The raw data was (post-operative biopsies and publicly available data) processed with nf-core/rnaseq (Nextflow pipeline) to quantify gene expression and isoform counts.

## Isoform Fractions

Isoform fractions were calculated as normalized isoform counts / normalized total gene counts 
* filtered if gene is not expression
* filtered if >95% of individuals express isoform as the only isoform of this gene 
* filtered if >95% of individuals do not express transcript 

## Differentially used Isoform Fractions (dIF)
mean IF in recurring disease - mean IF in non-recurring disease for each isoform
or 
mean IF in rectum - mean IF in ileum for each isoform 
* considered differentially used if dIF > +0.1 or < -0.1 with FDR adjusted p-value < 0.05

## Scripts directory

Contains scripts used to create sample sheet for publicly available data, run nf-core/rnaseq, and analysis of the post-op data. 
