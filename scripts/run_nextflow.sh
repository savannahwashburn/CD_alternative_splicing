#!/bin/bash

#activate conda env

source activate java

PROJ_DIR=/storage/home/swashburn30/splice_data/nextflow_test
NFC_PIPE=nf-core/rnaseq
NFC_VER=3.5
NFC_PROFILE=singularity
SAMPLESHEET=$PROJ_DIR/samplesheet_nextflow_test.csv
OUTDIR=$PROJ_DIR/nextflow_output_test
GENOME=GRCh38
ALIGNER=star_rsem

cd $PROJ_DIR

NXF_VER=22.10.6 nextflow run $NFC_PIPE \
	-r $NFC_VER \
	-profile $NFC_PROFILE \
        --input $SAMPLESHEET \
        --outdir $OUTDIR \
        --genome $GENOME \
        --aligner $ALIGNER


