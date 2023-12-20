#!/bin/bash

# Define variables
threads=8

source /home/ng02/miniconda3/etc/profile.d/conda.sh
conda activate flyeEnv

fastq=$1
out_dir=$2

flye --nano-hq $fastq \
--meta \
--extra-params output_gfa_before_rr=1 \
--keep-haplotypes \
--threads $threads \
--out-dir $out_dir \
--read-error 0.01
