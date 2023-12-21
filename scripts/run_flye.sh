#!/bin/bash

# Define variables
threads=8

fastq=$1
out_dir=$2

flye --nano-hq $fastq \
--meta \
--extra-params output_gfa_before_rr=1 \
--keep-haplotypes \
--threads $threads \
--out-dir $out_dir \
--read-error 0.01 \
--deterministic


