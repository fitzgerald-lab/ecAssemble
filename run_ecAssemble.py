#!/usr/bin/env python3
import os
import re
import pysam
import logging
import argparse
import subprocess
import pandas as pd 
from itertools import product
from multiprocessing import Pool


def extract_reads(bam, region_list, sample_name, read_file=False):
    '''Extract reads from a BAM based on a genomic region and write to a FASTQ file'''
    output_fastq = f'{sample_name}.fq'
    read_names = []

    with pysam.AlignmentFile(bam, 'rb') as bam_file:
        if read_file:
            with open(read_file, 'r') as read_buf:
                target_reads = [i.rstrip('\n') for i in read_buf.readlines()]

        for region in region_list.split(','):
            chrom, start, end = re.split(':|-', region)
            reads = bam_file.fetch(chrom, int(start), int(end))

            with open(output_fastq, "w") as fastq_out:
                for read in reads:
                    try:
                        if read_file:
                            if (read.query_name not in read_names) and (read.query_name in target_reads):
                                fastq_out.write(
                                    f"@{read.query_name}\n{read.seq}\n+\n{''.join(map(lambda x: chr(x+33), read.query_qualities))}\n")
                                read_names.append(read.query_name)
                        else:
                            if read.query_name not in read_names:
                                fastq_out.write(
                                    f"@{read.query_name}\n{read.seq}\n+\n{''.join(map(lambda x: chr(x+33), read.query_qualities))}\n")
                                read_names.append(read.query_name)
                    except Exception as e:
                        logging.error('ERROR : ' + str(e))


def generate_windows(region_bed, interval_size):
    '''Generate a BED file of genome windows for parallel CVLR'''
    bed_out = f'{region_bed.rstrip(".bed")}{int(interval_size/1000)}kb.bed'
    cmd = ['bedtools', 'makewindows', '-b', region_bed, '-w', str(interval_size)]
    bed_run = subprocess.run(cmd, shell=False, check=True, capture_output=True, text=True)

    with open(bed_out, 'w') as bed_file:
        bed_file.write(bed_run.stdout)

    format_cmd = ['awk', '{print$1 "\t" $2 "\t" $3 "\t" $1 ":" $2 "-" $3}', bed_out]
    awk_run = subprocess.run(format_cmd, shell=False, check=True, capture_output=True, text=True)

    with open(bed_out, 'w') as bed_file:
        bed_file.write(awk_run.stdout)

            
def run_cvlr(mod_bam, sample_name, region):
    '''Takes a modbam and a 10-50kb window to carry out methylation clustering'''
    cvlr_cmd = ['../scripts/runCVLR.sh', mod_bam, sample_name, region]
    cvlr_run = subprocess.run(cvlr_cmd, shell = False, check= True, capture_output=True, text=True)


def parse_cvlr(sample_name):
    '''Parse CVLR output and get a list of readnames of hypomethylated reads'''
    stat_files = os.listdir()
    stat_files =[i for i in stat_files if 'stat' in i]

    read_interest = [] 
    empty = [] 
    for s in stat_files: 
        try : 
            stats = pd.read_csv(s, sep ='\t', comment = '#', header=None) 
            if stats[3].mean() < stats[6].mean(): 
                group = 0 # group we want  
            else: 
                group = 1  
                clust = pd.read_csv(s.split('stats.txt')[0]+'clusters.txt', sep ='\t', comment = '#', header=None) 
                read_interest.append((clust.loc[clust[1] == group, 0 ])) 
        except Exception as e: 
            print(e) 
            empty.append(s) 
    
    readnames = pd.concat(read_interest).drop_duplicates()
    readnames.to_csv('%s_hypo-readnames.tsv' % sample_name , sep ='\t', header = False, index = False) 
    
input_bed     = 'data/regions.bed'
sample_name   = 'test'
interval_size = 50000
region_list   = 'chr7:116235772-131171604'
bam           = 'data/test.bam'
mod_bam       = 'data/test.bam'
threads       = 2

parser = argparse.ArgumentParser(description='Process some parameters.')

# Define arguments
parser.add_argument('--input_bed', type=str, default=input_bed, help='Input BED file')
parser.add_argument('--sample_name', type=str, default=sample_name, help='Sample name')
parser.add_argument('--interval_size', type=int, default=interval_size, help='Interval size')
parser.add_argument('--region_list', type=str, default=region_list, help='Region list')
parser.add_argument('--bam', type=str, default=bam, help='BAM file path')
parser.add_argument('--mod_bam', type=str, default=mod_bam, help='Modified BAM file path')
parser.add_argument('--threads', type=int, default=threads, help='Number of threads')

# Parse arguments from the command line
args = parser.parse_args()

# Assign parsed arguments to variables
input_bed = args.input_bed
sample_name = args.sample_name
interval_size = args.interval_size
region_list = args.region_list
bam = args.bam
mod_bam = args.mod_bam
threads = args.threads

# Rest of your script using these variables
print(f"Input BED: {input_bed}")
print(f"Sample Name: {sample_name}")
print(f"Interval Size: {interval_size}")
print(f"Region List: {region_list}")
print(f"BAM: {bam}")
print(f"Modified BAM: {mod_bam}")
print(f"Threads: {threads}")


def main(bam, region_list, sample_name, input_bed, interval_size, mod_bam = bam, threads = 2):
    window_bed = input_bed.rstrip('.bed') + str(int(int(interval_size)/1000)) + 'kb.bed'
    extract_reads(bam, region_list, sample_name)
    generate_windows(input_bed, interval_size)

    with open(window_bed, 'r') as bed_file:
        regions = bed_file.readlines()
        window_list = [i.split('\t')[-1].rstrip('\n') for i in regions]

    ## run cvlr if there is a mod bam for read clustering using methylation marks
    #with Pool(threads) as p:
    #    p.starmap(run_cvlr, product([mod_bam], [sample_name], window_list))
    #parse_cvlr(sample_name)
    hypo_file  = '%s_hypo-readnames.tsv' % sample_name
    #extract_reads(bam, region_list, sample_name + '_hypo', read_file=hypo_file)


    
if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
    main( bam, region_list, sample_name, input_bed, interval_size, mod_bam, threads)
