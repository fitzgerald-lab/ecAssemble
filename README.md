# ecAssemble
eAssemble is a tool to generate long read assemblies of extra-chromosomal DNA by providing a bam file with and list of genomic regions or bins amplified.

## Installation 
There are several dependencies for that can be installed using conda 
- flye assembler 
- bedtools 
- minimap2 (required for flye)

Optional (for CVLR) 
- matplotlib 

```
## install dependencies using conda 
conda env update --file environment.yml 
```

## Test data 
Currently ecAssemble can be run using a set of regions and a mini-bam in the data folder. It is designed to to take a mod bam ( with Mm and Ml tags) but currently the mod bam lacks methylation tags. 

## Test run
ecAssemble can be run using default parameters with the bam file and a set of genomic regions as input. 
```
## activate conda environment
conda activate ecAssemble

## run with default test bam and a set of amplified genomic regions
## a fastq file and a file with 50kb bins of the genomic regions will be generated for input into CVLR 
./run_ecAssemble.py 

## run the flye assembler using the extracted fastq
./scripts/run_flye.sh <test.fq> <test> ## output folder

## optional : CVLR can be run to generate a fastq with hypomethylated reads as input into flye.
```

## Check if fastqs generated are identical


```
diff test.fq output/test.fq
```

## Post processing
Annotation of the sequences in the assembly can be done using blast as performed in the manuscript, to identify genes and genomic regions in the assembly. Annotation from the Epigenome Road map can be used as well. 

## Visualization
Visualization of the final assembly.gfa graph can be done using BandangeNG
