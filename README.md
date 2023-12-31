# ecAssemble
ecAssemble is a tool to generate long read assemblies of extra-chromosomal DNA, by providing a bam file with and list of genomic regions or bins amplified.

## Installation 
There are several dependencies for that can be installed using conda 
- flye <https://github.com/fenderglass/Flye> 
- bedtools 
- minimap2 (required for flye)

The installation of the packages using conda via the environment file should take < 5 mins.

Optional (for CVLR <https://github.com/EmanueleRaineri/cvlr> ) 
- matplotlib 

```
## install dependencies using conda 
conda env update --file environment.yml

## install matplotlib for CVLR
pip install matplotlib

```

## Test data 
Currently ecAssemble can be run using a set of regions and a mini-bam in the data folder. It is designed to to take a mod bam ( with Mm and Ml tags) but currently the test bam file lacks methylation tags. 


## Test run
ecAssemble can be run using default parameters with the bam file and a set of genomic regions as input. ecAssemble has been tested in a Linux environment and MacOS with intel processors. Currently ecAssemble is not compatible with Apple sillicon hardware.
```
## activate conda environment
conda activate ecAssemble

## run with default test bam and a set of amplified genomic regions
## a fastq file and a file with 50kb bins of the genomic regions will be generated for input into CVLR 
./run_ecAssemble.py 

## run the flye assembler using the extracted fastq
./scripts/run_flye.sh <test.fq> <test> ## output folder

## optional : CVLR can be run to generate a fastq with hypomethylated reads as input into flye.
./scripts/run_cvlr.sh <mod_bam> <sample_name> <region>
```

## Check if fastqs generated and assembly output are identical

```
diff test.fq output/test.fq

diff test/assembly.fasta  output/flye/assembly.fasta
```


## Post processing
Annotation of the sequences in the assembly can be done using blast as performed in the manuscript, to identify genes and genomic regions in the assembly. Annotation from the Epigenome Road map can be used as well. 

## Visualization
Visualization of the final assembly.gfa graph can be done using BandangeNG
