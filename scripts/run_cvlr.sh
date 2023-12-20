#!/bin/bash

modBam=$1
sampleName=$2
region=$3
regionName=$(echo $region | sed -E 's/:|\-/_/g')
sampleName=${sampleName}-${regionName}

echo $region

./bin/cvlr-meth-of-bam ${modBam} ${region} > ${sampleName}-matrix.txt
./bin/cvlr-cluster ${sampleName}-matrix.txt 2 1 100 > ${sampleName}-clusters.txt 
./bin/cvlr-stats.py cluster ${sampleName}-clusters.txt ${sampleName}-matrix.txt > ${sampleName}-stats.txt

echo "Ended run"
date
