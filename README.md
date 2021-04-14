# RNASeq3USP
transcriptome-wide  profiling 3'UTR splicing events using conventional RNA sequencing data 


## Overview
This repository contains the pipline for transcirptome wide 3'UTR splicing events identification from the aligned RNA-seq data (bam files). Output of splicing junctions from the alighment is also required and the current version only accepts the format from STAR output (*SJ.out.tab)    


## Quick running 
To run the pipline, 6 parameters are required in the correct order, including

1, directory, include all the bam files, *SJ.out.tab should be also in the same directory and the prefix of file name should be consistent with bam files

2, gtf, could be downloaded from Gencode https://www.gencodegenes.org

3, directory of output

4, read_length, read length of the RNAseq data 

5, strand of the RNAseq data, it should be the value from 0,1,2. 0: not strand specific, 1: sense strand and 2: antisense strand 

6, read_end should be either single_end/SE or paired_end/PE  

bash run.sh ${1} ${2} ${3} ${4} ${5} ${6}


## Prerequisites
The pipline has been test with fowlloing dependencies: 

- bash 
- Bedtools 
- featureCounts 
- samtools 
- perl 
- R

The pipline is not guaranteed to work if different versions of the tools are used.
