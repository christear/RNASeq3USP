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

## Output 
The intermediate processed data and the final output are inside the output dierectory 

	- final outputs including: combined.UTR3.splicing.introns.bed and filtered.UTR3.splicing.events.tsv, could be found in the subdirectory "out"
	
		- combined.UTR3.splicing.introns.bed, a bed format file of the detcted introns within 3'UTR
		
		- filtered.UTR3.splicing.events.tsv, 3'UTR splicing events supported by at least two splicing junctions  
		
	- the detected introns are listed in subdirectory "intron"
	
	- counting of RNA-seq reads covering exon-intron bundary could be found in the sudirectory "count"
	
	- the intermediate processed data are stored in the subdirectory "processed_data"
	
	- log files are named as "*.log" and "*.elog"  

Please check the column information of filtered.UTR3.splicing.events.tsv are shown below

1: 3'UTR splicing events ID, named by genomic coordinates of introns

2: Gene symbol of the events located 

3: Types of the transcirpts where the 3'UTR splicing events come from 

columns named as SplicingLevel:* indicate the splicing level of each event in each sample

columns named as SplicingCount:* indicate the number of junctions reads supporting the splicing events in each samples    

## Prerequisites
The pipline has been test with fowlloing dependencies: 

- bash 
- bedtools 
- featureCounts 
- samtools 
- perl 
- R

The pipline is not guaranteed to work if different versions of the tools are used.
