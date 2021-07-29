#!/bin/bash
:<<USE
USE

data_path=${1}
read_length=${2}
out_dir=${3}
### if_stranded should be 0,1 or 2, which indicates the strandness of sequencing data  
# 0 means the data is not strand specific, 1 means sense strand and 2 means antisense strand 
strandness=${4}
readend=${5}

echo counting $data_path 
for BAM in $data_path/*.bam
do
	echo $BAM
	SAMPLEID=`echo $BAM |sed s/.*\\\/// | sed s/\\.bam//`
	if [[ $readend =~ [Pp].*[Ee] ]]; then
		echo paired end 
		echo left SS
	    featureCounts -a $out_dir/combined.intron.left.gtf.utr3 -o $out_dir/$SAMPLEID.intron.left $BAM -s $strandness --minOverlap $read_length -p -F GTF -t INTRON -g intron_id -T 20 -B -C -O
		echo right SS
	    featureCounts -a $out_dir/combined.intron.right.gtf.utr3 -o $out_dir/$SAMPLEID.intron.right $BAM -s $strandness --minOverlap $read_length -p -F GTF -t INTRON -g intron_id -T 20 -B -C -O
	else
		echo single end
		echo left SS
	    featureCounts -a $out_dir/combined.intron.left.gtf.utr3 -o $out_dir/$SAMPLEID.intron.left $BAM -s $strandness --minOverlap $read_length -F GTF -t INTRON -g intron_id -T 20 -B -C -O
		echo right SS
	    featureCounts -a $out_dir/combined.intron.right.gtf.utr3 -o $out_dir/$SAMPLEID.intron.right $BAM -s $strandness --minOverlap $read_length -F GTF -t INTRON -g intron_id -T 20 -B -C -O 
	fi
done

