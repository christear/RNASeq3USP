#!/bin/bash
:<<USE
extract introns based on star jucntion file ...
...
USE

cancer=${1}
echo $cancer

for BAM in /data/11000039/csizha/TCGA_data/01_gene_expression/bam/$cancer/*.bam
do
	echo $BAM
	SAMPLEID=`echo $BAM |sed s/.*\\\/// | sed s/\\.bam//`
	READLEN=`samtools view /data/11000039/csioan/labShare/TCGA/RNA/STAR_aligned/$SAMPLEID/*.bam | head -1 | awk '{print length($10)}'`
	OVER_HANG=$((READLEN - 6))
	echo junction processing: reads length of $SAMPLEID is $READLEN
	cat /data/11000039/csioan/labShare/TCGA/RNA/STAR_aligned/$SAMPLEID/*.SJ.out.tab | awk -v OVER_HANG=$OVER_HANG '$2 > OVER_HANG'| awk -v OVER_HANG=$OVER_HANG '{OFS="\t"} {
		if($4 == 1){
			print $1,"JUNC","INTRON",$2 - OVER_HANG + 1,$2 + OVER_HANG,".","+",".","intron_id "$1":"$2"-"$3" ; ss_type "$5" ; anno_type "$6" ; sout_count "$7
		}else{
			print $1,"JUNC","INTRON",$3 - OVER_HANG + 1,$3 + OVER_HANG,".","-",".","intron_id "$1":"$2"-"$3" ; ss_type "$5" ; anno_type "$6" ; sout_count "$7
		}
	}' | grep "^[1-9XY]" >  ./intron/$cancer/$SAMPLEID.intron.left.gtf
	###
	cat /data/11000039/csioan/labShare/TCGA/RNA/STAR_aligned/$SAMPLEID/*.SJ.out.tab | awk -v OVER_HANG=$OVER_HANG '$2 > OVER_HANG'| awk -v OVER_HANG=$OVER_HANG '{OFS="\t"} {
		if($4 == 1){
			print $1,"JUNC","INTRON",$3 - OVER_HANG + 1,$3 + OVER_HANG,".","+",".","intron_id "$1":"$2"-"$3" ; ss_type "$5" ; anno_type "$6" ; sout_count "$7
		}else{
			print $1,"JUNC","INTRON",$2 - OVER_HANG + 1,$2 + OVER_HANG,".","-",".","intron_id "$1":"$2"-"$3" ; ss_type "$5" ; anno_type "$6" ; sout_count "$7
		}
	}' | grep "^[1-9XY]" > ./intron/$cancer/$SAMPLEID.intron.right.gtf
	###
	cat /data/11000039/csioan/labShare/TCGA/RNA/STAR_aligned/$SAMPLEID/*.SJ.out.tab | awk -v OVER_HANG=$OVER_HANG '$2 > OVER_HANG' | awk -v OVER_HANG=$OVER_HANG '{OFS="\t"} {
		if($4 == 1){
			print $1,"JUNC","INTRON",$2,$3,".","+",".","intron_id "$1":"$2"-"$3" ; ss_type "$5" ; anno_type "$6" ; sout_count "$7
		}else{
			print $1,"JUNC","INTRON",$2,$3,".","-",".","intron_id "$1":"$2"-"$3" ; ss_type "$5" ; anno_type "$6" ; sout_count "$7
		}
	}' | grep "^[1-9XY]"  >  ./intron/$cancer/$SAMPLEID.intron.gtf
done
