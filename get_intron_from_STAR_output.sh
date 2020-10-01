#!/bin/bash
:<<USE
get introns by extracting splicing jucntions from STAR alignment output file ...
...
USE

data_path=${1}
echo $data_path

for BAM in $data_path/$cancer/*.bam
do
    echo $BAM
    SAMPLEID=`echo $BAM |sed s/.*\\\/// | sed s/\\.bam//`
    READLEN=`samtools view $BAM | head -1 | awk '{print length($10)}'`
    OVER_HANG=$((READLEN - 6))
    echo junction processing: reads length of $SAMPLEID is $READLEN
    cat $data_path/$SAMPLEID.SJ.out.tab | awk -v OVER_HANG=$OVER_HANG '$2 > OVER_HANG'| awk -v OVER_HANG=$OVER_HANG '{OFS="\t"} {
        if($4 == 1){
            print $1,"JUNC","INTRON",$2 - OVER_HANG + 1,$2 + OVER_HANG,".","+",".","intron_id "$1":"$2"-"$3" ; ss_type "$5" ; anno_type "$6" ; sout_count "$7
        }else{
            print $1,"JUNC","INTRON",$3 - OVER_HANG + 1,$3 + OVER_HANG,".","-",".","intron_id "$1":"$2"-"$3" ; ss_type "$5" ; anno_type "$6" ; sout_count "$7
        }
    }' | grep "^[1-9XY]" >  $SAMPLEID.intron.left.gtf
    ###
    cat $data_path/$SAMPLEID.SJ.out.tab | awk -v OVER_HANG=$OVER_HANG '$2 > OVER_HANG'| awk -v OVER_HANG=$OVER_HANG '{OFS="\t"} {
        if($4 == 1){
            print $1,"JUNC","INTRON",$3 - OVER_HANG + 1,$3 + OVER_HANG,".","+",".","intron_id "$1":"$2"-"$3" ; ss_type "$5" ; anno_type "$6" ; sout_count "$7
        }else{
            print $1,"JUNC","INTRON",$2 - OVER_HANG + 1,$2 + OVER_HANG,".","-",".","intron_id "$1":"$2"-"$3" ; ss_type "$5" ; anno_type "$6" ; sout_count "$7
        }
    }' | grep "^[1-9XY]" > $SAMPLEID.intron.left.gtf
    ###
    cat $data_path/$SAMPLEID.SJ.out.tab | awk -v OVER_HANG=$OVER_HANG '$2 > OVER_HANG' | awk -v OVER_HANG=$OVER_HANG '{OFS="\t"} {
        if($4 == 1){
            print $1,"JUNC","INTRON",$2,$3,".","+",".","intron_id "$1":"$2"-"$3" ; ss_type "$5" ; anno_type "$6" ; sout_count "$7
        }else{
            print $1,"JUNC","INTRON",$2,$3,".","-",".","intron_id "$1":"$2"-"$3" ; ss_type "$5" ; anno_type "$6" ; sout_count "$7
        }
    }' | grep "^[1-9XY]"  >  $SAMPLEID.intron.gtf
done

