#!/bin/bash
:<<USE
get introns by extracting splicing jucntions from STAR alignment output file ...
...
USE

data_path=${1}
read_length=${2}
out_dir=${3}
echo processing $data_path $read_length
if [ -f $out_dir/sampleID.list ]; then
	rm -f $out_dir/sampleID.list
fi
touch $out_dir/sampleID.list

for BAM in $data_path/$cancer/*.bam
do
    echo $BAM
    SAMPLEID=`echo $BAM |sed s/.*\\\/// | sed s/\\.bam//`
	echo $SAMPLEID >> $out_dir/sampleID.list
    OVER_HANG=$((read_length - 6))
    cat $data_path/$SAMPLEID.SJ.out.tab | awk -v OVER_HANG=$OVER_HANG '$2 > OVER_HANG'| awk -v OVER_HANG=$OVER_HANG '{OFS="\t"} {
        if($4 == 1){
            print $1,"JUNC","INTRON",$2 - OVER_HANG + 1,$2 + OVER_HANG,".","+",".","intron_id "$1":"$2"-"$3" ; ss_type "$5" ; anno_type "$6" ; sout_count "$7
        }else{
            print $1,"JUNC","INTRON",$3 - OVER_HANG + 1,$3 + OVER_HANG,".","-",".","intron_id "$1":"$2"-"$3" ; ss_type "$5" ; anno_type "$6" ; sout_count "$7
        }
    }' | grep "^[1-9XY]" >  $out_dir/$SAMPLEID.intron.left.gtf
    ###
    cat $data_path/$SAMPLEID.SJ.out.tab | awk -v OVER_HANG=$OVER_HANG '$2 > OVER_HANG'| awk -v OVER_HANG=$OVER_HANG '{OFS="\t"} {
        if($4 == 1){
            print $1,"JUNC","INTRON",$3 - OVER_HANG + 1,$3 + OVER_HANG,".","+",".","intron_id "$1":"$2"-"$3" ; ss_type "$5" ; anno_type "$6" ; sout_count "$7
        }else{
            print $1,"JUNC","INTRON",$2 - OVER_HANG + 1,$2 + OVER_HANG,".","-",".","intron_id "$1":"$2"-"$3" ; ss_type "$5" ; anno_type "$6" ; sout_count "$7
        }
    }' | grep "^[1-9XY]" >  $out_dir/$SAMPLEID.intron.left.gtf
    ###
    cat $data_path/$SAMPLEID.SJ.out.tab | awk -v OVER_HANG=$OVER_HANG '$2 > OVER_HANG' | awk -v OVER_HANG=$OVER_HANG '{OFS="\t"} {
        if($4 == 1){
            print $1,"JUNC","INTRON",$2,$3,".","+",".","intron_id "$1":"$2"-"$3" ; ss_type "$5" ; anno_type "$6" ; sout_count "$7
        }else{
            print $1,"JUNC","INTRON",$2,$3,".","-",".","intron_id "$1":"$2"-"$3" ; ss_type "$5" ; anno_type "$6" ; sout_count "$7
        }
    }' | grep "^[1-9XY]"  >  $out_dir/$SAMPLEID.intron.gtf
done

