#!/bin/bash
:<<USE
3'UTR splicing events identification ... 
USE


function usage () {
        echo
        echo
        echo Usage: bash run.sh [directory of bam] [gtf] [output directory] [RNAseq read_length] [RNAseq strand] [RNAseq reads_end]
        echo "      " directory of bam file should also include *SJ.out.tab files with identifcal prefix as bam files   
        echo "      " gtf could be downloaded from Gencode https://www.gencodegenes.org
		echo "      " output directory
		echo "      " read_length of RNAseq data
        echo "      " strand should be the value from 0,1,2. 0: not strand specific, 1: sense strand and 2: antisense strand 
        echo "      " read_end should be either single_end/SE or paired_end/PE 
        echo
		echo
}

if [ -z $6 ]; then
        usage
        exit
fi


data_path=${1}
gtf=${2}
out_dir=${3}
read_length=${4}
strandness=${5}
readend=${6}

echo processing annotation 
perl ./perl/get.utr.from.annotation.v2.pl UTR $gtf $out_dir/annotated.utr 1>$out_dir/extract.utr.log 2>$out_dir/extract.utr.elog
#perl ./perl/get.utr.from.annotation.v2.pl cds $gtf $out_dir/annotation.cds
awk '$3 == "UTR3"' $out_dir/annotated.utr > $out_dir/annotated.utr3
awk '$3 == "UTR5"' $out_dir/annotated.utr > $out_dir/annotated.utr5
awk '$3 == "CDS"' $gtf > $out_dir/annotated.CDS

echo extracting introns from bam file ...
./get_intron_from_STAR_output.sh $data_path $read_length $out_dir 1>$out_dir/extract.intron.log 2>$out_dir/extract.intron.elog
if [ -f ./$out_dir/intron ]; then
	mkdir $out_dir/intron
fi
mv $out_dir/*.gtf $out_dir/intron/

echo combining splicing events together ...
perl ./perl/combine.events.pl $out_dir/combined.intron.left.gtf $out_dir/intron/*.left.gtf > combined.splicing.out.fromleft.txt 
perl ./perl/combine.events.pl $out_dir/combined.intron.right.gtf $out_dir/intron/*.right.gtf > combined.splicing.out.fromright.txt
perl ./perl/combine.events.pl $out_dir/combined.intron.gtf $out_dir/intron/*.intron.gtf > combined.splicing.out.txt

echo extracting intron from UTR3
intersectBed -a $out_dir/combined.intron.left.gtf -b $out_dir/annotated.utr3 -f .5 -s -wo | awk '{OFS="\t"} {print $10,$13,$16,$19,$22,$29,substr(substr($0,index($0,"gene_name")),10,index(substr($0,index($0,"gene_name")),";") - 10),substr(substr($0,index($0,"transcript_type")),16,index(substr($0,index($0,"transcript_type")),";") - 16),$20":"$23"-"$24,$7}' | sed s/[\"\;]//g | sort | uniq  > $out_dir/combined.intron.left.gtf.utr3.wo
intersectBed -a $out_dir/combined.intron.left.gtf -b $out_dir/annotated.utr3 -f .5 -s -u > $out_dir/combined.intron.left.gtf.utr3
intersectBed -a $out_dir/combined.intron.right.gtf -b $out_dir/annotated.utr3 -f .5 -s -wo | awk '{OFS="\t"} {print $10,$13,$16,$19,$22,$29,substr(substr($0,index($0,"gene_name")),10,index(substr($0,index($0,"gene_name")),";") - 10),substr(substr($0,index($0,"transcript_type")),16,index(substr($0,index($0,"transcript_type")),";") - 16),$20":"$23"-"$24,$7}' | sed s/[\"\;]//g | sort | uniq  > $out_dir/combined.intron.right.gtf.utr3.wo
intersectBed -a $out_dir/combined.intron.right.gtf -b $out_dir/annotated.utr3 -f .5 -s -u > $out_dir/combined.intron.right.gtf.utr3
intersectBed -a $out_dir/combined.intron.gtf -b $out_dir/annotated.utr3 -f 1 -s -wo | awk '{OFS="\t"} {print $10,$13,$16,$19,$22,$29,substr(substr($0,index($0,"gene_name")),10,index(substr($0,index($0,"gene_name")),";") - 10),substr(substr($0,index($0,"transcript_type")),16,index(substr($0,index($0,"transcript_type")),";") - 16),$20":"$23"-"$24,$7}' | sed s/[\"\;]//g | sort | uniq  > $out_dir/combined.intron.gtf.utr3.wo
intersectBed -a $out_dir/combined.intron.gtf -b $out_dir/annotated.utr3 -f 1 -s -u > $out_dir/combined.intron.gtf.utr3
intersectBed -a $out_dir/combined.intron.gtf -b $out_dir/annotated.CDS -f 1 -s -wo | awk '{OFS="\t"} {print $10,$13,$16,$19,$22,$29,substr(substr($0,index($0,"gene_name")),10,index(substr($0,index($0,"gene_name")),";") - 10),substr(substr($0,index($0,"exon_id")),8,index(substr($0,index($0,"exon_id")),";") - 8),$20":"$23"-"$24,$7}' | sed s/[\"\;]//g | sort | uniq  > $out_dir/combined.intron.gtf.cds.wo
intersectBed -a $out_dir/combined.intron.gtf -b $out_dir/annotated.utr5 -f 0.1 -s -wo | awk '{OFS="\t"} {print $10,$13,$16,$23,$26,$33,substr(substr($0,index($0,"gene_name")),10,index(substr($0,index($0,"gene_name")),";") - 10),substr(substr($0,index($0,"exon_id")),8,index(substr($0,index($0,"exon_id")),";") - 8),$24":"$27"-"$28,$7}' | sed s/[\"\;]//g | sort | uniq  > $out_dir/combined.intron.gtf.UTR5.wo
perl perl/filter.cds.intron.pl $out_dir/combined.intron.gtf.utr3.wo $gtf $out_dir/combined.intron.gtf.utr3.wo.addcds 

./count_intron_bundary_reads.sh $data_path $read_length $out_dir $strandness $readend 1>$out_dir/count.reads.log 2>$out_dir/count.reads.elog
if [ -f ./$out_dir/count ]; then
	mkdir $out_dir/count
fi
mv $out_dir/*.intron.left $out_dir/intron/
mv $out_dir/*.intron.right $out_dir/intron/
paste $out_dir/*.intron.left |awk '$1 ~ ":"' | awk '{for(i=7;i<=NF;i+=7) {printf("%s\t",$i)} printf("\n")}' > $out_dir/combined.intron.left.count
paste $out_dir/*.intron.right |awk '$1 ~ ":"' | awk '{for(i=7;i<=NF;i+=7) {printf("%s\t",$i)} printf("\n")}' > $out_dir/combined.intron.right.count

lf=`ls -l $out_dir/count/*.intron.left | head -1 |awk '{print $9}'`
rf=`ls -l $out_dir/count/*.intron.right | head -1 |awk '{print $9}'`
awk '$1 ~ ":"' $lf | cut -f 1 > $out_dir/combined.intron.left.spid
awk '$1 ~ ":"' $rf | cut -f 1 > $out_dir/combined.intron.right.spid
errorline=`diff $out_dir/combined.intron.left.spid $out_dir/combined.intron.right.spid | wc -l`
echo $errorline line are different between left and right
 
paste $out_dir/combined.intron.left.spid $out_dir/combined.intron.left.count > $out_dir/combined.intron.left.count.withspid
rm -f $out_dir/combined.intron.left.count
paste $out_dir/combined.intron.right.spid $out_dir/combined.intron.right.count > $out_dir/combined.intron.right.count.withspid
rm -f $out_dir/combined.intron.right.count

awk '{OFS="\t"}{print $1,$4,$5,$10,$NF,$7}' $out_dir/combined.intron.gtf.utr3 | sortBed -i - > $out_dir/combined.intron.utr3.bed
Rscript R/filter.splicing.events.r $out_dir

# rm -f $out_dir/count
# rm -f $out_dir/intron
# rm -f $out_dir/*.spid
# rm -f *.log 
# rm -f *.elog 


 



