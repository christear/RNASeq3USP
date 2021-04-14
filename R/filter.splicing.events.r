###
# filtered splicing events, those
#1,junction linked annotated coding exons
#2,junction within annotated coding exons
#3,junction overlapped with 5' UTR 
#
###
percut = 0.1
juncut = 1
###
argv = commandArgs(TRUE)
out_dir = as.character(argv[1])

cat("filtering 3'UTR splicing events...\n")
utr3intr = read.table(paste(out_dir,"/combined.intron.gtf.utr3.wo.addcds",sep = ""))
cdsintr = read.table(paste(out_dir,"/combined.intron.gtf.cds.wo",sep = ""))
utr5intr = read.table(paste(out_dir,"/combined.intron.gtf.UTR5.wo",sep = ""))
cat(length(unique(utr3intr[,1]))," introns in total...\n")
### remove junction linked coding exons
cat("removing junctions linked conding exons...\n")
utr3intr = utr3intr[utr3intr[,1] %in% setdiff(utr3intr[,1],utr3intr[grep("CDS",utr3intr[,11]),1]),]
cat(length(unique(utr3intr[,1]))," introns left...\n")
### remove junction within coding exon
cat("removing junctions within conding exons...\n")
utr3intr_p1 = utr3intr[utr3intr[,1] %in% setdiff(utr3intr[,1],cdsintr[,1]),]
utr3intr_p2 = utr3intr[utr3intr[,1] %in% cdsintr[,1],]
utr3intr = rbind(utr3intr_p1,utr3intr_p2[grep("UTR:UTR",utr3intr_p2[,11]),])
cat(length(unique(utr3intr[,1]))," introns left...\n")
### remove junction overlap with 5' UTR
cat("removing junctions within 5' UTR...\n")
utr3intr = utr3intr[utr3intr[,1] %in% setdiff(utr3intr[,1],utr5intr[,1]),]
cat(length(unique(utr3intr[,1]))," introns left...\n")
#@utr3intr2 = utr3intr2[utr3intr2[,1] %in% setdiff(utr3intr2[,1],utr5intr[,1]),]
bed = read.table(paste(out_dir,"/combined.intron.utr3.bed",sep = ""))
subbed = bed[bed[,4] %in% utr3intr[,1],]
write.table(subbed,file = paste(out_dir,"/combined.intron.filtered.utr3.bed",sep = ""),col.names = F,row.names = F,sep = "\t",quote = F)
rm(bed)
### read splicing junction counting data
cat("reading splicng junctions counting data...\n")
utr3sout = read.table(paste(out_dir,"/combined.splicing.out.txt",sep = ""))
utr3lc = read.table(paste(out_dir,"/combined.intron.left.count.withspid",sep = ""))
utr3rc = read.table(paste(out_dir,"/combined.intron.right.count.withspid",sep = ""))
sublc = utr3lc[utr3lc[,1] %in% utr3intr[,1],2:ncol(utr3lc)]
subrc = utr3rc[utr3rc[,1] %in% utr3intr[,1],2:ncol(utr3rc)]
subsout = utr3sout[utr3sout[,1] %in% utr3intr[,1],2:ncol(utr3sout)]

rownames(sublc) = utr3lc[utr3lc[,1] %in% utr3intr[,1],1]
rownames(subrc) = utr3rc[utr3rc[,1] %in% utr3intr[,1],1]
rownames(subsout) = utr3sout[utr3sout[,1] %in% utr3intr[,1],1]

sampleid = read.table(paste(out_dir,"/sampleID.list",sep = ""))
colnames(sublc) = colnames(subrc) = colnames(subsout) = sampleid[,1]
subsums = sublc + subrc + subsout * 2
subratio = 2*subsout/(sublc + subrc + subsout * 2)

fillNA = function(x){
    subx = x
    subx[x == -1] = NA
    subx
}

### more than one supporting junctions
juncfilterd = apply(subsout2,1,sum) > juncut
sublc2 = sublc[juncfilterd,]
subrc2 = subrc[juncfilterd,]
subsout2 = subsout[juncfilterd,]
subratio2 = subratio[juncfilterd,]
cat("saving data ...\n")
save(list = c("out_dir","subsout2","utr3intr","subratio2","sublc2","subrc2"),file = paste(out_dir,"/filtered.splicing.count.Rdata",sep = ""))






