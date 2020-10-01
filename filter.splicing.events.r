#library(parallel)
###
# filtered splicing events, those
#1,junction linked annotated coding exons
#2,junction within annotated coding exons
#3,junction overlapped with 5' UTR 
#
###

argv=commandArgs(TRUE)
cancertype = as.character(argv[1])

cat("...start processing ",cancertype,"...\n")
utr3intr = read.table(paste(cancertype,".combined.intron.gtf.utr3.wo.addcds.adddis",sep = ""))
utr3intrl = read.table(paste(cancertype,".combined.intron.left.gtf.utr3.wo.addcds.adddis",sep = ""))
utr3intrr = read.table(paste(cancertype,".combined.intron.right.gtf.utr3.wo.addcds.adddis",sep = ""))

cdsintr = read.table(paste(cancertype,".combined.intron.gtf.cds.wo",sep = ""))
utr5intr = read.table(paste(cancertype,".combined.intron.gtf.UTR5.wo",sep = ""))
cat(length(unique(utr3intr[,1]))," introns in total...\n")
### remove junction linked coding exons
cat("removing junctions linked conding exons...\n")
utr3intr = utr3intr[utr3intr[,1] %in% setdiff(utr3intr[,1],utr3intr[grep("CDS",utr3intr[,11]),1]),]
cat(length(unique(utr3intr[,1]))," introns left...\n")
#@utr3intr2 = utr3intr2[utr3intr2[,1] %in% setdiff(utr3intr2[,1],utr3intr2[grep("CDS",utr3intr2[,11]),1]),]
### remove junction within coding exon
cat("removing junctions within conding exons...\n")
utr3intr_p1 = utr3intr[utr3intr[,1] %in% setdiff(utr3intr[,1],cdsintr[,1]),]
utr3intr_p2 = utr3intr[utr3intr[,1] %in% cdsintr[,1],]
utr3intr = rbind(utr3intr_p1,utr3intr_p2[grep("UTR:UTR",utr3intr_p2[,11]),])
cat(length(unique(utr3intr[,1]))," introns left...\n")
#@utr3intr2_p1 = utr3intr2[utr3intr2[,1] %in% setdiff(utr3intr2[,1],cdsintr[,1]),]
#@utr3intr2_p2 = utr3intr2[utr3intr2[,1] %in% cdsintr[,1],]
#@utr3intr2 = rbind(utr3intr2_p1,utr3intr2_p2[grep("UTR:UTR",utr3intr2_p2[,11]),])
### remove junction overlap with 5' UTR
cat("removing junctions within 5' UTR...\n")
utr3intr = utr3intr[utr3intr[,1] %in% setdiff(utr3intr[,1],utr5intr[,1]),]
cat(length(unique(utr3intr[,1]))," introns left...\n")
#@utr3intr2 = utr3intr2[utr3intr2[,1] %in% setdiff(utr3intr2[,1],utr5intr[,1]),]
bed = read.table(paste(cancertype,".combined.intron.addchr.bed",sep = ""))
subbed = bed[bed[,4] %in% utr3intr[,1],]
write.table(subbed,file = paste(cancertype,".filtered.intron.addchr.bed",sep = ""),col.names = F,row.names = F,sep = "\t",quote = F)
rm(bed)
### read splicing junction counting data
cat("reading splicng junctions counting data...\n")
utr3ratio = read.table(paste("./ratio/",cancertype,".utr3.ratio",sep = ""))
utr3sout = read.table(paste("./count/",cancertype,".combined.splicing.out.utr3.txt",sep = ""))
utr3lc = read.table(paste("./count_fragment/",cancertype,".left.count.withspid.utr3",sep = ""))
utr3rc = read.table(paste("./count_fragment/",cancertype,".right.count.withspid.utr3",sep = ""))
sublc = utr3lc[utr3lc[,1] %in% utr3intr[,1],2:ncol(utr3lc)]
subrc = utr3rc[utr3rc[,1] %in% utr3intr[,1],2:ncol(utr3rc)]
subsout = utr3sout[utr3sout[,1] %in% utr3intr[,1],2:ncol(utr3sout)]

rownames(sublc) = utr3lc[utr3lc[,1] %in% utr3intr[,1],1]
rownames(subrc) = utr3rc[utr3rc[,1] %in% utr3intr[,1],1]
rownames(subsout) = utr3sout[utr3sout[,1] %in% utr3intr[,1],1]

colnames(sublc) = gsub("\\..*","",dir(paste("./count_fragment/",cancertype,sep = ""))[grep("*.intron.left$",dir(paste("./count_fragment/",cancertype,sep = "")))])
colnames(subrc) = gsub("\\..*","",dir(paste("./count_fragment/",cancertype,sep = ""))[grep("*.intron.right$",dir(paste("./count_fragment/",cancertype,sep = "")))])
colnames(subsout) = gsub("\\..*","",dir(paste("./intron/",cancertype,sep = ""))[grep("*.intron.gtf",dir(paste("./intron/",cancertype,sep = "")))])

#subsums = sublc + subrc + subsout * 2
subratio = 2*subsout/(sublc + subrc + subsout * 2)

fillNA = function(x){
    subx = x
    subx[x == -1] = NA
    subx
}

subr2 = utr3ratio[utr3ratio[,1] %in% utr3intr[,1],2:ncol(utr3ratio)]
#cat("filling NAs....\n")
subr2 = apply(subr2,2,function(x) fillNA(x))
#cat("finished filling NA...\n")
rownames(subr2) = utr3ratio[utr3ratio[,1] %in% utr3intr[,1],1]
colnames(subr2) = colnames(subratio)
#subratio2 = subratio[rownames(subratio) %in% rownames(subr2),]

allspn = read.table(paste("./ratio/",cancertype,".total.spnum",sep = ""))
allspn = unlist(allspn[,2])
names(allspn) = colnames(subr2)
cat("saving data ...\n")
save(list = c("subsout","utr3intr","allspn","subratio","sublc","subrc"),file = paste("./Rdata/",cancertype,".filtered.splicing.count.v2.Rdata",sep = ""))





