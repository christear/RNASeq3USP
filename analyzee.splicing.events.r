library(ggpubr)
library(gplots)
library(circlize)

source("survive.analysis.functions.v2.r")

###
percut = 0.1
juncut = 1
cls = get_palette("lancet",7)
##
get_num = function(x,pcut = 0.1,subspid = subspid,utr3intr = utr3intr,annoutr3intr = annoutr3intr,unanutr3intr = unanutr3intr){
    subx = x[!(is.na(x))]
    eachid = subspid[!(is.na(x))]
    tn = sum(subx > pcut)
    itn = length(intersect(eachid[subx > pcut],utr3intr[,1]))
    gtn = length(unique(utr3intr[utr3intr[,1] %in% eachid[subx > pcut],7]))
    atn = length(intersect(eachid[subx > pcut],annoutr3intr[,1]))
    utn = length(intersect(eachid[subx > pcut],unanutr3intr[,1]))
    numvec = c(tn,itn,gtn,atn,utn)
    numvec
}

### calculate distance between splicing site (SS5 and SS3) and utr3:start-end ..
# r is the distance between utr3 end (PAS) and SS3
# l is the distance between utr3 start (stop codon) and SS5
caldis = function(x,method = "r"){
    spsp = unlist(strsplit(as.character(x[1]),"[:-]"))
    utrp = unlist(strsplit(as.character(x[9]),"[:-]"))
    if(method == "r"){
        dis = as.numeric(spsp[3]) - as.numeric(utrp[3])
        if(x[10] == "-"){
            dis = as.numeric(utrp[2]) - as.numeric(spsp[2])
        }
    }else{
        dis = as.numeric(spsp[2]) - as.numeric(utrp[2])
        if(x[10] == "-"){
            dis = as.numeric(utrp[3]) - as.numeric(spsp[3])
        }
    }
    dis
}

#cancertype="LIHC"
argv=commandArgs(TRUE)
cancertype = as.character(argv[1])

cat("loading the filtered splicing count data for",cancertype,"cohort...\n")
#load(paste("./Rdata/",cancertype,".filtered.splicing.count.Rdata",sep = ""))
#load(paste("./Rdata/",cancertype,".filtered.splicing.count.new.Rdata",sep = ""))
load(paste("./Rdata/",cancertype,".filtered.splicing.count.v2.Rdata",sep = ""))
##

##
tcgasummary = read.table("~/TCGA_RNA_summary.txt",head = T,sep = "\t")
eachsummary = tcgasummary[tcgasummary$analysis_id %in% colnames(subsout),]
eachsummary = eachsummary[order(eachsummary$analysis_id),]
eachntsamid = substr(eachsummary$barcode[eachsummary$sample_type == "NT"],0,12)
eachtpsamid = substr(eachsummary$barcode[eachsummary$sample_type == "TP"],0,12)
eachsamid = substr(eachsummary$barcode,0,12)
eachntsamid2 = names(table(eachntsamid))[table(eachntsamid) == 1]
eachtpsamid2 = names(table(eachtpsamid))[table(eachtpsamid) == 1]
pairedntsum = eachsummary[eachsummary$sample_type == "NT" & eachsamid %in% intersect(eachntsamid2,eachtpsamid2),]
pairedtpsum = eachsummary[eachsummary$sample_type == "TP" & eachsamid %in% intersect(eachntsamid2,eachtpsamid2),]
##
cat("summarizing splicing events for each sample...\n")
m1 = rep(0,nrow(utr3intr))
m1[grep("UTR",utr3intr[,11])] = 1
annoutr3intr = utr3intr[utr3intr[,3] == 1 | m1,]
unanutr3intr = utr3intr[utr3intr[,3] == 0,]

#subrsums = sublc + subrc + subsout * 2
#colnames(subrsums) = colnames(subratio)
#subrsums_tp = subrsums[,colnames(subrsums) %in% tcgasummary$analysis_id[tcgasummary$sample_type == "TP"]]
#subrsums_nt = subrsums[,colnames(subrsums) %in% tcgasummary$analysis_id[tcgasummary$sample_type == "NT"]]

### more than one supporting junctions
juncfilterd = apply(subsout,1,sum) > juncut
sublc2 = sublc[juncfilterd,]
subrc2 = subrc[juncfilterd,]
subsout2 = subsout[juncfilterd,]
subratio2 = subratio[juncfilterd,]

subrsums2 = sublc2 + subrc2 + subsout2 * 2
colnames(subrsums2) = colnames(subratio2)
subrsums2_tp = subrsums2[,colnames(subrsums2) %in% tcgasummary$analysis_id[tcgasummary$sample_type == "TP"]]
subrsums2_nt = subrsums2[,colnames(subrsums2) %in% tcgasummary$analysis_id[tcgasummary$sample_type == "NT"]]

### remove those lowly expressed splicing events ...
subratio2_m = as.matrix(subratio2)
subratio2_m[which(subrsums2 < 3)] = NA

###
spgenesnum = t(apply(subratio2,2,function(x) get_num(x,pcut = 0.1,rownames(subratio2,2),utr3intr,annoutr3intr,unanutr3intr)))
colnames(spgenesnum) = c("total_spnum","total_in_spnum","total_in_genenum","total_in_anno","total_in_unanno")
spgenesnum = spgenesnum[order(rownames(spgenesnum)),]
pariedntspgn = spgenesnum[rownames(spgenesnum) %in% pairedntsum$analysis_id,]
pariedtpspgn = spgenesnum[rownames(spgenesnum) %in% pairedtpsum$analysis_id,]

if(sum(rownames(pariedntspgn) == pairedntsum$analysis_id) == nrow(pairedntsum) & sum(rownames(pariedtpspgn) == pairedtpsum$analysis_id) == nrow(pariedtpspgn)){
    pariedntspgn = pariedntspgn[order(substr(pairedntsum$barcode,0,12)),]
    pariedtpspgn = pariedtpspgn[order(substr(pairedtpsum$barcode,0,12)),]
}else{
    cat("samples are not appropriately paired ...\n")
    pariedntspgn = pariedntspgn[order(rownames(pariedntspgn)),]
    pariedtpspgn = pariedtpspgn[order(rownames(pariedtpspgn)),]
    pariedntspgn = pariedntspgn[order(substr(pairedntsum$barcode,0,12)),]
    pariedtpspgn = pariedtpspgn[order(substr(pairedtpsum$barcode,0,12)),]
}

#subratio3 = subratio2[apply(subrsums2_nt,1,sum) > 1 & apply(subrsums2_tp,1,sum) > 1,]
subratio3 = subratio2_m[apply(subrsums2_nt,1,sum) > 1 & apply(subrsums2_tp,1,sum) > 1,]
subratio3_tp = subratio3[,colnames(subratio3) %in% eachsummary$analysis_id[eachsummary$sample_type == "TP"]]
subratio3_nt = subratio3[,colnames(subratio3) %in% eachsummary$analysis_id[eachsummary$sample_type == "NT"]]
#subratio3_ntn = apply(subratio3_nt,1,function(x) sum(x[!(is.na(x))] > 0.1))
#subratio3_tpn = apply(subratio3_tp,1,function(x) sum(x[!(is.na(x))] > 0.1))
#subratio3_tpn2 = apply(subratio3_tp[apply(subratio3_nt,1,function(x) sum(x[!(is.na(x))] > 0.1)) == 0,],1,function(x) sum(x[!(is.na(x))] > 0.1))
tpspid3 = rownames(subratio3_tp)[apply(subratio3_tp,1,function(x) sum(x[!(is.na(x))] > 0.1)) > 0]
ntspid3 = rownames(subratio3_nt)[apply(subratio3_nt,1,function(x) sum(x[!(is.na(x))] > 0.1)) > 0]
pairedtpspid3 = rownames(subratio3_tp)[apply(subratio3_tp[,colnames(subratio3_tp) %in% eachsummary$analysis_id[eachsamid %in% intersect(eachntsamid2,eachtpsamid2)]],1,function(x) sum(x[!(is.na(x))] > 0.1)) > 0]

spgenesnum3 = t(apply(subratio3,2,function(x) get_num(x,pcut = 0.1,rownames(subratio3,2),utr3intr,annoutr3intr,unanutr3intr)))
colnames(spgenesnum3) = c("total_spnum","total_in_spnum","total_in_genenum","total_in_anno","total_in_unanno")
spgenesnum3 = spgenesnum3[order(rownames(spgenesnum3)),]

## process 3'UTR splicing events for paired sample
cat("processing 3'UTR splicing events for paired sample ...\n")
pariedntspgn3 = spgenesnum3[rownames(spgenesnum3) %in% pairedntsum$analysis_id,]
pariedtpspgn3 = spgenesnum3[rownames(spgenesnum3) %in% pairedtpsum$analysis_id,]
if(sum(rownames(pariedntspgn3) == pairedntsum$analysis_id) == nrow(pairedntsum) & sum(rownames(pariedtpspgn3) == pairedtpsum$analysis_id) == nrow(pariedtpspgn3)){
    pariedntspgn3 = pariedntspgn3[order(substr(pairedntsum$barcode,0,12)),]
    pariedtpspgn3 = pariedtpspgn3[order(substr(pairedtpsum$barcode,0,12)),]
}else{
    cat("samples are not appropriately paired ...\n")
    pariedntspgn3 = pariedntspgn3[order(rownames(pariedntspgn3)),]
    pariedtpspgn3 = pariedtpspgn3[order(rownames(pariedtpspgn3)),]
    pariedntspgn3 = pariedntspgn3[order(substr(pairedntsum$barcode,0,12)),]
    pariedtpspgn3 = pariedtpspgn3[order(substr(pairedtpsum$barcode,0,12)),]
}

### process total splicing events for paired sample
cat("processing total splicing events for paired sample ...\n")
pariedntallspn = allspn[names(allspn) %in% pairedntsum$analysis_id]
pariedtpallspn = allspn[names(allspn) %in% pairedtpsum$analysis_id]
names(pariedntallspn) = names(allspn)[names(allspn) %in% pairedntsum$analysis_id]
names(pariedtpallspn) = names(allspn)[names(allspn) %in% pairedtpsum$analysis_id]
if(sum(names(pariedntallspn) == pairedntsum$analysis_id) == length(pariedntallspn) & sum(names(pariedtpallspn) == pairedtpsum$analysis_id) == length(pariedtpallspn)){
    pariedntallspn = pariedntallspn[order(substr(pairedntsum$barcode,0,12))]
    pariedtpallspn = pariedtpallspn[order(substr(pairedtpsum$barcode,0,12))]
}else{
    cat("samples are not appropriately paired ...\n")
    pariedntallspn = pariedntallspn[order(names(pariedntallspn)),]
    pariedtpallspn = pariedtpallspn[order(names(pariedtpallspn)),]
    pariedntallspn = pariedntallspn[order(substr(pairedntsum$barcode,0,12))]
    pariedtpallspn = pariedtpallspn[order(substr(pairedtpsum$barcode,0,12))]
}

spnda = data.frame(total_sp = allspn,utr3_sp = spgenesnum3[,2],type = factor(as.character(eachsummary$sample_type)))
spnda$utr3_ratio = spnda$utr3_sp/spnda$total_sp * 100
pairedspnda = data.frame(total_sp = c(pariedntallspn,pariedtpallspn),utr3_sp = c(pariedntspgn[,2],pariedtpspgn[,2]),type = c(rep("normal",nrow(pariedntspgn)),rep("tumor",nrow(pariedtpspgn))))
pairedspnda$utr3_ratio = pairedspnda$utr3_sp/pairedspnda$total_sp * 100

### association with OVS
cat("checking the association with OVS ...\n")
#spnummat = cbind(allspn,spgenesnum3[,c(1,3)],spgenesnum[,c(1,3)])
#spnummat = cbind(spnda$total_sp,spnda$utr3_sp,spnda$utr3_ratio)
spnummat = cbind(allspn,spgenesnum3[,2],spgenesnum3[,2]/allspn)
tpspnummat = spnummat[rownames(spnummat) %in% eachsummary$analysis_id[eachsummary$sample_type == "TP"],]
rownames(tpspnummat) = substr(eachsummary$barcode[eachsummary$sample_type == "TP"],0,12)

clinfile = paste("~/data/TCGA_clinical/",cancertype,".merged_only_clinical_clin_format.txt",sep = "")
clinical = as.data.frame(t(read.table(clinfile,head = T,sep = "\t",row.names = 1,quote = "",comment = "")))
clinicalout = processClinicalTab(clinical)
clinicalout = clinicalout[!(is.na(clinicalout[,5])),]
osam = intersect(rownames(tpspnummat),rownames(clinicalout))
subclinicalout = clinicalout[rownames(clinicalout) %in% osam,]
#subtpspnummat = tpspnummat[rownames(tpspnummat) %in% osam,]
calmean = function(vec){
    if(is.null(dim(vec))){
        mv = vec
    }else{
        mv = apply(vec,2,function(x) mean(x[!(is.na(x))]))
    }
    mv
}
subtpspnummat = t(sapply(1:length(osam),function(i) calmean(tpspnummat[rownames(tpspnummat) == osam[i],])))
colnames(subtpspnummat) = c("total_spnum","UTR3_spnum","UTR3_spratio")
rownames(subtpspnummat) = osam
subclinicalout = subclinicalout[order(rownames(subclinicalout)),]
subtpspnummat = subtpspnummat[order(rownames(subtpspnummat)),]
sfitl = apply(subtpspnummat,2,function(x) dosurv(subclinicalout$new_death,subclinicalout$death_event,x,quan = 1/3,ifordered = TRUE))
sfitl2 = apply(subtpspnummat,2,function(x) dosurv(subclinicalout$new_death,subclinicalout$death_event,x,quan = 1/4,ifordered = T))
sfitl3 = apply(subtpspnummat,2,function(x) dosurv(subclinicalout$new_death,subclinicalout$death_event,x,quan = 1/2,ifordered = T))

###

### check the splice site
caldinuf = function(mat,method = "freq"){
    fre = c(sum(mat[,2] == 1) + sum(mat[,2] == 2),sum(mat[,2] == 3) + sum(mat[,2] == 4),sum(mat[,2] == 5) + sum(mat[,2] == 6))
    if(method == "freq"){
        fre = round(fre/sum(fre) * 100,2)
    }
    names(fre) = c("GT/AG","GC/AG","AT/AC")
    fre
}

alldinucleotidef = c(405,963633,951712,36945,37003,7386,3121)
alldinucleotidef2 = c(alldinucleotidef[2] + alldinucleotidef[3],alldinucleotidef[4] + alldinucleotidef[5],alldinucleotidef[6] + alldinucleotidef[7])
names(alldinucleotidef2) = c("GT/AG","GC/AG","AT/AC")
alldinucleotidep = round(alldinucleotidef2/sum(alldinucleotidef2) * 100,2)
uniqutr3intr = unique(utr3intr[,1:2])

dnfda = data.frame(frequency = c(alldinucleotidep,caldinuf(uniqutr3intr),caldinuf(uniqutr3intr[uniqutr3intr[,1] %in% c(tpspid3,ntspid3),])),introntype = rep(c("all_introns","utr3_introns","utr3_introns(>10%)"),each = 3),sstype = rep(c("GT/AG","GC/AG","AT/AC"),3))
dnfda2 = data.frame(frequency = c(caldinuf(uniqutr3intr[uniqutr3intr[,1] %in% setdiff(ntspid3,tpspid3),]),caldinuf(uniqutr3intr[uniqutr3intr[,1] %in% setdiff(tpspid3,ntspid3),]),caldinuf(uniqutr3intr[uniqutr3intr[,1] %in% intersect(tpspid3,ntspid3),])),introntype = rep(c("normal","tumor","shared"),each = 3),sstype = rep(c("GT/AG","GC/AG","AT/AC"),3))


cp12 = chisq.test(cbind(alldinucleotidef2,caldinuf(uniqutr3intr,method = "n")))$p.value
cp13 = chisq.test(cbind(alldinucleotidef2,caldinuf(uniqutr3intr[uniqutr3intr[,1] %in% c(tpspid3,ntspid3),],method = "n")))$p.value
cp23 = chisq.test(cbind(caldinuf(uniqutr3intr,method = "n"),caldinuf(uniqutr3intr[uniqutr3intr[,1] %in% c(tpspid3,ntspid3),],method = "n")))$p.value

dnfmat = cbind(alldinucleotidef2,caldinuf(uniqutr3intr,"n"),caldinuf(uniqutr3intr[uniqutr3intr[,1] %in% c(tpspid3,ntspid3),],"n"))
dnfmat2 = cbind(caldinuf(uniqutr3intr[uniqutr3intr[,1] %in% setdiff(ntspid3,tpspid3),],"n"),caldinuf(uniqutr3intr[uniqutr3intr[,1] %in% setdiff(tpspid3,ntspid3),],"n"),caldinuf(uniqutr3intr[uniqutr3intr[,1] %in% intersect(ntspid3,tpspid3),],"n"))
colnames(dnfmat) = c("all_introns","utr3_introns","utr3_introns(>10%)")
colnames(dnfmat2) = c("normal","tumor","shared")

pvec1 = c(chisq.test(dnfmat[,c(1,2)])$p.value,chisq.test(dnfmat[,c(1,3)])$p.value,chisq.test(dnfmat[,c(2,3)])$p.value)
pvec2 = c(chisq.test(dnfmat2[,c(1,2)])$p.value,chisq.test(dnfmat2[,c(1,3)])$p.value,chisq.test(dnfmat2[,c(2,3)])$p.value)

pdf(paste("TCGA.",cancertype,".spliced.utr3.dinucleotide.v3.pdf",sep = ""))
ggbarplot(data = dnfda, x = "introntype", y = "frequency",fill = "sstype", palette = "npg", add = "jitter")  -> fp1
ggbarplot(data = dnfda2, x = "introntype", y = "frequency",fill = "sstype", palette = "npg", add = "jitter") -> fp2
ggarrange(fp1,fp2,ncol = 3,nrow = 2) -> fp
print(fp)
dev.off()

cat("plotting the results...\n")
pdf(paste("TCGA.",cancertype,".spliced.utr3.num.v3.pdf",sep = ""))
ggboxplot(data = spnda[spnda$type == "TP" | spnda$type == "NT",], x = "type", y = "utr3_sp",color = "type", palette = "lancet", add = "jitter") + stat_compare_means() -> fp1
ggpaired(data = pairedspnda,x = "type",y = "utr3_sp",color = "type", line.color = "gray", line.size = 0.4,palette = "lancet")+ stat_compare_means(paired = TRUE) -> fp2
ggboxplot(data = spnda[spnda$type == "TP" | spnda$type == "NT",], x = "type", y = "total_sp",color = "type", palette = "lancet", add = "jitter") + stat_compare_means() -> fp3
ggpaired(data = pairedspnda,x = "type",y = "total_sp",color = "type", line.color = "gray", line.size = 0.4,palette = "lancet")+ stat_compare_means(paired = TRUE) -> fp4
ggboxplot(data = spnda[spnda$type == "TP" | spnda$type == "NT",], x = "type", y = "utr3_ratio",color = "type", palette = "lancet", add = "jitter") + stat_compare_means() -> fp5
ggpaired(data = pairedspnda,x = "type",y = "utr3_ratio",color = "type", line.color = "gray", line.size = 0.4,palette = "lancet")+ stat_compare_means(paired = TRUE) -> fp6
ggarrange(fp1,fp2,fp3,fp4,fp5,fp6,ncol = 3,nrow = 2) -> fp
print(fp)
lwd = 2
#venn(spidl)
p1 = 1 - pchisq(sfitl[[1]][[2]]$chisq, length(sfitl[[1]][[2]]$n) - 1)
p11 = summary(coxph(Surv(subclinicalout$new_death,subclinicalout$death_event) ~ subtpspnummat[,1]))$logtest[3]
plot(sfitl[[1]][[1]],main = "#total splicing events",col = cls[2:1])
legend("topright",legend = paste(c("H=","L=","pval=","coxp="),c(sfitl[[1]][[2]]$n,p1,p11),sep = ""),col = cls[2:1],lty = c(1,1,NA,NA))
p2 = 1 - pchisq(sfitl[[2]][[2]]$chisq, length(sfitl[[2]][[2]]$n) - 1)
p12 = summary(coxph(Surv(subclinicalout$new_death,subclinicalout$death_event) ~ subtpspnummat[,2]))$logtest[3]
plot(sfitl[[2]][[1]],main = "#3'UTR splicing events",col = cls[2:1])
legend("topright",legend = paste(c("H=","L=","pval="),c(sfitl[[2]][[2]]$n,p2,p12),sep = ""),col = c(cls[2:1],NA,NA),lty = c(1,1,NA,NA))
p3 = 1 - pchisq(sfitl[[3]][[2]]$chisq, length(sfitl[[3]][[2]]$n) - 1)
p13 = summary(coxph(Surv(subclinicalout$new_death,subclinicalout$death_event) ~ subtpspnummat[,3]))$logtest[3]
plot(sfitl[[3]][[1]],main = "#3'UTR splicing ratios",col = cls[2:1])
legend("topright",legend = paste(c("H=","L=","pval="),c(sfitl[[3]][[2]]$n,p3,p13),sep = ""),col = cls[2:1],lty = c(1,1,NA,NA))
##
p21 = 1 - pchisq(sfitl2[[1]][[2]]$chisq, length(sfitl2[[1]][[2]]$n) - 1)
plot(sfitl2[[1]][[1]],main = "#total splicing events",col = cls[2:1])
legend("topright",legend = paste(c("H=","L=","pval=","coxp="),c(sfitl2[[1]][[2]]$n,p21,p11),sep = ""),col = cls[2:1],lty = c(1,1,NA,NA))
p22 = 1 - pchisq(sfitl2[[2]][[2]]$chisq, length(sfitl2[[2]][[2]]$n) - 1)
plot(sfitl2[[2]][[1]],main = "#3'UTR splicing events",col = cls[2:1])
legend("topright",legend = paste(c("H=","L=","pval="),c(sfitl2[[2]][[2]]$n,p22,p12),sep = ""),col = cls[2:1],lty = c(1,1,NA,NA))
p23 = 1 - pchisq(sfitl2[[3]][[2]]$chisq, length(sfitl2[[3]][[2]]$n) - 1)
plot(sfitl2[[3]][[1]],main = "#3'UTR splicing ratios",col = cls[2:1])
legend("topright",legend = paste(c("H=","L=","pval="),c(sfitl2[[3]][[2]]$n,p23,p13),sep = ""),col = cls[2:1],lty = c(1,1,NA,NA))
##
p31 = 1 - pchisq(sfitl3[[1]][[2]]$chisq, length(sfitl3[[1]][[2]]$n) - 1)
plot(sfitl3[[1]][[1]],main = "#total splicing events",col = cls[2:1])
legend("topright",legend = paste(c("H=","L=","pval=","coxp="),c(sfitl3[[1]][[2]]$n,p31,p11),sep = ""),col = cls[2:1],lty = c(1,1,NA,NA))
p32 = 1 - pchisq(sfitl3[[2]][[2]]$chisq, length(sfitl3[[2]][[2]]$n) - 1)
plot(sfitl3[[2]][[1]],main = "#3'UTR splicing events",col = cls[2:1])
legend("topright",legend = paste(c("H=","L=","pval="),c(sfitl3[[2]][[2]]$n,p32,p12),sep = ""),col = cls[2:1],lty = c(1,1,NA,NA))
p33 = 1 - pchisq(sfitl3[[3]][[2]]$chisq, length(sfitl3[[3]][[2]]$n) - 1)
plot(sfitl3[[3]][[1]],main = "#3'UTR splicing ratios",col = cls[2:1])
legend("topright",legend = paste(c("H=","L=","pval="),c(sfitl3[[3]][[2]]$n,p33,p13),sep = ""),col = cls[2:1],lty = c(1,1,NA,NA))
dev.off()

cat("saving the data ...\n")
save(list = ls(),file = paste("./Rdata/",cancertype,".splicing.events.v2.Rdata",sep = ""))
#save(list = ls(),file = paste("./Rdata/",cancertype,".splicing.events.v3.Rdata",sep = ""))

#events = "3:41280846-41281309"
#tmp =


### old
if(FALSE){
    ggboxplot(data = tnda[tnda$type == "TP" | tnda$type == "NT",], x = "type", y = "total_in_spnum",color = "type", palette = "lancet", add = "jitter") + stat_compare_means() -> fp1
    ggboxplot(data = tnda[tnda$type == "TP" | tnda$type == "NT",], x = "type", y = "total_in_genenum",color = "type", palette = "lancet", add = "jitter") + stat_compare_means() -> fp2
    ggboxplot(data = tnda[tnda$type == "TP" | tnda$type == "NT",], x = "type", y = "total_in_anno",color = "type", palette = "lancet", add = "jitter") + stat_compare_means() -> fp3
    ggboxplot(data = tnda[tnda$type == "TP" | tnda$type == "NT",], x = "type", y = "total_in_unanno",color = "type", palette = "lancet", add = "jitter") + stat_compare_means() -> fp4
    ggpaired(data = pairedda,x = "type",y = "number",color = "type", line.color = "gray", line.size = 0.4,palette = "lancet")+ stat_compare_means(paired = TRUE) -> fp5
    ggboxplot(data = alltnda[alltnda$type == "TP" | alltnda$type == "NT",], x = "type", y = "number",color = "type", palette = "lancet", add = "jitter", ylim = c(0,220000))  + stat_compare_means() -> fp71
    ggarrange(fp1,fp2,fp3,fp4,fp5,fp71,ncol = 3,nrow = 2) -> fp
    print(fp)
    ggboxplot(data = tnda3[tnda3$type == "TP" | tnda3$type == "NT",], x = "type", y = "total_in_spnum",color = "type", palette = "lancet", add = "jitter") + stat_compare_means() -> fp1
    ggboxplot(data = tnda3[tnda3$type == "TP" | tnda3$type == "NT",], x = "type", y = "total_in_genenum",color = "type", palette = "lancet", add = "jitter") + stat_compare_means() -> fp2
    ggboxplot(data = tnda3[tnda3$type == "TP" | tnda3$type == "NT",], x = "type", y = "total_in_anno",color = "type", palette = "lancet", add = "jitter") + stat_compare_means() -> fp3
    ggboxplot(data = tnda3[tnda3$type == "TP" | tnda3$type == "NT",], x = "type", y = "total_in_unanno",color = "type", palette = "lancet", add = "jitter") + stat_compare_means() -> fp4
    ggpaired(data = pairedda3,x = "type",y = "number",color = "type", line.color = "gray", line.size = 0.4,palette = "lancet")+ stat_compare_means(paired = TRUE) -> fp5
    ggpaired(data = pairedalltnda,x = "type",y = "number",color = "type", line.color = "gray", line.size = 0.4,palette = "lancet",ylim = c(0,220000))+ stat_compare_means(paired = TRUE) -> fp72
    ggarrange(fp1,fp2,fp3,fp4,fp5,fp72,ncol = 3,nrow = 2) -> fp
    print(fp)
}


