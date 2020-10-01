library(ggpubr,quietly = T)
library(gplots,quietly = T)

###
cls = get_palette("lancet",7)
#pcut = 0.05
pcut = 0.1
rpkmcut = 1
###
dowilcoxtest = function(x,y,paired = FALSE){
    if(sum(!(is.na(x))) < 10 | sum(!(is.na(y))) < 10){
        pval = NA
    }else{
        if(paired){
            if(sum(!(is.na(x + y))) > 5){
                p1 = wilcox.test(x,y,paired = paired,alternative = 'less')$p.value
                p2 = wilcox.test(x,y,paired = paired,alternative = 'greater')$p.value
                pval = min(p1,p2)
            }else{
                pval = NA
            }
        }else{
            p1 = wilcox.test(x,y,paired = paired,alternative = 'less')$p.value
            p2 = wilcox.test(x,y,paired = paired,alternative = 'greater')$p.value
            pval = min(p1,p2)
        }
    }
    pval
}

###
dottest = function(x,y,paired = FALSE){
    if(sum(!(is.na(x))) < 10 | sum(!(is.na(y))) < 10){
        pval = NA
    }else{
        if(sd(x[!(is.na(x))]) == 0 & sd(y[!(is.na(y))]) == 0){
            pval = NA
        }else{
            if(paired){
                if(sum(!(is.na(x + y))) > 5){
                    p1 = t.test(x,y,paired = paired,alternative = 'less')$p.value
                    p2 = t.test(x,y,paired = paired,alternative = 'greater')$p.value
                    pval = min(p1,p2)
                }else{
                    pval = NA
                }
            }else{
                p1 = t.test(x,y,paired = paired,alternative = 'less')$p.value
                p2 = t.test(x,y,paired = paired,alternative = 'greater')$p.value
                pval = min(p1,p2)
            }
        }
    }
    pval
}


###
argv=commandArgs(TRUE)
cancertype = as.character(argv[1])

cat("loading the data from",cancertype,"...\n")
#load(paste("./Rdata/",cancertype,".splicing.events.features.v2.Rdata",sep = ""))
load(paste("./Rdata/",cancertype,".splicing.events.v2.Rdata",sep = ""))
#load(paste("./Rdata/",cancertype,".splicing.events.v3.Rdata",sep = ""))

subratio3_tpn = apply(subratio3_tp,1,function(x) sum(x[!(is.na(x))] > 0.1))
subratio3_ntn = apply(subratio3_nt,1,function(x) sum(x[!(is.na(x))] > 0.1))
oid3 = c(tpspid3,ntspid3)
oid3 = names(subratio3_tpn)[subratio3_tpn > ncol(subratio3_tp) * 0.5]
cat("there are",ncol(subratio3_nt),"NT samples and ",ncol(subratio3_tp)," TP samples...\n")

## for each event which was detected in both tumor and normal and in more than half of tumor samples
## splicing levels in tumor were compared to that in adjancet normal
cat("comparing tumor to normal....\n")
subratio3_nttp = cbind(subratio3_nt[rownames(subratio3_nt) %in% oid3,],subratio3_tp[rownames(subratio3_tp) %in% oid3,])
subratio3_nttp = subratio3_nttp[order(rownames(subratio3_nttp)),]
subratio3_nt = subratio3_nt[order(rownames(subratio3_nt)),]
subratio3_tp = subratio3_tp[order(rownames(subratio3_tp)),]
write.table(subratio3_nt[rownames(subratio3_nt) %in% oid3,],file = paste("./txt/",cancertype,".normal.c3usp.ratio.txt",sep = ""),sep = "\t",quote = F)
write.table(subratio3_tp[rownames(subratio3_nt) %in% oid3,],file = paste("./txt/",cancertype,".tumor.c3usp.ratio.txt",sep = ""),sep = "\t",quote = F)

psppval1 = apply(subratio3_nttp,1,function(x) dowilcoxtest(x[1:ncol(subratio3_nt)],x[(ncol(subratio3_nt) + 1):length(x)]))
psppval2 = apply(subratio3_nttp,1,function(x) dottest(x[1:ncol(subratio3_nt)],x[(ncol(subratio3_nt) + 1):length(x)]))
psppadj1 = p.adjust(psppval1,'fdr')
psppadj2 = p.adjust(psppval2,'fdr')

subratio3_ntm1 = apply(subratio3_nt,1,function(x) median(x[!(is.na(x))]))
subratio3_tpm1 = apply(subratio3_tp,1,function(x) median(x[!(is.na(x))]))
subratio3_ntm2 = apply(subratio3_nt,1,function(x) mean(x[!(is.na(x))]))
subratio3_tpm2 = apply(subratio3_tp,1,function(x) mean(x[!(is.na(x))]))

fres1 = cbind(subratio3_ntm1[names(subratio3_ntm1) %in% oid3],subratio3_tpm1[names(subratio3_tpm1) %in% oid3],psppval1,psppadj1,subratio3_ntm2[names(subratio3_ntm2) %in% oid3],subratio3_tpm2[names(subratio3_tpm2) %in% oid3],psppval2,psppadj2)
#subfres1 = fres1[!(is.na(fres1[,3])) | !(is.na(fres1[,7])),]

#subfres2 = subfres1[!(is.na(subfres1[,5])),]
#sigfres1 = subfres1[subfres1[,4] < pcut,]
#sigfres2 = subfres2[subfres2[,6] < pcut,]

cat("comparing paired tumor to normal....\n")
subratio3_nt_paired = subratio3_nttp[,colnames(subratio3_nttp) %in% pairedntsum$analysis_id]
subratio3_tp_paired = subratio3_nttp[,colnames(subratio3_nttp) %in% pairedtpsum$analysis_id]
subratio3_nt_paired = subratio3_nt_paired[,order(pairedntsum$barcode)]
subratio3_tp_paired = subratio3_tp_paired[,order(pairedtpsum$barcode)]

psppval21 = apply(cbind(subratio3_nt_paired,subratio3_tp_paired),1,function(x) dowilcoxtest(x[1:ncol(subratio3_nt_paired)],x[(ncol(subratio3_nt_paired) + 1):(ncol(subratio3_nt_paired) + ncol(subratio3_tp_paired))],paired = T))
psppval22 = apply(cbind(subratio3_nt_paired,subratio3_tp_paired),1,function(x) dottest(x[1:ncol(subratio3_nt_paired)],x[(ncol(subratio3_nt_paired) + 1):(ncol(subratio3_nt_paired) + ncol(subratio3_tp_paired))],paired = T))
psppadj21 = p.adjust(psppval21,'fdr')
psppadj22 = p.adjust(psppval22,'fdr')

fres2 = cbind(apply(subratio3_nt_paired,1,function(x) median(x[!(is.na(x))])),apply(subratio3_tp_paired,1,function(x) median(x[!(is.na(x))])),apply(subratio3_tp_paired - subratio3_nt_paired,1,function(x) median(x[!(is.na(x))])),psppval21,psppadj21,apply(subratio3_nt_paired,1,function(x) mean(x[!(is.na(x))])),apply(subratio3_tp_paired,1,function(x) mean(x[!(is.na(x))])),apply(subratio3_tp_paired - subratio3_nt_paired,1,function(x) mean(x[!(is.na(x))])),psppval22,psppadj22)

## to keep the same dimenion of subres1 and subres2
#subfres2 = fres2[!(is.na(fres1[,4])) | !(is.na(fres1[,9])),]
##subfres21 = fres2[!(is.na(fres2[,5])),]
##subfres22 = fres2[!(is.na(fres2[,7])),]
##sigfres21 = subfres21[subfres21[,5] < pcut,]
##sigfres22 = subfres22[subfres22[,7] < pcut,]

if(cancertype == "BRCA"){
    cat("loading the data from GTEX normal breast tissue ...\n")
    load("../09_GTEX_data/Rdata/Breast.filtered.splicing.count.Rdata")
}else if(cancertype == "COAD"){
    cat("loading the data from GTEX normal colon tissue ...\n")
    load("../09_GTEX_data/Rdata/Colon.filtered.splicing.count.Rdata")
}else if(cancertype == "HNSC"){
    cat("no head and neck normal data found in Gtex ...\n")
    gtexr2q9 = t(apply(subratio3_nt,1,function(x) c(quantile(x[!(is.na(x))],c(0.05,0.1,0.5,0.9,0.95)),mean(x[!(is.na(x))]))))
    #    load("../09_GTEX_data/Rdata/Colon.filtered.splicing.count.Rdata")
}else if(cancertype == "KIRC"){
    cat("loading the data from GTEX normal kidney tissue ...\n")
    load("../09_GTEX_data/Rdata/Kidney.filtered.splicing.count.Rdata")
}else if(cancertype == "LIHC"){
    cat("loading the data from GTEX normal liver tissue ...\n")
    load("../09_GTEX_data/Rdata/Liver.filtered.splicing.count.Rdata")
}else if(cancertype == "LUAD"){
    cat("loading the data from GTEX normal lung tissue ...\n")
    load("../09_GTEX_data/Rdata/Lung.filtered.splicing.count.Rdata")
}else if(cancertype == "LUSC"){
    cat("loading the data from GTEX normal lung tissue ...\n")
    load("../09_GTEX_data/Rdata/Lung.filtered.splicing.count.Rdata")
}else if(cancertype == "PRAD"){
    cat("loading the data from GTEX normal prostate tissue ...\n")
    load("../09_GTEX_data/Rdata/Prostate.filtered.splicing.count.Rdata")
}else if(cancertype == "STAD"){
    cat("loading the data from GTEX normal stomach tissue ...\n")
    load("../09_GTEX_data/Rdata/Stomach.filtered.splicing.count.Rdata")
}else if(cancertype == "THCA"){
    cat("loading the data from GTEX normal thyroid tissue ...\n")
    load("../09_GTEX_data/Rdata/Thyroid.filtered.splicing.count.Rdata")
}

if(cancertype != "HNSC"){
    cat("there are ",nrow(gtexmat)," samples from GTEX ...\n")
}

## calculate overspliced sample numbers ...
tcgantrq9 = t(apply(subratio3_nt,1,function(x) c(quantile(x[!(is.na(x))],c(0.05,0.1,0.5,0.9,0.95)),mean(x[!(is.na(x))]))))
colnames(tcgantrq9) = c("q05","q10","q50","q90","q95","ave")

normalspq9 = sapply(1:nrow(fres1),function(i) max(tcgantrq9[rownames(tcgantrq9) == rownames(fres1)[i],4],gtexr2q9[rownames(gtexr2q9) == rownames(fres1)[i],4]))
normalspq1 = sapply(1:nrow(fres1),function(i) min(tcgantrq9[rownames(tcgantrq9) == rownames(fres1)[i],2],gtexr2q9[rownames(gtexr2q9) == rownames(fres1)[i],2]))
names(normalspq9) = names(normalspq1) = rownames(fres1)

#normalsp2q9 = sapply(1:nrow(fres2),function(i) max(tcgantrq9[rownames(tcgantrq9) == rownames(fres2)[i],4],gtexr2q9[rownames(gtexr2q9) == rownames(fres2)[i],4]))
#normalsp2q1 = sapply(1:nrow(fres2),function(i) min(tcgantrq9[rownames(tcgantrq9) == rownames(fres2)[i],2],gtexr2q9[rownames(gtexr2q9) == rownames(fres2)[i],2]))
#names(normalsp2q9) = names(normalsp2q1) = rownames(fres2)

tmptpr = subratio3_tp[rownames(subratio3_tp) %in% rownames(fres1),]
tmptpr = tmptpr[order(rownames(tmptpr)),]
tmpntr = subratio3_nt[rownames(subratio3_nt) %in% rownames(fres1),]
tmpntr = tmpntr[order(rownames(tmpntr)),]

eachtpsummary = eachsummary[eachsummary$sample_type == "TP",]
eachntsummary = eachsummary[eachsummary$sample_type == "NT",]

if(sum(colnames(tmptpr) == eachtpsummary$analysis_id) == ncol(tmptpr)){
    colnames(tmptpr) = substr(eachtpsummary$barcode,0,12)
}

### foreach event, calculate number of samples which are over/under-splicng
calnum = function(x,cut,compare = "upper",diffcut = 0){
    if(compare == "upper"){
        num = sum(x[!(is.na(x))] - cut > diffcut)
    }else if(compare == "lower"){
        num = sum(x[!(is.na(x))] - cut < -diffcut)
    }
    num
}
### foreach event, get the sample list which are over/under-splicing
getsam = function(x,cut,compare = "upper",samid,diffcut = 0){
    subsamid = samid[!(is.na(x))]
    subx = x[!(is.na(x))]
    if(compare == "upper"){
        sam = subsamid[subx - cut > diffcut]
    }else if(compare == "lower"){
        sam = subsamid[subx - cut < -diffcut]
    }
    sam
}

if(sum(rownames(tmptpr) == rownames(fres1)) == nrow(fres1)){
    overspnum = sapply(1:nrow(tmptpr),function(i) calnum(unlist(tmptpr[i,]),normalspq9[i]))
    downspnum = sapply(1:nrow(tmptpr),function(i) calnum(unlist(tmptpr[i,]),normalspq1[i],compare = "lower"))
    overspsaml = sapply(1:nrow(tmptpr),function(i) getsam(unlist(tmptpr[i,]),normalspq9[i],samid = colnames(tmptpr)))
    downspsaml = sapply(1:nrow(tmptpr),function(i) getsam(unlist(tmptpr[i,]),normalspq1[i],compare = "lower",samid = colnames(tmptpr)))
}

## add splicing levels from GTEX
if(length(setdiff(rownames(tmptpr),rownames(gtexr2q9))) > 0){
    addgtex = as.data.frame(matrix(NA,ncol = ncol(gtexr2q9),nrow = length(setdiff(rownames(tmptpr),rownames(gtexr2q9)))))
    rownames(addgtex) = setdiff(rownames(tmptpr),rownames(gtexr2q9))
    colnames(addgtex) = colnames(gtexr2q9)
    gtexr2q9 = rbind(gtexr2q9,addgtex)
}
sub1gtexr2q9 = gtexr2q9[rownames(gtexr2q9) %in% rownames(tmptpr),]
sub1gtexr2q9 = sub1gtexr2q9[order(rownames(sub1gtexr2q9)),]
##
gns = sapply(1:nrow(fres1),function(i) paste0(unique(as.character(utr3intr[utr3intr[,1] == rownames(fres1)[i],7])),collapse = ":"))

subfres1addn = cbind(fres1,sub1gtexr2q9[,c(3,6)],overspnum,downspnum,apply(tmptpr,1,function(x) sum(x[!(is.na(x))] > percut)),apply(tmpntr,1,function(x) sum(x[!(is.na(x))] > percut)),gns)
subfres2addn = cbind(fres2,sub1gtexr2q9[,c(3,6)],overspnum,downspnum,apply(tmptpr,1,function(x) sum(x[!(is.na(x))] > percut)),apply(tmpntr,1,function(x) sum(x[!(is.na(x))] > percut)),gns)

colnames(subfres1addn) = c(paste(c("NT","TP","pval","padj"),"median",sep = "_"),paste(c("NT","TP","pval","padj"),"mean",sep = "_"),"GTEX_median","GTEX_mean","over_spnum","down_spnum","tp_spnum","nt_spnum","gene")
colnames(subfres2addn) = c(paste(c("NT","TP","diff","pval","padj"),"median",sep = "_"),paste(c("NT","TP","diff","pval","padj"),"mean",sep = "_"),"GTEX_median","GTEX_mean","over_spnum","down_spnum","tp_spnum","nt_spnum","gene")
cat("outputing differentially splicing analysis results ...\n")
write.table(subfres1addn,file = paste("./txt/",cancertype,".sub.splicing.events.v1.txt",sep = ""),sep = "\t",quote = F)
write.table(subfres2addn,file = paste("./txt/",cancertype,".sub.splicing.events.v2.txt",sep = ""),sep = "\t",quote = F)

#oversppatl = lapply(overspsaml,function(x) unique(substr(eachsummary$barcode[eachsummary$analysis_id %in% x],0,12)))
#downsppatl = lapply(downspsaml,function(x) unique(substr(eachsummary$barcode[eachsummary$analysis_id %in% x],0,12)))
#oversppatl2 = lapply(overspsaml2,function(x) unique(substr(eachsummary$barcode[eachsummary$analysis_id %in% x],0,12)))
#downsppatl2 = lapply(downspsaml2,function(x) unique(substr(eachsummary$barcode[eachsummary$analysis_id %in% x],0,12)))

### checking spicing level across differernt stages ...
cat("checking splicing level across stages ...\n")
ntstage = rep(0,ncol(subratio3_nt))
tpstage = rep(NA,ncol(subratio3_tp))
tpstage[colnames(subratio3_tp) %in% eachtpsummary$analysis_id[substr(eachtpsummary$barcode,0,12) %in% rownames(clinicalout)[grep("stage i$|stage i[abcd]",clinicalout$stage)]]] = 1
tpstage[colnames(subratio3_tp) %in% eachtpsummary$analysis_id[substr(eachtpsummary$barcode,0,12) %in% rownames(clinicalout)[grep("stage ii$|stage ii[abcd]",clinicalout$stage)]]] = 2
tpstage[colnames(subratio3_tp) %in% eachtpsummary$analysis_id[substr(eachtpsummary$barcode,0,12) %in% rownames(clinicalout)[grep("stage iii$|stage iii[abcd]",clinicalout$stage)]]] = 3
tpstage[colnames(subratio3_tp) %in% eachtpsummary$analysis_id[substr(eachtpsummary$barcode,0,12) %in% rownames(clinicalout)[grep("stage iv$|stage iv[abcd]",clinicalout$stage)]]] = 4
#lmpvall = apply(subratio3_tp[rownames(subratio3_tp) %in% oid3,],1,function(x) summary(lm(x ~ tpstage))$coefficients[8])

numstages = c(ntstage,tpstage)
cal_sm = function(x){
    subx = x[!(is.na(numstages))]
    subnumstages = numstages[!(is.na(numstages))]
    subx2 = subx[!(is.na(subx))]
    subnumstages2 = subnumstages[!(is.na(subx))]
    sm = sapply(0:4,function(i) median(subx2[subnumstages2 == i]))
    sm
}

lmpvall = apply(subratio3_nttp[rownames(subratio3_nttp) %in% oid3,],1,function(x) summary(lm(x ~ numstages))$coefficients[8])
lmpadjl = p.adjust(lmpvall,'fdr')
stage_ms = t(apply(subratio3_nttp[rownames(subratio3_nttp) %in% oid3,],1,function(x) cal_sm(x)))
#colnames(stage_ms) = c("normal",paste("stage",1:4,sep = ""))
## gradually increase or decrease across stages
#monoup = apply(stage_ms,1,function(x) x[1] < x[2] & x[2] < x[3] & x[3] < x[4] & x[4] < x[5])
#monodn = apply(stage_ms,1,function(x) x[1] > x[2] & x[2] > x[3] & x[3] > x[4] & x[4] > x[5])
#monoupsig = cbind(stage_ms[monoup & lmpadjl < pcut,],lmpvall[monoup & lmpadjl < pcut],lmpadjl[monoup & lmpadjl < pcut])
#monodnsig = cbind(stage_ms[monodn & lmpadjl < pcut,],lmpvall[monodn & lmpadjl < pcut],lmpadjl[monodn & lmpadjl < pcut])
monores = cbind(stage_ms,lmpvall,lmpadjl)
colnames(monores) = c("NT",paste("stage",1:4,sep = ""),"stage_pval","stage_padj")

# dysregulated event at the first stage or metastasis
s1pval = apply(subratio3_nttp[rownames(subratio3_nttp) %in% oid3,],1,function(x) dowilcoxtest(x[numstages == 1],x[numstages == 0]))
s4pval = apply(subratio3_nttp[rownames(subratio3_nttp) %in% oid3,],1,function(x) wilcox.test(x[numstages == 4],x[numstages < 4 & numstages > 0])$p.value)

s1padj = p.adjust(s1pval,'fdr')
s4padj = p.adjust(s4pval,'fdr')

s1s4res = cbind(s1pval,s1padj,s4pval,s4padj)
colnames(s1s4res) = c("stage1_pval","stage1_padj","stage4_pval","stage4_padj")

#lmpval = unlist(lapply(lmcoefl,function(x) x[8]))
## early and advanced marker
eapval = apply(subratio3_nttp[rownames(subratio3_nttp) %in% oid3,],1,function(x) dowilcoxtest(x[numstages == 1 | numstages == 2],x[numstages == 0]))
adpval = apply(subratio3_nttp[rownames(subratio3_nttp) %in% oid3,],1,function(x) dowilcoxtest(x[numstages == 3 | numstages == 4],x[numstages == 0]))
eaadpval = apply(subratio3_nttp[rownames(subratio3_nttp) %in% oid3,],1,function(x) dowilcoxtest(x[numstages == 3 | numstages == 4],x[numstages == 1 | numstages == 2]))

cal_eaadm = function(x){
    subx = x[!(is.na(numstages))]
    subnumstages = numstages[!(is.na(numstages))]
    subx2 = subx[!(is.na(subx))]
    subnumstages2 = subnumstages[!(is.na(subx))]
    eaadm = c(median(subx2[subnumstages2 == 0]),median(subx2[subnumstages2 == 1 | subnumstages2 == 2]),median(subx2[subnumstages2 == 3 | subnumstages2 == 4]))
    eaadm
}
eaadms = t(apply(subratio3_nttp[rownames(subratio3_nttp) %in% oid3,],1,function(x) cal_eaadm(x)))
eaadres = cbind(eaadms,eapval,p.adjust(eapval,'fdr'),adpval,p.adjust(adpval,'fdr'),eaadpval,p.adjust(eaadpval,'fdr'))
colnames(eaadres) = c("NT","TP_early","TP_advanced","early_pval","early_padj","advanced_pval","advanced_padj","advanced_early_pval","advanced_early_padj")
write.table(cbind(monores,eaadres[,-1],s1s4res),file = paste("./txt/",cancertype,".stage.splicing.level.txt",sep = ""),sep = "\t",quote = F)


cat("overlapping 3'UTR splicing events with somatic mutation...\n")
mutdat = read.table(paste("./cbioportal/cbiop_",tolower(cancertype),"/data_mutations_extended.txt",sep = ""),head = T,sep = "\t",quote = "")
submutdat = mutdat[grep("Frame|Mutation|Site|Splice",mutdat[,10]),]
tn = length(intersect(substr(eachsummary$barcode,0,12),substr(submutdat[,17],0,12)))
tsam = intersect(substr(eachsummary$barcode,0,12),substr(submutdat[,17],0,12))
## only keep the sample profiled Mutation ...
oversppatl = lapply(overspsaml,function(x) intersect(x,tsam))
downsppatl = lapply(downspsaml,function(x) intersect(x,tsam))
##oversppatl2 = lapply(overspsaml2,function(x) intersect(x,tsam))
##downsppatl2 = lapply(downspsaml2,function(x) intersect(x,tsam))
ctnnb1os = cbind(tsam,rep("No",length(tsam)))
ctnnb1os[tsam %in% oversppatl[rownames(tmptpr) == "3:41280846-41281309"][[1]],2] = "Yes"
if(cancertype == "LIHC"){
    write.table(ctnnb1os,file = "LIHC.ctnnb1.utr3.splicing.status.txt",quote = F,sep = "\t",row.names = F,col.names = F)
}

##
mutgenel = read.table(paste("./cbioportal/Mutated_Genes_",cancertype,".txt",sep = ""),skip = 1,sep = "\t")
if(sum(!(is.na(mutgenel[,2])) > 0)){
    submutgenel = mutgenel[!(is.na(mutgenel[,2])),]
}else{
    submutgenel = mutgenel[mutgenel[,6] == "Yes",]
}
submutgenel2 = submutgenel[as.numeric(gsub("\\%","",submutgenel[,5]))/100 * tn > 10,]
submutgenel2 = submutgenel2[submutgenel2[,1] %in% mutdat[,1],]

submutgenelpat = sapply(1:nrow(submutgenel2),function(k) unique(substr(submutdat[submutdat[,1] == as.character(submutgenel2[k,1]),17],0,12)))
submutgenelpat2 = lapply(submutgenelpat,function(x) x[x %in% substr(eachsummary$barcode,0,12)])
### l1 should be overspliced sample, implement by distribition
pseudoc = 1
overlapmut = function(l1,l2,output = "p"){
    q = length(intersect(l1,l2))
    m = length(l1)
    n = tn - m
    k = length(l2)
    if(output == "p"){
        if(length(l1) > 10 & length(l1) < tn * 0.9){
            p1 = sum(dhyper(q:min(m,k),m,n,k))
            p2 = phyper(q,m,n,k,lower.tail = TRUE)
            if(q == 0){
                q = q + pseudoc
            }
            if(m - q == 0){
                m = m + pseudoc
            }
            if(k - q == 0){
                k = k + pseudoc
            }
            if(tn - m - k + q == 0){
                tn = tn + pseudoc
            }
            or = ((tn - m - k + q)/(m - q))/((k - q)/q)
            res = paste(p1,p2,or,sep = ":")
        }else{
            res = "NA:NA:NA"
        }
    }else{
        res = paste(q,m,n,k,sep = ":")
    }
    res
}

overspmutpval = t(sapply(1:length(oversppatl),function(k) unlist(lapply(submutgenelpat2,function(x) overlapmut(oversppatl[[k]],x)))))
downspmutpval = t(sapply(1:length(downsppatl),function(k) unlist(lapply(submutgenelpat2,function(x) overlapmut(oversppatl[[k]],x)))))

overspmutspn = t(sapply(1:length(oversppatl),function(k) unlist(lapply(submutgenelpat2,function(x) overlapmut(oversppatl[[k]],x,output = "num")))))
downspmutspn = t(sapply(1:length(downsppatl),function(k) unlist(lapply(submutgenelpat2,function(x) overlapmut(oversppatl[[k]],x,output = "num")))))

colnames(overspmutpval) = colnames(overspmutpval) = as.character(submutgenel2[,1])
rownames(overspmutpval) = rownames(downspmutpval) = paste(rownames(tmptpr),gns,sep = "|")

colnames(overspmutspn) = colnames(overspmutspn) = as.character(submutgenel2[,1])
rownames(overspmutspn) = rownames(downspmutspn) = paste(rownames(tmptpr),gns,sep = "|")

write.table(overspmutpval,file = paste("./txt/",cancertype,".u3spup.overlap.mutation.v2.txt",sep = ""),sep = "\t",quote = F)
write.table(downspmutpval,file = paste("./txt/",cancertype,".u3spdn.overlap.mutation.v2.txt",sep = ""),sep = "\t",quote = F)
write.table(overspmutspn,file = paste("./txt/",cancertype,".u3spup.overlap.mutation.num.v2.txt",sep = ""),sep = "\t",quote = F)
write.table(downspmutspn,file = paste("./txt/",cancertype,".u3spdn.overlap.mutation.num.v2.txt",sep = ""),sep = "\t",quote = F)

cat("correlatting 3'UTR splicing levels with somatic mutation...\n")

## dowilcox test to check associatiob between mutation and splicing
compareR = function(rvec,mutl,tl){
    #subrvec = rvec[substr(eachtpsummary$barcode,0,12) %in% tsam]
    #    subrvec = subrvec[!(is.na(subrvec))]
    mutrvec = rvec[tl %in% mutl]
    wtrvec = rvec[!(tl %in% mutl)]
    if(sum(!(is.na(mutrvec))) > 10 & sum(!(is.na(wtrvec))) > 10){
        p = wilcox.test(mutrvec,wtrvec)$p.value
    }else{
        p = NA
    }
    diff = median(mutrvec[!(is.na(mutrvec))]) - median(wtrvec[!(is.na(wtrvec))])
    res = paste(p,diff,sep = ":")
    res
}

subtmptpr = tmptpr[,colnames(tmptpr) %in% tsam]
#subtmptpr2 = tmptpr2[,colnames(tmptpr2) %in% tsam]

spmutwpval = t(apply(subtmptpr,1,function(x) unlist(lapply(submutgenelpat2,function(y) compareR(x,y,colnames(subtmptpr))))))
rownames(spmutwpval) = paste(rownames(subtmptpr),gns,sep = "|")
colnames(spmutwpval) = as.character(submutgenel2[,1])
write.table(spmutwpval,file = paste("./txt/",cancertype,".u3sp.correlate.mutation.v2.txt",sep = ""),sep = "\t",quote = F)
#write.table(spmutwpval2,file = paste("./txt/",cancertype,".u3sp.correlate.mutation.v2.txt",sep = ""),sep = "\t",quote = F)

cat("saving the data ...\n")
save(list = ls(),file = paste("./Rdata/",cancertype,".tp.vs.nt.v2.Rdata",sep = ""))

### plot example of co-occurence and mutually exclusive
if(cancertype == "COAD"){
    sevent = "11:532523-532630"
    smutge = "TP53"
    eachtpr = unlist(subtmptpr[rownames(subtmptpr) == sevent,])
    mutgn = (1:nrow(submutgenel2))[submutgenel2[,1] == smutge]
    mutgsam = submutgenelpat2[[mutgn]]
}

## implement by fisher test (too slow to run ...)
overlapmut2 = function(l1,l2){
    if(length(l1) > 10 & length(l1) < tn * 0.9){
        n1 = tn - (length(l1) + length(l2) - length(intersect(l1,l2)))
        n2 = length(setdiff(l1,l2))
        n3 = length(setdiff(l2,l1))
        n4 = length(intersect(l2,l1))
        if(n2 == 0 | n3 == 0 | n4 == 0 | n1 == 0){
            n1 = n1 + pseudoc
            n2 = n2 + pseudoc
            n3 = n3 + pseudoc
            n4 = n4 + pseudoc
        }
        f = fisher.test(matrix(c(n1,n2,n3,n4),2),alternative="greater")
        p2 = fisher.test(matrix(c(n1,n2,n3,n4),2),alternative="less")$p.value
        res = paste(f$p.value,p2,f$estimate,sep = ":")
    }else{
        res = "NA:NA:NA"
    }
    res
}


if(FALSE){
    pdf(paste("TCGA.",cancertype,".volcano.analysis.v2.pdf",sep = ""))
    plot(subfres21[,3],-log10(subfres21[,5]),pch = 19,xlim = c(-0.5,0.5),xlab = "splicing changes (TP vs NT)",ylab = "-log10 (FDR)",main = "paired samples")
    points(sigfres21[abs(sigfres21[,3]) > 0.05,3],-log10(sigfres21[abs(sigfres21[,3]) > 0.05,5]),pch = 19,col = cls[2])
    labels = sapply(1:nrow(sigfres21),function(n) unique(as.character(utr3intr[utr3intr[,1] == rownames(sigfres21)[n],7])))
    filtered = abs(sigfres21[,3]) > 0.05 & sigfres21[,5] < 0.01
    text(sigfres21[filtered,3],-log10(sigfres21[filtered,5]) + .1,label = labels[filtered],cex = .6)
    ##
    plot(subfres1[,2] - subfres1[,1],-log10(subfres1[,4]),pch = 19,xlim = c(-0.5,0.5),xlab = "splicing changes (TP vs NT)",ylab = "-log10 (FDR)",main = "all samples")
    points(sigfres1[abs(sigfres1[,2] - sigfres1[,1]) > 0.05,2] - sigfres1[abs(sigfres1[,2] - sigfres1[,1]) > 0.05,1],-log10(sigfres1[abs(sigfres1[,2] - sigfres1[,1]) > 0.05,4]),pch = 19,col = cls[2])
    labels = sapply(1:nrow(sigfres1),function(n) unique(as.character(utr3intr[utr3intr[,1] == rownames(sigfres1)[n],7])))
    filtered = abs(sigfres1[,2] - sigfres1[,1]) > 0.05 & sigfres1[,4] < 0.01
    text(sigfres1[filtered,2] - sigfres1[filtered,1],-log10(sigfres1[filtered,4]) + .1,label = labels[filtered],cex = .6)

    plot(subfres22[,3],-log10(subfres22[,7]),pch = 19,xlim = c(-0.5,0.5),xlab = "splicing changes (TP vs NT)",ylab = "-log10 (FDR)",main = "paired samples")
    points(sigfres22[abs(sigfres22[,3]) > 0.05,3],-log10(sigfres22[abs(sigfres22[,3]) > 0.05,7]),pch = 19,col = cls[2])
    labels = sapply(1:nrow(sigfres22),function(n) unique(as.character(utr3intr[utr3intr[,1] == rownames(sigfres22)[n],7])))
    filtered = abs(sigfres22[,3]) > 0.05 & sigfres22[,5] < 0.01
    text(sigfres22[filtered,3],-log10(sigfres22[filtered,7]) + .1,label = labels[filtered],cex = .6)
    ##
    plot(subfres2[,2] - subfres2[,1],-log10(subfres2[,6]),pch = 19,xlim = c(-0.5,0.5),xlab = "splicing changes (TP vs NT)",ylab = "-log10 (FDR)",main = "all samples")
    points(sigfres2[abs(sigfres2[,2] - sigfres2[,1]) > 0.05,2] - sigfres2[abs(sigfres2[,2] - sigfres2[,1]) > 0.05,1],-log10(sigfres2[abs(sigfres2[,2] - sigfres2[,1]) > 0.05,6]),pch = 19,col = cls[2])
    labels = sapply(1:nrow(sigfres2),function(n) unique(as.character(utr3intr[utr3intr[,1] == rownames(sigfres2)[n],7])))
    filtered = abs(sigfres2[,2] - sigfres2[,1]) > 0.05 & sigfres2[,5] < 0.01
    text(sigfres2[filtered,2] - sigfres2[filtered,1],-log10(sigfres2[filtered,6]) + .1,label = labels[filtered],cex = .6)
    dev.off()
}

### select samples for IGV plot
events = c("CHEK1","CTNNB1","MAPK1","THUMPD1","WDR55")
eventsid = c("11:125525243-125526100","3:41280846-41281309","22:22118530-22123483","16:20747681-20748181","5:140049276-140050853")
ctnnb1_tpr = unlist(subratio3_tp[rownames(subratio3_tp) == eventsid[2],])
ctnnb1_ntr = unlist(subratio3_nt[rownames(subratio3_nt) == eventsid[2],])
ctnnb1_tpr = ctnnb1_tpr[!(is.na(ctnnb1_tpr))]
ctnnb1_ntr = ctnnb1_ntr[!(is.na(ctnnb1_ntr))]
# top 3 over-splcied tumor samples
ctnnb1_tpr = ctnnb1_tpr[order(ctnnb1_tpr,decreasing = T)]
ctnnb1_ntr = ctnnb1_ntr[order(names(ctnnb1_ntr))]
ctnnb1_rsum = unlist(subrsums2[rownames(subrsums2) == eventsid[2],])
ctnnb1_tpr2 = ctnnb1_tpr[names(ctnnb1_tpr) %in% names(ctnnb1_rsum)[ctnnb1_rsum > 200]]
ctnnb1_ntr2 = ctnnb1_ntr[names(ctnnb1_ntr) %in% names(ctnnb1_rsum)[ctnnb1_rsum > 200]]
write.table(cbind(names(ctnnb1_tpr2)[1:3],names(ctnnb1_ntr2)[1:3]),file = paste("./txt/",cancertype,".ctnnb1.selected.event.samples.txt",sep = ""),sep = "\t",quote = F,row.names = F,col.names = F)


if(FALSE){
    if(cancertype == "LIHC"){

        eventstpsam = vector()
        eventsntsam = vector()
        for(i in 1:length(events)){
            eachtpratio = unlist(subratio3_tp[rownames(subratio3_tp) == eventsid[i],])
            eachtprsums = unlist(subrsums2_tp[rownames(subrsums2_tp) == eventsid[i],])
            eachntratio = unlist(subratio3_nt[rownames(subratio3_nt) == eventsid[i],])
            eachntrsums = unlist(subrsums2_nt[rownames(subrsums2_nt) == eventsid[i],])
            eachtprsums = eachtprsums[order(eachtprsums,decreasing = T)]
            eachntrsums = eachntrsums[order(eachntrsums,decreasing = T)]
            eachtpratio = eachtpratio[!(is.na(eachtpratio))]
            eachntratio = eachntratio[!(is.na(eachntratio))]
            etpsam = c(names(eachtprsums)[names(eachtprsums) %in% names(eachtpratio)[eachtpratio == max(eachtpratio)]][1],names(eachtprsums)[names(eachtprsums) %in% names(eachtpratio)[eachtpratio == min(eachtpratio)]][1])
            entsam = c(names(eachntrsums)[names(eachntrsums) %in% names(eachntratio)[eachntratio == max(eachntratio)]][1],names(eachntrsums)[names(eachntrsums) %in% names(eachntratio)[eachntratio == min(eachntratio)]][1])
            eventstpsam = c(eventstpsam,etpsam)
            eventsntsam = c(eventsntsam,entsam)
        }
        write.table(cbind(eventstpsam,eventsntsam),file = "./txt/LIHC.selected.event.samples.txt",sep = "\t",quote = F,row.names = F,col.names = F)
    }

}
### event list

