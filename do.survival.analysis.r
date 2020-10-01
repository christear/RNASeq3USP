library(ggpubr)
#### load survival analysis functions
cat("loading survival analysis functions ...\n")

##
pcut = 0.1
##

argv=commandArgs(TRUE)
cancertype = as.character(argv[1])
cls = get_palette("lancet",7)
##
cat("loading the data from",cancertype,"...\n")
load(paste("./Rdata/",cancertype,".tp.vs.nt.v2.Rdata",sep = ""))
source("survive.analysis.functions.v2.r")
###
cat("renaming stages and adding splicing event numbers...\n")
cstage = rep(NA,nrow(clinicalout))
#cstage[cstage == ""] = "unclear"
cstage[grep("stage 0",clinicalout$stage)] = 0
cstage[grep("stage i[abcd]$|stage i$",clinicalout$stage)] = 1
cstage[grep("stage ii[abcd]$|stage ii$",clinicalout$stage)] = 2
cstage[grep("stage iii[abcd]$|stage iii",clinicalout$stage)] = 3
cstage[grep("stage iv[abcd]$|stage iv",clinicalout$stage)] = 4

clinicalout2 = cbind(clinicalout,cstage)
uispnum = spgenesnum3[,2]
names(uispnum) = eachsummary$barcode
uispnum_tp = uispnum[names(uispnum) %in% eachsummary$barcode[eachsummary$sample_type == "TP"]]
alspnum = allspn
names(alspnum) = eachsummary$barcode
alspnum_tp = alspnum[names(alspnum) %in% eachsummary$barcode[eachsummary$sample_type == "TP"]]
uispnum_tp = uispnum_tp[order(names(uispnum_tp))]
alspnum_tp = alspnum_tp[order(names(alspnum_tp))]

subuispnum_tp = sapply(1:nrow(clinicalout),function(i) mean(uispnum_tp[substr(names(uispnum_tp),0,12) == rownames(clinicalout)[i]]))
subalspnum_tp = sapply(1:nrow(clinicalout),function(i) mean(alspnum_tp[substr(names(alspnum_tp),0,12) == rownames(clinicalout)[i]]))
clinicalout2$utr3spnum = subuispnum_tp
clinicalout2$totalspnum = subalspnum_tp

### if the input is a matrix/data.frame, return vector of the mean of each row
getmean = function(x){
    if(is.null(dim(x))){
        mx = x
    }else{
        mx = apply(x,1,function(x) mean(x[!(is.na(x))]))
    }
    mx
}

eachtpsummary = eachsummary[eachsummary$sample_type == "TP",]
eachntsummary = eachsummary[eachsummary$sample_type == "NT",]

if(ncol(subratio3_tp) == sum(colnames(subratio3_tp) == eachtpsummary$analysis_id)){
    patientid = intersect(rownames(clinicalout),substr(eachtpsummary$barcode,0,12))
    subratio4_tp = sapply(1:length(patientid),function(i) getmean(subratio3_tp[,substr(eachtpsummary$barcode,0,12) == patientid[i]]))
    rownames(subratio4_tp) = rownames(subratio3_tp)
    colnames(subratio4_tp) = patientid
}else{
    cat("error:sample ids don't match\n")
}

#if(ncol(subgexpda2_tp) == sum(colnames(subgexpda2_tp) == eachtpsummary$analysis_id)){
#    subgexpda3_tp = sapply(1:length(patientid),function(i) getmean(subgexpda2_tp[,substr(eachtpsummary$barcode,0,12) == patientid[i]]))
#    rownames(subgexpda3_tp) = rownames(subgexpda2_tp)
#    colnames(subgexpda3_tp) = patientid
#}

#gn = rownames(subgexpda)
#oid4 = intersect(oid3,utr3intr[utr3intr[,7] %in% gn,1])
#oid4 = intersect(oid4,rownames(subratio4_tp)[apply(subratio4_tp,1,function(x) sd(x[!(is.na(x))])) > 0])
oid4 = intersect(oid3,rownames(subratio4_tp)[apply(subratio4_tp,1,function(x) sd(x[!(is.na(x))])) > 0])
subratio4_tp = subratio4_tp[rownames(subratio4_tp) %in% oid4,]
### tp samples are collapsed while nt samples don't need
subratio4_nt = subratio3_nt[rownames(subratio3_nt) %in% oid4,]
if(sum(colnames(subratio4_nt) == eachntsummary$analysis_id) == ncol(subratio4_nt)){
    colnames(subratio4_nt) = substr(eachntsummary$barcode,0,12)
}

subclinicalout2 = clinicalout2[rownames(clinicalout2) %in% patientid,]
tmp = as.matrix(cbind(subclinicalout2,t(subratio4_tp)))
write.table(tmp,file = paste("./txt/",cancertype,".data4surv.txt",sep = ""),sep = "\t",quote = F)

docoxph = function(tvec,dvec,x){
    subt = tvec[!(is.na(x))]
    subd = dvec[!(is.na(x))]
    subx = x[!(is.na(x))]
    subx2 = subx[subx > 0 & subx < 1]
    subt2 = subt[subx > 0 & subx < 1]
    subd2 = subd[subx > 0 & subx < 1]
    if(max(subx2) - min(subx2) < 0.1 | sum(subd2) < 2){
        cfit = NA
    }else{
        cfit = coxph(Surv(subt2,subd2) ~ subx2)
    }
    cfit
}

dokmsurv = function(tvec,dvec,x){
    subt = tvec[!(is.na(x))]
    subd = dvec[!(is.na(x))]
    subx = x[!(is.na(x))]
    resl = list()
    percetile = 20:80/100
    for(j in 1:length(percetile)){
        segv = quantile(subx,percetile[j])
        if(segv == max(subx)){
            resl[[j]] = dosurv(subt,subd,subx,segv = segv)
            break
        }else{
            resl[[j]] = dosurv(subt,subd,subx,segv = segv)
        }
    }
    pval = unlist(lapply(resl,function(x) 1 - pchisq(x[[2]]$chisq, length(x[[2]]$n) - 1)))
    hdif = unlist(lapply(resl,function(x) x[[2]]$obs[1] - x[[2]]$exp[1]))
    hn = unlist(lapply(resl,function(x) x[[2]]$n[1]))
    ln = unlist(lapply(resl,function(x) x[[2]]$n[2]))
    ## sample size in either group is too small (<20% of total avaible sample), will define p equal to 1 ...
    pval[hn < length(subx)/20 | ln < length(subx)/20] = 1
    ##
    fitl = lapply(resl,function(x) x[[1]])
    finalres = list(fitl[pval == min(pval)],min(percetile[pval == min(pval)]),min(pval),mean(hdif[pval == min(pval)]),min(hn[pval == min(pval)]),max(ln[pval == min(pval)]))
    finalres
}

dokmsurv2 = function(tvec,dvec,x,quan = 1/2){
    subt = tvec[!(is.na(x))]
    subd = dvec[!(is.na(x))]
    subx = x[!(is.na(x))]
    subx2 = subx[subx > 0 & subx < 1]
    subt2 = subt[subx > 0 & subx < 1]
    subd2 = subd[subx > 0 & subx < 1]
    if(sum(subx2 <= quantile(subx2,quan)) == length(subx2) | max(subx2) - min(subx2) < 0.1){
        res = NA
    }else{
        res = dosurv(subt2,subd2,subx2,quan = quan)
    }
    res
}


cat("doing survival analysis on splicing levels ...\n")
### 3'UTR splicing levels
u3spcoxphl = sapply(1:nrow(subratio4_tp),function(i) docoxph(subclinicalout$new_death,subclinicalout$death_event,unlist(subratio4_tp[i,])),simplify = F)
u3spcoxphpval = unlist(lapply(u3spcoxphl,function(x) ifelse(length(x) == 1,NA,summary(x)[[12]][[3]])))
u3spcoxphhr = unlist(lapply(u3spcoxphl,function(x) ifelse(length(x) == 1,NA,summary(x)$coef[2])))
u3spcoxphres = cbind(u3spcoxphpval,u3spcoxphhr)
gns = sapply(1:nrow(subratio4_tp),function(i) paste0(unique(utr3intr[utr3intr[,1] == rownames(subratio4_tp)[i],7]),collapse = ":"))
rownames(u3spcoxphres) = paste(rownames(subratio4_tp),gns,sep = "|")
colnames(u3spcoxphres) = c("pvalue","HR")
write.table(u3spcoxphres,file = paste("./txt/",cancertype,".utr3sp.coxph.result.txt",sep = ""),quote = F,sep = "\t")

### K-M
#u3spkmsurvl = sapply(1:10,function(i) dokmsurv(subclinicalout$new_death,subclinicalout$death_event,unlist(subratio4_tp[i,])),simplify = F)
u3spkmsurvl = sapply(1:nrow(subratio4_tp),function(i) dokmsurv(subclinicalout$new_death,subclinicalout$death_event,unlist(subratio4_tp[i,])),simplify = F)
u3spkmsurvl2 = sapply(1:nrow(subratio4_tp),function(i) dokmsurv2(subclinicalout$new_death,subclinicalout$death_event,unlist(subratio4_tp[i,]),quan = 1/2),simplify = F)

u3spkmpval = unlist(lapply(u3spkmsurvl,function(x) x[[3]]))
u3spkmpval2 = unlist(lapply(u3spkmsurvl2,function(x) ifelse(length(x) == 1,NA,1 - pchisq(x[[2]]$chisq,length(x[[2]]$n) - 1))))
u3spkmhdif = unlist(lapply(u3spkmsurvl,function(x) x[[4]]))
u3spkmhdif2 = unlist(lapply(u3spkmsurvl2,function(x) ifelse(length(x) == 1,NA,x[[2]]$obs[1] - x[[2]]$exp[1])))
u3spkmhn = unlist(lapply(u3spkmsurvl,function(x) x[[5]]))
u3spkmhn2 = unlist(lapply(u3spkmsurvl2,function(x) ifelse(length(x) == 1,NA,x[[2]]$n[1])))
u3spkmln = unlist(lapply(u3spkmsurvl,function(x) x[[6]]))
u3spkmln2 = unlist(lapply(u3spkmsurvl2,function(x) ifelse(length(x) == 1,NA,x[[2]]$n[2])))


cex = .8
lwd = 1.5
pdf(paste("TCGA.",cancertype,".survival.analysis.v2.pdf",sep = ""))
par(mfrow = c(2,2))
for(k in 1:nrow(subratio4_tp)){
    if(k %% 100 == 0){
        cat("processed K-M plot for",k," events ...\n")
    }
    eachkmsurv = u3spkmsurvl[[k]]
    if(eachkmsurv[[3]] < pcut/2){
        plot(eachkmsurv[[1]][[1]],col = cls[2:1],main = paste(rownames(subratio4_tp)[k],"\n",gns[k],":",eachkmsurv[[4]],sep = ""),cex.main = cex)
        eachspr = unlist(subratio4_tp[k,])
        eachspr = eachspr[!(is.na(eachspr))]
        legend("topright",col = cls[2:1],lty = c(1,1,NA),legend = c(paste(c("H=","L="),c(sum(eachspr > quantile(eachspr,eachkmsurv[[2]])),sum(eachspr <= quantile(eachspr,eachkmsurv[[2]])))),eachkmsurv[[3]]),lwd = lwd,cex = cex)
    }
}
dev.off()

kmsurvres = cbind(u3spkmpval,u3spkmhdif,u3spkmhn,u3spkmln)
kmsurvres2 = cbind(u3spkmpval2,u3spkmhdif2,u3spkmhn2,u3spkmln2)
#kmsurvres2 = kmsurvres2[!(is.na(kmsurvres2[,1])),]
rownames(kmsurvres) = rownames(kmsurvres2) = paste(rownames(subratio4_tp),gns,sep = "|")
colnames(kmsurvres) = colnames(kmsurvres2) = c("pvalue","H_ob-exp","high","low")
write.table(kmsurvres,file = paste("./txt/",cancertype,".utr3sp.kmsurv.result.txt",sep = ""),quote = F,sep = "\t")
write.table(kmsurvres2[!(is.na(kmsurvres2[,1])),],file = paste("./txt/",cancertype,".utr3sp.kmsurv.result.v2.txt",sep = ""),quote = F,sep = "\t")

#save(list = ls()[grep("surv|km|cox|clini|ratio")],file = paste("./Rdata/",cancertype,".surv.v2.Rdata",sep = ""))






