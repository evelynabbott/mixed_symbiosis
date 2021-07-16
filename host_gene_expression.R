#coral host gene expression
library(DESeq2)
library(vegan)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
theme_set(theme_cowplot())
library(data.table)
library(limma)
library(viridis)
library(hexbin)

#deseq================================================================================================
rm(list=ls())
setwd("~/Dropbox/codominant_symbiosis/")

ll=load('counts_metadata/coldata.Rdata') 

# calculating symbiont proportions
cc=10^(coldata$logCcount)
dd=10^(coldata$logDcount)
coldata$Cprop=cc/(cc+dd)
coldata$Dprop=1-coldata$Cprop

#remove  Barshis outliers 
outliers=c("SRR629129",
"SRR629150",
"SRR629152",
"SRR629133",
"SRR629130",
"SRR629151",
"SRR629134",
"SRR629153",
"SRR629155")
coldata=coldata[!(coldata$Run %in% outliers),]
row.names(coldata)=coldata$Run

# loading counts, aligning them with coldata
ll=load('counts_metadata/host_counts')
rownames(coldata)=coldata$Run
keeps = row.names(coldata)
counts = counts[,keeps]

# normalizing "pool" factors
coldata$pool[coldata$pool=="300"]="HV"
coldata$pool[coldata$pool=="400"]="MV"

# removing low-count genes
cut=3
means=apply(counts,1,mean)
table(means>=cut) #how many genes have mean below the cutoff?
# FALSE  TRUE 
# 593 18987 
counts=counts[means>=cut,]

coldata$hostCount=apply(counts,2,sum)
coldata$symCount=10^coldata$logCcount+10^coldata$logDcount
coldata$symRatio=coldata$symCount/coldata$hostCount
coldata$logSymRatio=log(coldata$symRatio,10)

pdf("logCDratio_vs_logSymRatio.pdf",width=4, height=4)
plot(logCDratio~logSymRatio,coldata,xlab="log(symbiont counts/host counts)",ylab="log(C:D ratio)",mgp=c(2.3,1,0))
dev.off()

# ----- supplemental figures of D proportion by-genotype etc

coldata$study="Seneca and Palumbi 2015"
coldata$study[coldata$my_title=="L_Barshis_bleachResillience_PRJNA177515"]="Barshis et al 2013"
coldata$treatment="heat"
coldata$treatment[coldata$c_t=="c"]="control"
pdf("byTreatment.pdf",width=8,height=2.7)
ggplot(coldata,aes(geno,Dprop,color=treatment))+
  geom_boxplot()+
  ylab("Durusdinium proportion")+
  xlab("source colony")+
  theme(legend.position="top")
dev.off()
pdf("byStudy.pdf",width=8,height=2.7)
ggplot(coldata,aes(geno,Dprop,color=study))+geom_boxplot()+ylab("Durusdinium proportion")+xlab("source colony")+  theme(legend.position="top")
dev.off()
pdf("byPool.pdf",width=8,height=2.7)
ggplot(coldata,aes(geno,Dprop,color=pool))+geom_boxplot()+ylab("Durusdinium proportion")+xlab("source colony")+  theme(legend.position="top")
dev.off()

#----------------------------

# resampling data (samples) for all experimental groups to have equal number (10) of mixed, all-D and all-C samples

# defining groups to resample
evens.h=coldata$Run[coldata$Dprop>0.1 & coldata$Dprop<0.9 & coldata$c_t=="h"]
evens.c=coldata$Run[coldata$Dprop>0.1 & coldata$Dprop<0.9 & coldata$c_t=="c"]
Ds.h=coldata$Run[coldata$Dprop>0.9 & coldata$c_t=="h"]
Ds.c=coldata$Run[coldata$Dprop>0.9 & coldata$c_t=="c"]
Cs.h=coldata$Run[coldata$Cprop>0.9 & coldata$c_t=="h"]
Cs.c=coldata$Run[coldata$Cprop>0.9 & coldata$c_t=="c"]
length(evens.h) # 17
length(evens.c) # 15
length(Ds.c) # 50
length(Ds.h) # 48
length(Cs.c) # 22
length(Cs.h) # 20

# Dixon red module
ll=load("moduleAssignment.Rdata")
table(moduleColors) # 634 genes in the original red module
# extracting red module genes (only kME>0.25)
reds=row.names(geneModuleMembership[moduleColors=="red" & geneModuleMembership$MMred>0.25,])
reds=intersect(reds,row.names(counts))
length(reds)
# 621
# calculating red eigengene across all samples (to track its sign during resampling)
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~pool+time)
vsd=assay(vst(dds))

# pca before batch effect removal
rawpca=rda(t(vsd)~1)
coldata$rawPC1=scores(rawpca,choices=1,dis="sites",scaling="sites")
coldata$rawPC2=scores(rawpca,choices=2,dis="sites",scaling="sites")

# pca after batch effect removal
vsd=removeBatchEffect(vsd,batch=coldata$pool,batch2=coldata$time)
pca1=rda(t(vsd)~1)
coldata$PC1=scores(pca1,choices=1,dis="sites",scaling="sites")
coldata$PC2=scores(pca1,choices=2,dis="sites",scaling="sites")

coldata$mixed="no"
coldata$mixed[coldata$Dprop>0.1 & coldata$Dprop<0.9]="yes"
coldata$treatment=factor(coldata$treatment,levels=c("heat","control"))

# plotting PCA versions
ggplot(coldata,aes(rawPC1,rawPC2))+geom_point(aes(shape=treatment, color=study))+coord_equal()+scale_shape_manual(values=c(1,19))
ggplot(coldata,aes(PC1,PC2))+geom_point(aes(shape=treatment, color=study))+coord_equal()+scale_shape_manual(values=c(1,19))
ggplot(coldata,aes(PC1,PC2))+geom_point(aes(color=treatment, shape=mixed))+coord_equal()+scale_shape_manual(values=c(1,19))

# red module eigengene and PCA
red.vsd=vsd[reds,]
redpc=rda(t(red.vsd)~1)
# importance of axes:
plot(sqrt(redpc$CA$eig/sum(redpc$CA$eig)))
plot(redpc,main="RED",scaling="sites")
coldata$redpc=scores(redpc,choices=1,dis="sites",scaling="sites")
# plotting red module PCA with different coloring
coldata1=cbind(coldata,scores=scores(redpc,dis="sites",scaling="sites"))
ggplot(coldata1,aes(scores.PC1,scores.PC2,color=c_t))+geom_point()+coord_equal()
ggplot(coldata1,aes(scores.PC1,scores.PC2,color=pool))+geom_point()+coord_equal()
ggplot(coldata1,aes(scores.PC1,scores.PC2,color=time))+geom_point()+coord_equal()
ggplot(coldata1,aes(scores.PC1,scores.PC2,color=Dprop))+geom_point()+coord_equal()

# loading rose modules to examine correlation with the red module
load("rosemodules.RData")
rosemodules$GSR=coldata[row.names(rosemodules),"redpc"]
pdf("pairs_rose_dixon_modules.pdf",width=5,height=5)
pairs(rosemodules)
dev.off()

# resampling and calculating
nreps=100
nsamp=10
deseq=1 # set to 0 if you just want to recalculate Fig 2 data (boxplots of module expression)
cdd=cdcd=cdcdd=cdd.c=cdcd.c=cdcdd.c=rep(0,nrow(counts))
redm.c=redm.h=c()
for ( i in 1:nreps) {
  message(i)
  eec=as.character(sample(evens.c,nsamp))
  ddc=as.character(sample(Ds.c,nsamp))
  ccc=as.character(sample(Cs.c,nsamp))
  eeh=as.character(sample(evens.h,nsamp))
  cch=as.character(sample(Cs.h,nsamp))
  ddh=as.character(sample(Ds.h,nsamp))
  meta=coldata[c(eec,eeh,ccc,cch,ddc,ddh),]
  meta$CD=c(rep("CD",length(eec)+length(eeh)),rep("C",length(ccc)+length(cch)),rep("D",length(ddc)+length(ddh)))
  cts=counts[,c(eec,eeh,ccc,cch,ddc,ddh)]
  ddscd = DESeqDataSetFromMatrix(countData=cts, colData=meta, design=~pool+time)
  vsd=assay(vst(ddscd))
  vsd=removeBatchEffect(vsd,batch=meta$pool,batch2=meta$time)

  if(i==1) { vsda=vsd/nreps} else { vsda=vsda+vsd/nreps }
  # red module eigengene expression
  red.vsd=vsd[reds,]
  redpc=rda(t(red.vsd)~1)
  redscores=scores(redpc,choices=1,dis="sites",scaling="sites")

# flipping PC1 axis if does not match dataset-wide direction
  if(coef(lm(meta$redpc~redscores))[2]<0) { redscores=redscores*(-1)}
  meta$redpc=redscores

# recording linear model coefficients (separately for heated and control)
  redm.c=rbind(redm.c,coef(lm(redpc~CD,subset(meta,c_t=="c"))))
  redm.h=rbind(redm.h,coef(lm(redpc~CD,subset(meta,c_t=="h"))))

# DESeq for heated samples
  if (deseq==1) {
     cts.h=cts[,meta$c_t=="h"]
     meta.h=meta[meta$c_t=="h",]
     ddscd = DESeqDataSetFromMatrix(countData=cts.h, colData=meta.h, design=~CD+pool+time)
     ddscd=DESeq(ddscd)
     cdd=cdd+results(ddscd,name="CD_D_vs_C")$log2FoldChange/nreps
     cdcd=cdcd+results(ddscd,name="CD_CD_vs_C")$log2FoldChange/nreps
     cdcdd=cdcdd+results(ddscd,contrast=c("CD","CD","D"))$log2FoldChange/nreps

# Deseq for control samples

      cts.c=cts[,meta$c_t=="c"]
      meta.c=meta[meta$c_t=="c",]
      ddscd = DESeqDataSetFromMatrix(countData=cts.c, colData=meta.c, design=~CD+pool+time)
      ddscd=DESeq(ddscd)
      cdd.c=cdd.c+results(ddscd,name="CD_D_vs_C")$log2FoldChange/nreps
      cdcd.c=cdcd.c+results(ddscd,name="CD_CD_vs_C")$log2FoldChange/nreps
      cdcdd.c=cdcdd.c+results(ddscd,contrast=c("CD","CD","D"))$log2FoldChange/nreps
  }
}
if(deseq==1) {
  resps=data.frame(cbind(cdd,cdcd,cdcdd,cdd.c,cdcd.c,cdcdd.c))
  row.names(resps)=row.names(counts)
  save(resps,redm.c,redm.h,file="resampled_host_CD_log2fc_cleaned.RData")
}

# assembling data table with lm-predicted expression of the red module
redm.c=data.frame(redm.c)
redm.h=data.frame(redm.h)
redm.c$CDCD=redm.c[,1]+redm.c$CDCD
redm.c$CDD=redm.c[,1]+redm.c$CDD
names(redm.c)=c("C","C+D","D")
redm.h$CDCD=redm.h[,1]+redm.h$CDCD
redm.h$CDD=redm.h[,1]+redm.h$CDD
names(redm.h)=c("C","C+D","D")
rh=stack(redm.h)
rh$treatment="heat"
rc=stack(redm.c)
rc$treatment="control"
rch=data.frame(rbind(rc,rh))
rch$treatment=factor(rch$treatment,levels=c("heat","control"))

# computing jackknifing support
sum(redm.h$C>redm.h$D)
# 92
sum(redm.h[,"C+D"]<redm.h$D)
# 80
sum(redm.h[,"C+D"]<redm.h$C)
# 100

sum(redm.c$C<redm.c$D)
# 86
sum(redm.c[,"C+D"]<redm.c$D)
# 83
sum(redm.c[,"C+D"]>redm.c$C)
# 58


# plotting and Tukey-testing
gg=ggplot(rch,aes(ind,values,color=treatment))+
  geom_boxplot()+
#  geom_point(position=position_jitterdodge(),shape=1)+
#  geom_boxplot(outlier.shape = NA)+
  xlab("symbiont dominance")+
  ylab("expression")+
  ggtitle("GSR module")
pdf("redmod_resampled_30_10_cleaned2_scaled_100.pdf",width=3.35,height=3)
gg
dev.off()

# ---- Dixon blobs

ll=load("resampled_host_CD_log2fc_cleaned.RData")
#r1=stack(resps)
library(viridis)
head(resps)
summary(lm(cdcd~cdd,resps))
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 0.010767   0.003148   3.421 0.000626 ***
#   cdd         0.610310   0.006666  91.550  < 2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4337 on 18989 degrees of freedom
# Multiple R-squared:  0.3062,	Adjusted R-squared:  0.3062 

summary(lm(cdcd.c~cdd.c,resps))
# Coefficients:
# Estimate Std. Error t value Pr(>|t|)    
# (Intercept) -0.030282   0.003217  -9.414   <2e-16 ***
#   cdd.c        0.492329   0.007605  64.737   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Residual standard error: 0.4422 on 18989 degrees of freedom
# Multiple R-squared:  0.1808,	Adjusted R-squared:  0.1808 
# F-statistic:  4191 on 1 and 18989 DF,  p-value: < 2.2e-16

# response to C+D against response to D
pdf("cdcd_vs_cdd_cleaned2+10r.pdf",height=3.5,width=4.2)
gg=ggplot(resps,aes(cdd,cdcd))+
  geom_hex(bins=30)+
  scale_fill_viridis(trans="log")+
  coord_equal()+
#  geom_smooth(method="lm",color="red")+
  geom_abline(intercept=0,slope=1,color="grey80")+
  geom_abline(intercept=0,slope=-1,color="grey80")+
  ylim(-5,5)+xlim(-5,5)+
  geom_vline(xintercept=0,linetype="dotted")+
  geom_hline(yintercept=0,linetype="dotted")+
  xlab("all-D vs all-C")+
  ylab("C+D vs all-C")
gg
dev.off()

# same under control condition
pdf("cdcd.c_vs_cdd.c_cleaned2_10r.pdf",height=3.5,width=4.2)
gg=ggplot(resps,aes(cdd.c,cdcd.c))+
  geom_hex(bins=30)+
  scale_fill_viridis(trans="log")+
  coord_equal()+
  #  geom_smooth(method="lm",color="red")+
  geom_abline(intercept=0,slope=1,color="grey80")+
  geom_abline(intercept=0,slope=-1,color="grey80")+
  ylim(-5,5)+xlim(-5,5)+
  geom_vline(xintercept=0,linetype="dotted")+
  geom_hline(yintercept=0,linetype="dotted")+
  xlab("all-D vs all-C")+
  ylab("C+D vs all-C")
gg
dev.off()

# ------ writing input tables for GO_MWU

write.csv(data.frame(cbind(gene=row.names(resps),lfc=resps$cdd)),file="~/Dropbox/mega2019/mega2019_clean/TagSeq/GO_MWU/cdd_10r.csv",row.names=F, quote=F)
write.csv(data.frame(cbind(gene=row.names(resps),lfc=resps$cdcd)),file="~/Dropbox/mega2019/mega2019_clean/TagSeq/GO_MWU/cdcd_10r.csv",row.names=F, quote=F)
write.csv(data.frame(cbind(gene=row.names(resps),lfc=resps$cdcdd)),file="~/Dropbox/mega2019/mega2019_clean/TagSeq/GO_MWU/cdcdd_10r.csv",row.names=F, quote=F)
write.csv(data.frame(cbind(gene=row.names(resps),lfc=resps$cdcdd.c)),file="~/Dropbox/mega2019/mega2019_clean/TagSeq/GO_MWU/cdcdd.c_10r.csv",row.names=F, quote=F)

