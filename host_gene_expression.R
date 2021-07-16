#coral host gene expression
library(DESeq2)
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

coldata0=read.csv("seneca2015_fullRTEdesign.csv")
coldata0$sample=sub("_.+","",coldata0$X)
coldata0$X=sub("^\\d+_","",coldata0$X)

# loading our data table, recalculating C and D proportions
ll=load('counts_metadata/coldata.Rdata') 
dim(coldata)
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

# this is for matching with Seneca data
coldata$X=paste(coldata$geno,coldata$pool,coldata$c_t,coldata$time,sep="_")

coldata01=merge(coldata0,coldata,by="X",sort=F)

# loading and formatting counts table
counts0=read.csv("seneca2015_Raw_counts_data.csv")
annot=counts0[,154:172]
counts0=counts0[,-c(154:172)]
genes=counts0$ContigName
counts0$ContigName=NULL
row.names(counts0)=genes

names(counts0)=sub("_.+","",names(counts0))
row.names(coldata01)=paste("X",coldata01$sample,sep="")

# aligning Seneca data
goods=intersect(names(counts0), row.names(coldata01))
length(goods)
counts=counts0[,goods]
coldata=coldata01[goods,]
row.names(coldata)=coldata$Run
names(counts)=coldata$Run

# removing low-counts genes
cut=3
means=apply(counts,1,mean)
table(means>=cut) #how many genes have mean below the cutoff?
# FALSE  TRUE 
# 14928 18568
counts=counts[means>=cut,]

# defininng experimental groups for resampling
evens.h=coldata$Run[coldata$Dprop>0.1 & coldata$Dprop<0.9 & coldata$c_t=="h"]
evens.c=coldata$Run[coldata$Dprop>0.1 & coldata$Dprop<0.9 & coldata$c_t=="c"]
Ds.h=coldata$Run[coldata$Dprop>0.9 & coldata$c_t=="h"]
Ds.c=coldata$Run[coldata$Dprop>0.9 & coldata$c_t=="c"]
Cs.h=coldata$Run[coldata$Cprop>0.9 & coldata$c_t=="h"]
Cs.c=coldata$Run[coldata$Cprop>0.9 & coldata$c_t=="c"]
length(evens.h) # 16
length(evens.c) # 13
length(Ds.c) # 44
length(Ds.h) # 43
length(Cs.c) # 19
length(Cs.h) # 17

# loading modules from Noah Rose: 10 and 12 (associated with bleahcing outcome)
# module 10: negative assoc with bleach score
# module 12: positive assoc with bleach score
rosemod=read.table("rose2015_moduleAssignments.txt",sep="\t",header=T)
rosemod10=rosemod$contig[rosemod$module==10]
length(rosemod10) # 181
rosemod12=rosemod$contig[rosemod$module==12]
length(rosemod12) # 201
coldata$timepoint=as.factor(coldata$timepoint)

# computing module eigengenes for the whole dataset (only to keep track of the eigengene's sign during resampling)
dds = DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~c_t*Dprop+pool+timepoint)
vsd=assay(vst(dds))
vsd=removeBatchEffect(vsd,batch=coldata$pool,batch2=coldata$timepoint)
# mod10 eigengene expression
m10.vsd=vsd[rosemod10,]
m10pc=rda(t(m10.vsd)~1)
#plot(m10pc,main="m10")
coldata$m10pc=scores(m10pc,choices=1,dis="sites")
# mod12 eigengene expression
m12.vsd=vsd[rosemod12,]
m12pc=rda(t(m12.vsd)~1)
#plot(m12pc,main="m12")
coldata$m12pc=scores(m12pc,choices=1,dis="sites")
rosemodules=coldata[,c("m10pc","m12pc")]
names(rosemodules)=c("m10","m12")
save(rosemodules,file="rosemodules.RData")

# resampling
nreps=100
nsamp=10
coldata$id=paste(coldata$geno,coldata$pool,sep=".")

m10.c=m10.h=m12.c=m12.h=c()
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
  dds = DESeqDataSetFromMatrix(countData=cts, colData=meta, design=~c_t*Dprop+pool+timepoint)
  vsd=assay(vst(dds))
  vsd=removeBatchEffect(vsd,batch=meta$pool,batch2=meta$timepoint)
  #  str(vsd)
  
  # mod10 eigengene expression
  m10.vsd=vsd[rosemod10,]
  m10pc=rda(t(m10.vsd)~1)
  m10scores=scores(m10pc,choices=1,dis="sites",scaling="sites")
  # flipping axis if needed
  if(coef(lm(meta$m10pc~m10scores))[2]<0) { m10scores=m10scores*(-1)}
  meta$m10pc=m10scores  
  m10.c=rbind(m10.c,coef(lm(m10pc~CD,subset(meta,c_t=="c"))))
  m10.h=rbind(m10.h,coef(lm(m10pc~CD,subset(meta,c_t=="h"))))

  # mod12 eigengene expression
  m12.vsd=vsd[rosemod12,]
  m12pc=rda(t(m12.vsd)~1)
  m12scores=scores(m12pc,choices=1,dis="sites",scaling="sites")
  # flipping axis if needed
  if(coef(lm(meta$m12pc~m12scores))[2]<0) { m12scores=m12scores*(-1)}
  meta$m12pc=m12scores  
  m12.c=rbind(m12.c,coef(lm(m12pc~CD,subset(meta,c_t=="c"))))
  m12.h=rbind(m12.h,coef(lm(m12pc~CD,subset(meta,c_t=="h"))))
  
}

# putting together data frames of lm-predicted values
m10.c=data.frame(m10.c)
m10.h=data.frame(m10.h)
m10.c$CDCD=m10.c[,1]+m10.c$CDCD
m10.c$CDD=m10.c[,1]+m10.c$CDD
names(m10.c)=c("C","C+D","D")
m10.h$CDCD=m10.h[,1]+m10.h$CDCD
m10.h$CDD=m10.h[,1]+m10.h$CDD
names(m10.h)=c("C","C+D","D")
m10h=stack(m10.h)
m10h$treatment="heat"
m10c=stack(m10.c)
m10c$treatment="control"
m10ch=data.frame(rbind(m10c,m10h))
m10ch$treatment=factor(m10ch$treatment,levels=c("heat","control"))

#---- module 10 plotting and jackknifing support 

gg=ggplot(m10ch,aes(ind,values,color=treatment))+
    geom_boxplot()+
  #  geom_point(position=position_jitterdodge(),shape=1)+
 # geom_boxplot(outlier.shape = NA)+
  xlab("symbiont dominance")+
  ylab("expression")+
  ggtitle("module 10")
pdf("module10_resampled_30_10_cleaned2_scaled_100.pdf",width=3.5,height=3)
gg
dev.off()

#  jackknifing support

sum(m10.h$C>m10.h$D)
# 89
sum(m10.h[,"C+D"]<m10.h$D)
# 85
sum(m10.h[,"C+D"]<m10.h$C)
# 100

sum(m10.c$C<m10.c$D)
# 86
sum(m10.c[,"C+D"]<m10.c$D)
# 84
sum(m10.c[,"C+D"]>m10.c$C)
# 62

#---- module 12 plotting and jackknifing support

m12.c=data.frame(m12.c)
m12.h=data.frame(m12.h)
m12.c$CDCD=m12.c[,1]+m12.c$CDCD
m12.c$CDD=m12.c[,1]+m12.c$CDD
names(m12.c)=c("C","C+D","D")
m12.h$CDCD=m12.h[,1]+m12.h$CDCD
m12.h$CDD=m12.h[,1]+m12.h$CDD
names(m12.h)=c("C","C+D","D")
m12h=stack(m12.h)
m12h$treatment="heat"
m12c=stack(m12.c)
m12c$treatment="control"
m12ch=data.frame(rbind(m12c,m12h))
m12ch$treatment=factor(m12ch$treatment,levels=c("heat","control"))
gg=ggplot(m12ch,aes(ind,values,color=treatment))+
  #  geom_boxplot()+
  #  geom_point(position=position_jitterdodge(),shape=1)+
  geom_boxplot(outlier.shape = NA)+
  xlab("symbiont dominance")+
  ylab("expression")+
  ggtitle("module 12")
pdf("module12_resampled_30_10_cleaned2_scaled_100.pdf",width=3.5,height=3)
gg
dev.off()

sum(m12.h$C>m12.h$D)
# 94
sum(m12.h[,"C+D"]<m12.h$D)
# 66
sum(m12.h[,"C+D"]<m12.h$C)
# 100

sum(m12.c$C<m12.c$D)
# 58
sum(m12.c[,"C+D"]<m12.c$D)
# 80
sum(m12.c[,"C+D"]<m12.c$C)
# 81
