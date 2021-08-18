#Durusdinium gene expression

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
setwd("~/Dropbox/mixed_symbiosis/")
rm(list=ls())

#ll=load('counts_metadata/coldata_redME_evenness.Rdata') 
ll=load('coldata.Rdata') 
dim(coldata)

# calculating proportions
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


ll=load('cladeD_counts')
rownames(coldata)=coldata$Run
keeps = row.names(coldata)
keeps=intersect(keeps,names(counts))
table(keeps %in% names(counts))
counts = counts[,keeps]
coldata=coldata[keeps,]
dim(coldata)

# normalizing pool factors
coldata$pool[coldata$pool=="300"]="HV"
coldata$pool[coldata$pool=="400"]="MV"

#choose read count to subsample
zooxsums=apply(counts,2,sum)
plot(coldata$evenness~log(zooxsums,10))
mincount=25000
abline(v=log(mincount,10))
table(zooxsums<mincount)
# this is in v3: removing samples with counts lower than 25k
goods=which(zooxsums>=mincount)
counts=counts[,goods]
coldata=coldata[goods,]
targetcount = rep(mincount,ncol(counts))

resampled=list()
nrep=10
for (i in 1:nrep) {
  message(i)
  res=c()
  for (s in 1:ncol(counts)) {
    probs= counts[,s]/zooxsums[s]
    cts=hist(sample(c(1:nrow(counts)), targetcount[s],prob=probs,replace=TRUE),breaks=c(0:nrow(counts)),plot=F)$counts
    res=data.frame(cbind(res, cts))
  }
  row.names(res)=row.names(counts)
  resampled[[i]]=res
}
str(resampled)
coldata.re=coldata
coldata.re$logDcount=log(apply(resampled[[1]],2,sum),10)
plot(coldata.re$evenness~coldata.re$logDcount)

means09=apply(resampled[[1]][,coldata.re$Dprop>0.9],1,mean)
means01=apply(resampled[[1]][,coldata.re$Dprop<0.1],1,mean)
table(coldata.re$Dprop<0.1)
table(coldata.re$Dprop>0.9)

himean=0.5
lomean=0.25

table(names(means01[means01>lomean]) %in% names(means09[means09>himean]))
table(names(means01[means09>himean]) %in% names(means09[means01>lomean]))
goodmeans=union(names(means09[means09>himean]),names(means01[means01>lomean]))
ngenes=length(goodmeans)
ngenes

# averaging vsds among reps, collecting pvalues and signs
vsd=c();lps=c();signs=c();ups=c();downs=c();lfc=c()
for (i in 1:nrep) {
  message(i)
  counts.re=resampled[[i]][goodmeans,]
  dds.re = DESeqDataSetFromMatrix(countData=counts.re, colData=coldata.re, design=~Dprop+pool+time)
  dds.re=DESeq(dds.re)
  vsd.re=assay(varianceStabilizingTransformation(dds.re))
  vsd.re=removeBatchEffect(vsd.re,batch=coldata.re$time,batch2=coldata.re$pool)
  res=results(dds.re,name="Dprop")
  if(i==1) { 
    lfc=res$log2FoldChange/nrep
    vsd=vsd.re/nrep 
  } else { 
    vsd=vsd+vsd.re/nrep 
    lfc=lfc+res$log2FoldChange/nrep
  }
 ups=append(ups,sum(!(is.na(res$padj)) & res$padj<=0.1 & res$log2FoldChange>0)/ngenes)
  downs=append(downs,sum(!(is.na(res$padj)) & res$padj<=0.1 & res$log2FoldChange<0)/ngenes)
    lps=data.frame(cbind(lps,res$pvalue))
  signs=data.frame(cbind(signs,res$log2FoldChange>0))
}

updown.d=stack(data.frame(cbind(ups,downs)))
boxplot(values~ind,updown.d)
plot(density(ups),xlim=c(0,1),col="red",bty="n",yaxt="n",ylab="",ylim=c(0,200),xlab="% DEGs")
lines(density(downs),col="blue")
vsd.re.d=vsd.re
coldata.re.d=coldata.re
save(vsd.re.d,coldata.re.d,updown.d,file="resample25k_D_v3b.RData")

#  writing data for GO
lfc=data.frame(cbind(gene=row.names(vsd),lfc=lfc))
write.csv(data.frame(lfc),file="~/Dropbox/mixed_symbiosis/resampledD_go_Dprop_v3b_lfc.csv",row.names=F,quote=F)

# heatmap of the averaged vsd
library(pheatmap)
library(RColorBrewer)
heat.color=colorRampPalette(c("navyblue","navyblue","navyblue","navyblue",rev(brewer.pal(n = 7, name ="RdYlBu")),"firebrick","firebrick","firebrick","firebrick"))(100)
pdf("vsd_D_resampled_heatmap_k100_v3.pdf",width=5,height=3)
pheatmap(vsd.re[,order(coldata.re$Dprop)],cluster_cols=F,kmeans_k=250,scale="row",cex=0,color=heat.color)
dev.off()

# evenness and D proportion profile
pdf("evenness_propD_v3_d.pdf",width=10,height=6)
par(mfrow=c(2,1))
plot(coldata.re$evenness[order(coldata.re$Dprop)],type="l",xlab="",ylab="evenness",bty="n",xaxt="n",mgp=c(2.3,1,0))
plot(coldata.re$Dprop[order(coldata.re$Dprop)],type="l",xlab="",ylab="proportion of D",bty="n",xaxt="n",mgp=c(2.3,1,0))
dev.off()

