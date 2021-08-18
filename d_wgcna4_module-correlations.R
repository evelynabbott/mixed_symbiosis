#wgcna4_module-correlations.R
#This script calculates and plots relationships between 
#module eigengenes and sample traits. Inputs come from:
#get_variance_stabilized_counts.R and wgcna3b_step-wise_network_construction.R
#code is adapted from examples given here: https://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/

rm(list=ls())

#setwd("~/Dropbox/project/bleww_final/scratchwork")
setwd("~/Dropbox/mixed_symbiosis/")
#=====================================================================================
#
#  Code chunk 1
#
#=====================================================================================


# Load the WGCNA package
library(WGCNA)
library(tidyverse)
source('WGCNA/wgcna_functions.R')
# The following setting is important, do not omit.
options(stringsAsFactors = FALSE);


# Load the expression and trait data saved in the first part
#lnames = load(file = "/home/evelyn/Dropbox/project/bleww_final/Zoox_expr/fightzone/C/wgcna/test2/anothertry/tenremC_varCut0.75_expCut0.75.Rdata");
#lnames = load(file = "/home/evelyn/Dropbox/project/bleww_final/Zoox_expr/fightzone/D/wgcna/test/anothertry2/tenremD_varCut0.75_expCut0.75.Rdata");
lnames = load(file = "resample25k_D_v3b.RData")

#The variable lnames contains the names of loaded variables.
lnames #(when you load a .Rdata object, if can have multiple variables assigned, the list of them is saved in this lnames variable, in this case the two objects you loaded are "datExpr" and "datTraits")


# Load network data saved in the second part.
lnames = load(file = "WGCNA/D_resampled_wgcna.RData")
#lnames = load(file = "/home/evelyn/Dropbox/project/bleww_final/Zoox_expr/fightzone/D/wgcna/test/anothertry2/wgcna3b_manual_sft12_minModSize30_cutDpth0.4_signed.Rdata") #D outliers removed
lnames


datExpr=t(vsd.re.d)
datTraits=coldata.re.d
# check if they match
table(row.names(datExpr)==row.names(datTraits))


#=====================================================================================
#
#  Code chunk 2
#
#=====================================================================================

# Define numbers of genes and samples
nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
# Recalculate module eigengenes with color labels
#This will provide the first principal component for expression
#behavior for the genes in each module. On that principal component,
#each sample will have a loading value. These values can be corelated
#with our known sample traits to get an idea what biological mechanism
#the co-reguated genes in a given module might be responding to
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs1 = MEs0[rownames(datTraits),]


########################################
############ SUBSET BY SIZE ############
########################################

#plot all modules
minModSizeToPlot = 0;outputMM=TRUE

# plot only large modules
#minModSizeToPlot = 600;outputMM=FALSE

########################################
########################################
########################################


module.sizes = table(moduleColors)
passing=names(module.sizes)[module.sizes>minModSizeToPlot]
MEs = orderMEs(MEs1[,paste('ME', passing, sep='')])

# 
########################################
########################################
#look at logCcount, CD ratio, minorLR
datTraits = datTraits %>%
  select(c_t,Dprop, evenness)
datTraits$t = ifelse(datTraits$c_t== "h",1,0)
 datTraits = datTraits %>%
   select(Dprop, evenness,t)

#use the cor() function to get the correlations between the module eigengenes and the trait data
moduleTraitCor = cor(MEs, datTraits, use = "p");
#get p values as well
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#write out module loadings
mEigs=MEs
rownames(mEigs) = rownames(datTraits)
#save(mEigs, file='moduleEigengenes_redo.Rdata')



#=====================================================================================
#
#  Code chunk 3
#
#=====================================================================================

# sizeGrWindow(10,6)
# 
# #now replot the heatmap
# sizeGrWindow(10,6)
# 
# # Will display correlations and their p-values
# textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
#                     signif(moduleTraitPvalue, 1), ")", sep = "");
# par(mar = c(6, 8.5, 3, 3));
# 
# # Display the correlation values within a heatmap plot
# COLORS = greenWhiteRed(50)
# blueWhiteRed(50)
# 
# 
# # Will display correlations and their p-values
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# 
# # Display the correlation values within a heatmap plot
# COLORS = greenWhiteRed(50)
# blueWhiteRed(50)
# blueWhiteRed(50)
#set up additional formatting variables
rows = rownames(moduleTraitCor)
sub.colors = substr(rownames(moduleTraitCor), 3, 50) #trim away the ME from the
module.sizes = paste(sub.colors, table(moduleColors)[sub.colors], sep = "\n")
dev.off()
#coul <- coul <- rev(colorRampPalette(brewer.pal(8, "RdYlBu"))(25))
pdf("D_module_trait_heatmap_resampled_3b.pdf",height=4,width=3)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = module.sizes,
               ySymbols = module.sizes,
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               cex.lab.y = .8,
               cex.lab.x = .8,
               zlim = c(-1,1),
               font.lab.y = 1,
               main = paste("Module-trait relationships"))
dev.off()
#====================================================================================
#make a heatmap showing sample vs module correlation

library(pheatmap)
library(RColorBrewer)
#Build heatmap of MEs - samples as rows and MEs as columns
MEsm = as.matrix(MEs)
#make rownames in order of most C to least C
datTraits1 <- datTraits[order(datTraits$Dprop),]
MEsm = MEsm[match(rownames(datTraits1), rownames(MEsm)),]
table(moduleColors)
coul <- coul <- rev(colorRampPalette(brewer.pal(8, "RdYlBu"))(25))
pdf("evelyn_plot_D_resampled_v3b.pdf",height=1,width=5)
pheatmap(t(MEsm),scale="row",cluster_cols=FALSE, labels_col = "",color = blueWhiteRed(50),border_color=NA)
dev.off()
pdf("evelyn_plot_D_resampled_v3b_Dprop.pdf",height=3,width=5)
plot(datTraits$Dprop[order(datTraits$Dprop)],xaxt="na",bty="n",xlab="",mgp=c(2.1,1,0),ylab="proportion of\nDurusdinium",type="l")
dev.off()


#=====================================================================================
#REPLOT THE MODULE EIGENGENE CLUSTERING TREE
#This shows you how the modules are similar to one another.
#Pushing the merging theshold (argument 3 for wgcna3b_step-wise_network_construction.R) will join modules together
#This is like moving a horizontal line up this tree figure, if the line is above a node, modules below that node will be joined into one
MEDiss = 1-cor(MEs);
# Cluster module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
sizeGrWindow(7, 6)

#plot them
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")
abline(h=0.4, lty=2, col='red')

#=====================================================================================

#Plot line sample x MEs by color
#greenyellow, blue, logCDratio
#plot(datTraits[rownames(datTraits),"logDCratio"],type="l")

MEs1 = MEs
MEs1$sample = rownames(MEs1)
#MEs1 = MEs1[match(rownames(datTraits),rownames(MEs1)),]
# 
# #Turn your 'treatment' column into a character vector
MEs1$sample <- as.character(MEs1$sample)
# #Then turn it back into a factor with the levels in the correct order
MEs1$sample <- factor(MEs1$sample, levels=unique(MEs1$sample)) 
datTraits$sample = rownames(datTraits)
datTraits$sample = as.factor(datTraits$sample)


# t = ggplot(data=MEs1, aes(x=sample, y=MEyellow, group=1)) +
#   geom_line(color="gold2")+
#   geom_point(color="gold2")+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())
# t
# b = ggplot(data=MEs1, aes(x=sample, y=MEgrey60, group=1)) +
#   geom_line(color="grey69")+
#   geom_point(color="grey69")+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())  
# b
# p = ggplot(data=MEs1, aes(x=sample, y=MEpurple, group=1)) +
#   geom_line(color="purple1")+
#   geom_point(color="purple1")+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())  
# 
# r = ggplot(data=MEs1, aes(x=sample, y=MEroyalblue, group=1)) +
#   geom_line(color="steelblue3")+
#   geom_point(color="steelblue3")+
#   theme(axis.title.x=element_blank(),
#         axis.text.x=element_blank(),
#         axis.ticks.x=element_blank())  

datTraits1 = datTraits
datTraits1$sample = rownames(datTraits1)
datTraits1$sample <- as.character(datTraits1$sample)
datTraits1$sample <- factor(datTraits1$sample, levels=unique(datTraits1$sample))  

ggplot(data=datTraits1, aes(x=sample, y=Dprop, group=1))+
  geom_point(color="black")+
  #geom_point()+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())  

ggplot(data=datTraits1, aes(x=sample,group=1))+
  geom_jitter(aes(y=evenness),color="red")+
  geom_jitter(aes(y=Dprop),color = "black")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

logCDratio = datTraits$logCDratio
MEs1 = cbind(MEs1,logCDratio)
ggplot(MEs1, aes(x=sample, group=1)) +
  geom_smooth(aes(y = MEmagenta), color = "magenta") +
  geom_smooth(aes(y = MEorange), color="orange")+
  #geom_smooth(aes(y = MElightgreen), color = "seagreen3") +
  #geom_smooth(aes(y = MEroyalblue), color = "steelblue3") +
  #geom_smooth(aes(y = MEpurple), color="purple1")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())

ggplot(MEs1, aes(x=sample, group=1)) + 
  # geom_line(aes(y = MEgrey60), color = "grey60")+
  #geom_line(aes(y = MEyellow), color="sienna4")+
  geom_smooth(aes(y = MEyellow), color = "gold2") +
  #geom_smooth(aes(y = MEgreen), color="springgreen4")+
  #geom_smooth(aes(y = MEcyan), color="turquoise3")+
  #geom_smooth(aes(y = MEdarkorange), color="orangered3")+
  geom_smooth(aes(y = MEsaddlebrown), color="grey80")+
  geom_smooth(aes(y = MEmagenta), color = "deeppink3")+
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())  





#=====================================================================================
#
#  Code chunk 4 Gather module membership data for the genes in each module
#
#=====================================================================================
#"module membership" is a measure of how strongly individual genes correlate with the 
#module eigengene. The gene that matches best with the module eigengene can be thought
#of as a hub gene, (ie it's variation in expression across the samples is most exemplary 
#of the module)

#SELECT A TRAIT AND P-VALUE CUTOFF
TRAIT='minorLR'
PCUT=1

#SUBSET FOR MODULES SIGNIFICANT FOR THAT TRAIT
traitPvalues = moduleTraitPvalue[,TRAIT]
keep = traitPvalues < PCUT
subCor = moduleTraitCor[keep,]
subP = moduleTraitPvalue[keep,]
rows = rownames(moduleTraitCor)[keep]
sub.colors = substr(rownames(subCor), 3, 50) #trim away the ME from the
module.sizes = paste(sub.colors, table(moduleColors)[sub.colors], sep = "\n")
length(sub.colors)



# # Define dataframe trait.df containing the a trait of interest from datTraits3
trait.df = as.data.frame(datTraits[,TRAIT], row.names=rownames(datTraits));
names(trait.df) = TRAIT
# names (colors) of the modules
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr[rownames(MEs),], MEs, use = "p"));
head(geneModuleMembership)
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr[rownames(trait.df),], trait.df, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(trait.df), sep="");
names(GSPvalue) = paste("p.GS.", names(trait.df), sep="");
modules = sub.colors
plot.cols = modules 

#save the module assignment data for other plotting and GO_MWU
if (outputMM){
  save(geneModuleMembership, moduleColors, file='moduleAssignment.Rdata')
}

mediumpurple3 = data.frame(gene=rownames(geneModuleMembership),
                           mediumpurple3Membership = geneModuleMembership[,'MMmediumpurple3'])

mediumpurple3$mediumpurple3Membership[!(moduleColors == "mediumpurple3")] = 0

table(mediumpurple3$mediumpurple3Membership > 0)


head(mediumpurple3)
mediumpurple3 %>% 
  write_csv(path='~/Dropbox/mixed_symbiosis/sft12_merge4_mediumpurple3Membership_ForMWU.csv')

#stopped here

ll=load('moduleAssignment.Rdata')
ll
head(geneModuleMembership)


#=====================================================================================
#
#  Code chunk 5 plot scatterplots of module membership and trait correlations
#
#=====================================================================================   
#The code below loops through each of the modules that were significant for the TRAIT of 
#interest that you chose above. Each point in the scatterpot is a gene in the module. 
#For each one of the it plots the the genes' correlation
#with the trait of interest against the genes' module memberships,
#(the correlation between the gene's variation accross samples and the module eigengene)
#A tight correlation here is suggestive that the correlated variation in gene expression 
#captured by the module is truely associated with the trait of interest.
#A tight correlation basically means, the better a gene fits into this module, the more strongly
#it correlates with the trait.

# quartz()
length(modules)
par(mfrow=c(2,2))
ggplotList = list()
for (m in modules){
  column = match(m, modNames);
  moduleGenes = moduleColors==m;
  
  # sizeGrWindow(7, 7);
  #par(mfrow = c(1,2));
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste("Module Membership in", m, "module"),
                     ylab = paste("Correlation with", TRAIT),
                     main = paste("Module membership vs. gene significance\n"),
                     cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = 'black', bg = m, pch = 21, cex = 1.5)
  ggplotList[[m]]=ggVerboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                                       abs(geneTraitSignificance[moduleGenes, 1]),
                                       xlab = paste("membership in", m, "module"),
                                       ylab = paste("correlation with", TRAIT),
                                       main = paste(""),
                                       cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = 'black', bg = m, pch = 21, cex = 1.5)
}
#ggplotList[['red']]
#plot_grid(plotlist = ggplotList, nrow=1)

