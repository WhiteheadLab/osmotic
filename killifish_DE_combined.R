library(DESeq2)
library(RColorBrewer)
library(gplots)
source('~/Documents/scripts/plotPCAWithSampleNames.R')
source('~/Documents/scripts/overLapper_original.R')
setwd("~/Documents/UCDavis/Whitehead/osmotic/osmotic")

data<-read.csv("killifish_allcounts.csv")
colnames(data)
id<-data$GeneID
rownames(data)<-id
data.1<-data[,c(2:129)]

combined_BW_FW<-c()
combined_transfer_FW<-c()
combined_transfer_BW<-c()

col.names<-colnames(data.1)
conditions = sapply(strsplit(col.names,"_"),`[`,4)
genus = sapply(strsplit(col.names,"_"),`[`,1)
species = sapply(strsplit(col.names,"_"),`[`,2)
genus_species = paste(genus,species,sep="_")
pop = sapply(strsplit(col.names,"_"),`[`,3)
genus_species_pop = paste(genus_species,pop,sep=".")
genus_species = gsub(".NA", "", genus_species_pop)
ExpDesign <- data.frame(row.names=colnames(data.1), condition = conditions,genus_species = genus_species)
ExpDesign
cds<-DESeqDataSetFromMatrix(countData=data.1, 
                            colData=ExpDesign,design= ~ condition+genus_species)
cds<-DESeq(cds, betaPrior=FALSE)
log_cds<-rlog(cds)
plotPCAWithSampleNames(log_cds,intgroup="condition",ntop=40000)
res.1<-results(cds,contrast=c("condition","BW","FW"))
dim(res.1)
res.2<-results(cds,contrast=c("condition","transfer","FW"))
res.3<-results(cds,contrast=c("condition","transfer","BW"))
resultsNames(cds)
res1_ordered <-as.data.frame(res.1[order(res.1$padj),])
dim(res1_ordered)
res1_filtered <-subset(res1_ordered,res1_ordered$padj<0.05)
res1_filtered <-subset(res1_filtered,res1_filtered$log2FoldChange>1 | res1_filtered$log2FoldChange< -1)
id<-rownames(res1_filtered)
res1_filtered<-cbind(res1_filtered,id)
dim(res1_filtered)
res2_ordered <-as.data.frame(res.2[order(res.2$padj),])
res2_filtered<-subset(res2_ordered,res2_ordered$padj<0.05)
res2_filtered <-subset(res2_filtered,res2_filtered$log2FoldChange>1 | res2_filtered$log2FoldChange< -1)
id<-rownames(res2_filtered)
res2_filtered<-cbind(res2_filtered,id)
res3_ordered<-as.data.frame(res.3[order(res.3$padj),])
res3_filtered<-subset(res3_ordered,res3_ordered$padj<0.05)
res3_filtered <-subset(res3_filtered,res3_filtered$log2FoldChange>1 | res3_filtered$log2FoldChange< -1)
id<-rownames(res3_filtered)
res3_filtered<-cbind(res3_filtered,id)
F_heteroclitus.MDPL_norm_counts<-counts(cds,normalized=TRUE)
plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (BW vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (transfer vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.3$baseMean), res.3$log2FoldChange, 
     col=ifelse(res.3$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (transfer vs. BW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
m<-res1_filtered$id
length(m)
n<-res2_filtered$id
length(n)
o<-res3_filtered$id
length(o)
setlist <- list(BW_FW=as.vector(m),transfer_FW=as.vector(n),transfer_BW=as.vector(o))
OLlist <- overLapper(setlist=setlist, sep="", type="vennsets")
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts)
# extract intersections:
names(OLlist$Venn_List)
overlap_BW_FWtransfer_FW<-OLlist$Venn_List$BW_FWtransfer_FW
length(overlap_BW_FWtransfer_FW)
overlap_BW_FWtransfer_BW<-OLlist$Venn_List$BW_FWtransfer_BW
length(overlap_BW_FWtransfer_BW)
overlap_transfer_FWtransfer_BW<-OLlist$Venn_List$transfer_FWtransfer_BW
length(overlap_transfer_FWtransfer_BW)
overlap_BW_FW<-OLlist$Venn_List$BW_FW
length(overlap_BW_FW)
overlap_transfer_FW<-OLlist$Venn_List$transfer_FW
length(overlap_transfer_FW)
overlap_transfer_BW<-OLlist$Venn_List$transfer_BW
length(overlap_transfer_BW)
combined_BW_FW<-union(overlap_BW_FW,combined_BW_FW)
length(combined_BW_FW)
combined_transfer_FW<-union(overlap_transfer_FW,combined_transfer_FW)
length(combined_transfer_FW)
combined_transfer_BW<-union(overlap_transfer_BW,combined_transfer_BW)
length(combined_transfer_BW)

#######
##F_heteroclitus.MDPP
#######
F_heteroclitus.MDPP<-data.1[,c(50:58)]
colnames(F_heteroclitus.MDPP)
col.names<-colnames(F_heteroclitus.MDPP)
head(F_heteroclitus.MDPP)
conditions = sapply(strsplit(col.names,"_"),`[`,4)
genus = sapply(strsplit(col.names,"_"),`[`,1)
species = sapply(strsplit(col.names,"_"),`[`,2)
genus_species = paste(genus,species,sep="_")
pop = sapply(strsplit(col.names,"_"),`[`,3)
genus_species_pop = paste(genus_species,pop,sep=".")
genus_species = gsub(".NA", "", genus_species_pop)
ExpDesign <- data.frame(row.names=colnames(F_heteroclitus.MDPP), condition = conditions,genus_species = genus_species)
ExpDesign
cds<-DESeqDataSetFromMatrix(countData=F_heteroclitus.MDPP, 
                            colData=ExpDesign,design= ~ condition)
cds<-DESeq(cds, betaPrior=FALSE)
log_cds<-rlog(cds)
plotPCAWithSampleNames(log_cds,intgroup="condition",ntop=40000)
res.1<-results(cds,contrast=c("condition","BW","FW"))
dim(res.1)
res.2<-results(cds,contrast=c("condition","transfer","FW"))
res.3<-results(cds,contrast=c("condition","transfer","BW"))
resultsNames(cds)
res1_ordered <-as.data.frame(res.1[order(res.1$padj),])
dim(res1_ordered)
res1_filtered <-subset(res1_ordered,res1_ordered$padj<0.05)
res1_filtered <-subset(res1_filtered,res1_filtered$log2FoldChange>1 | res1_filtered$log2FoldChange< -1)
id<-rownames(res1_filtered)
res1_filtered<-cbind(res1_filtered,id)
dim(res1_filtered)
res2_ordered <-as.data.frame(res.2[order(res.2$padj),])
res2_filtered<-subset(res2_ordered,res2_ordered$padj<0.05)
res2_filtered <-subset(res2_filtered,res2_filtered$log2FoldChange>1 | res2_filtered$log2FoldChange< -1)
id<-rownames(res2_filtered)
res2_filtered<-cbind(res2_filtered,id)
res3_ordered<-as.data.frame(res.3[order(res.3$padj),])
res3_filtered<-subset(res3_ordered,res3_ordered$padj<0.05)
res3_filtered <-subset(res3_filtered,res3_filtered$log2FoldChange>1 | res3_filtered$log2FoldChange< -1)
id<-rownames(res3_filtered)
res3_filtered<-cbind(res3_filtered,id)
F_heteroclitus.MDPP_norm_counts<-counts(cds,normalized=TRUE)
plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (BW vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (transfer vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.3$baseMean), res.3$log2FoldChange, 
     col=ifelse(res.3$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (transfer vs. BW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
m<-res1_filtered$id
length(m)
n<-res2_filtered$id
length(n)
o<-res3_filtered$id
length(o)
setlist <- list(BW_FW=as.vector(m),transfer_FW=as.vector(n),transfer_BW=as.vector(o))
OLlist <- overLapper(setlist=setlist, sep="", type="vennsets")
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts)
# extract intersections:
names(OLlist$Venn_List)
overlap_BW_FWtransfer_FW<-OLlist$Venn_List$BW_FWtransfer_FW
length(overlap_BW_FWtransfer_FW)
overlap_BW_FWtransfer_BW<-OLlist$Venn_List$BW_FWtransfer_BW
length(overlap_BW_FWtransfer_BW)
overlap_transfer_FWtransfer_BW<-OLlist$Venn_List$transfer_FWtransfer_BW
length(overlap_transfer_FWtransfer_BW)
overlap_BW_FW<-OLlist$Venn_List$BW_FW
length(overlap_BW_FW)
overlap_transfer_FW<-OLlist$Venn_List$transfer_FW
length(overlap_transfer_FW)
overlap_transfer_BW<-OLlist$Venn_List$transfer_BW
length(overlap_transfer_BW)
combined_BW_FW<-union(overlap_BW_FW,combined_BW_FW)
length(combined_BW_FW)
combined_transfer_FW<-union(overlap_transfer_FW,combined_transfer_FW)
length(combined_transfer_FW)
combined_transfer_BW<-union(overlap_transfer_BW,combined_transfer_BW)
length(combined_transfer_BW)

#######
##F_chrysotus
#######
F_chrysotus<-data.1[,c(18:25)]
colnames(F_chrysotus)
col.names<-colnames(F_chrysotus)
head(F_chrysotus)
conditions = sapply(strsplit(col.names,"_"),`[`,4)
genus = sapply(strsplit(col.names,"_"),`[`,1)
species = sapply(strsplit(col.names,"_"),`[`,2)
genus_species = paste(genus,species,sep="_")
pop = sapply(strsplit(col.names,"_"),`[`,3)
genus_species_pop = paste(genus_species,pop,sep=".")
genus_species = gsub(".NA", "", genus_species_pop)
ExpDesign <- data.frame(row.names=colnames(F_chrysotus), condition = conditions,genus_species = genus_species)
ExpDesign
cds<-DESeqDataSetFromMatrix(countData=F_chrysotus, colData=ExpDesign,design= ~ condition)
cds<-DESeq(cds, betaPrior=FALSE)
log_cds<-rlog(cds)
plotPCAWithSampleNames(log_cds,intgroup="condition",ntop=40000)
res.1<-results(cds,contrast=c("condition","BW","FW"))
dim(res.1)
res.2<-results(cds,contrast=c("condition","transfer","FW"))
res.3<-results(cds,contrast=c("condition","transfer","BW"))
resultsNames(cds)
res1_ordered <-as.data.frame(res.1[order(res.1$padj),])
dim(res1_ordered)
res1_filtered <-subset(res1_ordered,res1_ordered$padj<0.05)
res1_filtered <-subset(res1_filtered,res1_filtered$log2FoldChange>1 | res1_filtered$log2FoldChange< -1)
id<-rownames(res1_filtered)
res1_filtered<-cbind(res1_filtered,id)
dim(res1_filtered)
res2_ordered <-as.data.frame(res.2[order(res.2$padj),])
res2_filtered<-subset(res2_ordered,res2_ordered$padj<0.05)
res2_filtered <-subset(res2_filtered,res2_filtered$log2FoldChange>1 | res2_filtered$log2FoldChange< -1)
id<-rownames(res2_filtered)
res2_filtered<-cbind(res2_filtered,id)
res3_ordered<-as.data.frame(res.3[order(res.3$padj),])
res3_filtered<-subset(res3_ordered,res3_ordered$padj<0.05)
res3_filtered <-subset(res3_filtered,res3_filtered$log2FoldChange>1 | res3_filtered$log2FoldChange< -1)
id<-rownames(res3_filtered)
res3_filtered<-cbind(res3_filtered,id)
F_chrysotus_norm_counts<-counts(cds,normalized=TRUE)
plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (BW vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (transfer vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.3$baseMean), res.3$log2FoldChange, 
     col=ifelse(res.3$padj < 0.05, "red","gray67"),
     main="F. chrysotus (transfer vs. BW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
m<-res1_filtered$id
length(m)
n<-res2_filtered$id
length(n)
o<-res3_filtered$id
length(o)
setlist <- list(BW_FW=as.vector(m),transfer_FW=as.vector(n),transfer_BW=as.vector(o))
OLlist <- overLapper(setlist=setlist, sep="", type="vennsets")
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts)
# extract intersections:
names(OLlist$Venn_List)
overlap_BW_FWtransfer_FW<-OLlist$Venn_List$BW_FWtransfer_FW
length(overlap_BW_FWtransfer_FW)
overlap_BW_FWtransfer_BW<-OLlist$Venn_List$BW_FWtransfer_BW
length(overlap_BW_FWtransfer_BW)
overlap_transfer_FWtransfer_BW<-OLlist$Venn_List$transfer_FWtransfer_BW
length(overlap_transfer_FWtransfer_BW)
overlap_BW_FW<-OLlist$Venn_List$BW_FW
length(overlap_BW_FW)
overlap_transfer_FW<-OLlist$Venn_List$transfer_FW
length(overlap_transfer_FW)
overlap_transfer_BW<-OLlist$Venn_List$transfer_BW
length(overlap_transfer_BW)
combined_BW_FW<-union(overlap_BW_FW,combined_BW_FW)
length(combined_BW_FW)
combined_transfer_FW<-union(overlap_transfer_FW,combined_transfer_FW)
length(combined_transfer_FW)
combined_transfer_BW<-union(overlap_transfer_BW,combined_transfer_BW)
length(combined_transfer_BW)

#######
##F_diaphanus
#######

F_diaphanus<-data.1[,c(26:31)]
colnames(F_diaphanus)
col.names<-colnames(F_diaphanus)
head(F_diaphanus)
conditions = sapply(strsplit(col.names,"_"),`[`,4)
genus = sapply(strsplit(col.names,"_"),`[`,1)
species = sapply(strsplit(col.names,"_"),`[`,2)
genus_species = paste(genus,species,sep="_")
pop = sapply(strsplit(col.names,"_"),`[`,3)
genus_species_pop = paste(genus_species,pop,sep=".")
genus_species = gsub(".NA", "", genus_species_pop)
ExpDesign <- data.frame(row.names=colnames(F_diaphanus), condition = conditions,genus_species = genus_species)
ExpDesign
cds<-DESeqDataSetFromMatrix(countData=F_diaphanus, colData=ExpDesign,design= ~ condition)
cds<-DESeq(cds, betaPrior=FALSE)
log_cds<-rlog(cds)
plotPCAWithSampleNames(log_cds,intgroup="condition",ntop=40000)
res.1<-results(cds,contrast=c("condition","BW","FW"))
dim(res.1)
res.2<-results(cds,contrast=c("condition","transfer","FW"))
res.3<-results(cds,contrast=c("condition","transfer","BW"))
resultsNames(cds)
res1_ordered <-as.data.frame(res.1[order(res.1$padj),])
dim(res1_ordered)
res1_filtered <-subset(res1_ordered,res1_ordered$padj<0.05)
res1_filtered <-subset(res1_filtered,res1_filtered$log2FoldChange>1 | res1_filtered$log2FoldChange< -1)
id<-rownames(res1_filtered)
res1_filtered<-cbind(res1_filtered,id)
dim(res1_filtered)
res2_ordered <-as.data.frame(res.2[order(res.2$padj),])
res2_filtered<-subset(res2_ordered,res2_ordered$padj<0.05)
res2_filtered <-subset(res2_filtered,res2_filtered$log2FoldChange>1 | res2_filtered$log2FoldChange< -1)
id<-rownames(res2_filtered)
res2_filtered<-cbind(res2_filtered,id)
res3_ordered<-as.data.frame(res.3[order(res.3$padj),])
res3_filtered<-subset(res3_ordered,res3_ordered$padj<0.05)
res3_filtered <-subset(res3_filtered,res3_filtered$log2FoldChange>1 | res3_filtered$log2FoldChange< -1)
id<-rownames(res3_filtered)
res3_filtered<-cbind(res3_filtered,id)
F_diaphanus_norm_counts<-counts(cds,normalized=TRUE)
plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="F. diaphanus (BW vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (transfer vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.3$baseMean), res.3$log2FoldChange, 
     col=ifelse(res.3$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (transfer vs. BW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
m<-res1_filtered$id
length(m)
n<-res2_filtered$id
length(n)
o<-res3_filtered$id
length(o)
setlist <- list(BW_FW=as.vector(m),transfer_FW=as.vector(n),transfer_BW=as.vector(o))
OLlist <- overLapper(setlist=setlist, sep="", type="vennsets")
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts)
# extract intersections:
names(OLlist$Venn_List)
overlap_BW_FWtransfer_FW<-OLlist$Venn_List$BW_FWtransfer_FW
length(overlap_BW_FWtransfer_FW)
overlap_BW_FWtransfer_BW<-OLlist$Venn_List$BW_FWtransfer_BW
length(overlap_BW_FWtransfer_BW)
overlap_transfer_FWtransfer_BW<-OLlist$Venn_List$transfer_FWtransfer_BW
length(overlap_transfer_FWtransfer_BW)
overlap_BW_FW<-OLlist$Venn_List$BW_FW
length(overlap_BW_FW)
overlap_transfer_FW<-OLlist$Venn_List$transfer_FW
length(overlap_transfer_FW)
overlap_transfer_BW<-OLlist$Venn_List$transfer_BW
length(overlap_transfer_BW)
combined_BW_FW<-union(overlap_BW_FW,combined_BW_FW)
length(combined_BW_FW)
combined_transfer_FW<-union(overlap_transfer_FW,combined_transfer_FW)
length(combined_transfer_FW)
combined_transfer_BW<-union(overlap_transfer_BW,combined_transfer_BW)
length(combined_transfer_BW)

#######
##F_grandis
#######

F_grandis<-data.1[,c(32:40)]
colnames(F_grandis)
col.names<-colnames(F_grandis)
head(F_grandis)
conditions = sapply(strsplit(col.names,"_"),`[`,4)
genus = sapply(strsplit(col.names,"_"),`[`,1)
species = sapply(strsplit(col.names,"_"),`[`,2)
genus_species = paste(genus,species,sep="_")
pop = sapply(strsplit(col.names,"_"),`[`,3)
genus_species_pop = paste(genus_species,pop,sep=".")
genus_species = gsub(".NA", "", genus_species_pop)
ExpDesign <- data.frame(row.names=colnames(F_grandis), condition = conditions,genus_species = genus_species)
ExpDesign
cds<-DESeqDataSetFromMatrix(countData=F_grandis, colData=ExpDesign,design= ~ condition)
cds<-DESeq(cds, betaPrior=FALSE)
log_cds<-rlog(cds)
plotPCAWithSampleNames(log_cds,intgroup="condition",ntop=40000)
res.1<-results(cds,contrast=c("condition","BW","FW"))
dim(res.1)
res.2<-results(cds,contrast=c("condition","transfer","FW"))
res.3<-results(cds,contrast=c("condition","transfer","BW"))
resultsNames(cds)
res1_ordered <-as.data.frame(res.1[order(res.1$padj),])
dim(res1_ordered)
res1_filtered <-subset(res1_ordered,res1_ordered$padj<0.05)
res1_filtered <-subset(res1_filtered,res1_filtered$log2FoldChange>1 | res1_filtered$log2FoldChange< -1)
id<-rownames(res1_filtered)
res1_filtered<-cbind(res1_filtered,id)
dim(res1_filtered)
res2_ordered <-as.data.frame(res.2[order(res.2$padj),])
res2_filtered<-subset(res2_ordered,res2_ordered$padj<0.05)
res2_filtered <-subset(res2_filtered,res2_filtered$log2FoldChange>1 | res2_filtered$log2FoldChange< -1)
id<-rownames(res2_filtered)
res2_filtered<-cbind(res2_filtered,id)
res3_ordered<-as.data.frame(res.3[order(res.3$padj),])
res3_filtered<-subset(res3_ordered,res3_ordered$padj<0.05)
res3_filtered <-subset(res3_filtered,res3_filtered$log2FoldChange>1 | res3_filtered$log2FoldChange< -1)
id<-rownames(res3_filtered)
res3_filtered<-cbind(res3_filtered,id)
F_grandis_norm_counts<-counts(cds,normalized=TRUE)
plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (BW vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (transfer vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.3$baseMean), res.3$log2FoldChange, 
     col=ifelse(res.3$padj < 0.05, "red","gray67"),
     main="F. grandis (transfer vs. BW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
m<-res1_filtered$id
length(m)
n<-res2_filtered$id
length(n)
o<-res3_filtered$id
length(o)
setlist <- list(BW_FW=as.vector(m),transfer_FW=as.vector(n),transfer_BW=as.vector(o))
OLlist <- overLapper(setlist=setlist, sep="", type="vennsets")
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts)
# extract intersections:
names(OLlist$Venn_List)
overlap_BW_FWtransfer_FW<-OLlist$Venn_List$BW_FWtransfer_FW
length(overlap_BW_FWtransfer_FW)
overlap_BW_FWtransfer_BW<-OLlist$Venn_List$BW_FWtransfer_BW
length(overlap_BW_FWtransfer_BW)
overlap_transfer_FWtransfer_BW<-OLlist$Venn_List$transfer_FWtransfer_BW
length(overlap_transfer_FWtransfer_BW)
overlap_BW_FW<-OLlist$Venn_List$BW_FW
length(overlap_BW_FW)
overlap_transfer_FW<-OLlist$Venn_List$transfer_FW
length(overlap_transfer_FW)
overlap_transfer_BW<-OLlist$Venn_List$transfer_BW
length(overlap_transfer_BW)
combined_BW_FW<-union(overlap_BW_FW,combined_BW_FW)
length(combined_BW_FW)
combined_transfer_FW<-union(overlap_transfer_FW,combined_transfer_FW)
length(combined_transfer_FW)
combined_transfer_BW<-union(overlap_transfer_BW,combined_transfer_BW)
length(combined_transfer_BW)

#######
##F_olivaceous
#######

F_olivaceous<-data.1[,c(70:77)]
colnames(F_olivaceous)
col.names<-colnames(F_olivaceous)
head(F_olivaceous)
conditions = sapply(strsplit(col.names,"_"),`[`,4)
genus = sapply(strsplit(col.names,"_"),`[`,1)
species = sapply(strsplit(col.names,"_"),`[`,2)
genus_species = paste(genus,species,sep="_")
pop = sapply(strsplit(col.names,"_"),`[`,3)
genus_species_pop = paste(genus_species,pop,sep=".")
genus_species = gsub(".NA", "", genus_species_pop)
ExpDesign <- data.frame(row.names=colnames(F_olivaceous), condition = conditions,genus_species = genus_species)
ExpDesign
cds<-DESeqDataSetFromMatrix(countData=F_olivaceous, colData=ExpDesign,design= ~ condition)
cds<-DESeq(cds, betaPrior=FALSE)
log_cds<-rlog(cds)
plotPCAWithSampleNames(log_cds,intgroup="condition",ntop=40000)
res.1<-results(cds,contrast=c("condition","BW","FW"))
dim(res.1)
res.2<-results(cds,contrast=c("condition","transfer","FW"))
res.3<-results(cds,contrast=c("condition","transfer","BW"))
resultsNames(cds)
res1_ordered <-as.data.frame(res.1[order(res.1$padj),])
dim(res1_ordered)
res1_filtered <-subset(res1_ordered,res1_ordered$padj<0.05)
res1_filtered <-subset(res1_filtered,res1_filtered$log2FoldChange>1 | res1_filtered$log2FoldChange< -1)
id<-rownames(res1_filtered)
res1_filtered<-cbind(res1_filtered,id)
dim(res1_filtered)
res2_ordered <-as.data.frame(res.2[order(res.2$padj),])
res2_filtered<-subset(res2_ordered,res2_ordered$padj<0.05)
res2_filtered <-subset(res2_filtered,res2_filtered$log2FoldChange>1 | res2_filtered$log2FoldChange< -1)
id<-rownames(res2_filtered)
res2_filtered<-cbind(res2_filtered,id)
res3_ordered<-as.data.frame(res.3[order(res.3$padj),])
res3_filtered<-subset(res3_ordered,res3_ordered$padj<0.05)
res3_filtered <-subset(res3_filtered,res3_filtered$log2FoldChange>1 | res3_filtered$log2FoldChange< -1)
id<-rownames(res3_filtered)
res3_filtered<-cbind(res3_filtered,id)
F_olivaceous_norm_counts<-counts(cds,normalized=TRUE)
plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="F. olivaceous (BW vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (transfer vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.3$baseMean), res.3$log2FoldChange, 
     col=ifelse(res.3$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (transfer vs. BW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
m<-res1_filtered$id
length(m)
n<-res2_filtered$id
length(n)
o<-res3_filtered$id
length(o)
setlist <- list(BW_FW=as.vector(m),transfer_FW=as.vector(n),transfer_BW=as.vector(o))
OLlist <- overLapper(setlist=setlist, sep="", type="vennsets")
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts)
# extract intersections:
names(OLlist$Venn_List)
overlap_BW_FWtransfer_FW<-OLlist$Venn_List$BW_FWtransfer_FW
length(overlap_BW_FWtransfer_FW)
overlap_BW_FWtransfer_BW<-OLlist$Venn_List$BW_FWtransfer_BW
length(overlap_BW_FWtransfer_BW)
overlap_transfer_FWtransfer_BW<-OLlist$Venn_List$transfer_FWtransfer_BW
length(overlap_transfer_FWtransfer_BW)
overlap_BW_FW<-OLlist$Venn_List$BW_FW
length(overlap_BW_FW)
overlap_transfer_FW<-OLlist$Venn_List$transfer_FW
length(overlap_transfer_FW)
overlap_transfer_BW<-OLlist$Venn_List$transfer_BW
length(overlap_transfer_BW)
combined_BW_FW<-union(overlap_BW_FW,combined_BW_FW)
length(combined_BW_FW)
combined_transfer_FW<-union(overlap_transfer_FW,combined_transfer_FW)
length(combined_transfer_FW)
combined_transfer_BW<-union(overlap_transfer_BW,combined_transfer_BW)
length(combined_transfer_BW)

#######
##F_parvapinis
#######

F_parvapinis<-data.1[,c(78:85)]
colnames(F_parvapinis)
col.names<-colnames(F_parvapinis)
head(F_parvapinis)
conditions = sapply(strsplit(col.names,"_"),`[`,4)
genus = sapply(strsplit(col.names,"_"),`[`,1)
species = sapply(strsplit(col.names,"_"),`[`,2)
genus_species = paste(genus,species,sep="_")
pop = sapply(strsplit(col.names,"_"),`[`,3)
genus_species_pop = paste(genus_species,pop,sep=".")
genus_species = gsub(".NA", "", genus_species_pop)
ExpDesign <- data.frame(row.names=colnames(F_parvapinis), condition = conditions,genus_species = genus_species)
ExpDesign
cds<-DESeqDataSetFromMatrix(countData=F_parvapinis, colData=ExpDesign,design= ~ condition)
cds<-DESeq(cds, betaPrior=FALSE)
log_cds<-rlog(cds)
plotPCAWithSampleNames(log_cds,intgroup="condition",ntop=40000)
res.1<-results(cds,contrast=c("condition","BW","FW"))
dim(res.1)
res.2<-results(cds,contrast=c("condition","transfer","FW"))
res.3<-results(cds,contrast=c("condition","transfer","BW"))
resultsNames(cds)
res1_ordered <-as.data.frame(res.1[order(res.1$padj),])
dim(res1_ordered)
res1_filtered <-subset(res1_ordered,res1_ordered$padj<0.05)
res1_filtered <-subset(res1_filtered,res1_filtered$log2FoldChange>1 | res1_filtered$log2FoldChange< -1)
id<-rownames(res1_filtered)
res1_filtered<-cbind(res1_filtered,id)
dim(res1_filtered)
res2_ordered <-as.data.frame(res.2[order(res.2$padj),])
res2_filtered<-subset(res2_ordered,res2_ordered$padj<0.05)
res2_filtered <-subset(res2_filtered,res2_filtered$log2FoldChange>1 | res2_filtered$log2FoldChange< -1)
id<-rownames(res2_filtered)
res2_filtered<-cbind(res2_filtered,id)
res3_ordered<-as.data.frame(res.3[order(res.3$padj),])
res3_filtered<-subset(res3_ordered,res3_ordered$padj<0.05)
res3_filtered <-subset(res3_filtered,res3_filtered$log2FoldChange>1 | res3_filtered$log2FoldChange< -1)
id<-rownames(res3_filtered)
res3_filtered<-cbind(res3_filtered,id)
F_parvapinis_norm_counts<-counts(cds,normalized=TRUE)
plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="F. parvapinis (BW vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (transfer vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.3$baseMean), res.3$log2FoldChange, 
     col=ifelse(res.3$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (transfer vs. BW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
m<-res1_filtered$id
length(m)
n<-res2_filtered$id
length(n)
o<-res3_filtered$id
length(o)
setlist <- list(BW_FW=as.vector(m),transfer_FW=as.vector(n),transfer_BW=as.vector(o))
OLlist <- overLapper(setlist=setlist, sep="", type="vennsets")
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts)
# extract intersections:
names(OLlist$Venn_List)
overlap_BW_FWtransfer_FW<-OLlist$Venn_List$BW_FWtransfer_FW
length(overlap_BW_FWtransfer_FW)
overlap_BW_FWtransfer_BW<-OLlist$Venn_List$BW_FWtransfer_BW
length(overlap_BW_FWtransfer_BW)
overlap_transfer_FWtransfer_BW<-OLlist$Venn_List$transfer_FWtransfer_BW
length(overlap_transfer_FWtransfer_BW)
overlap_BW_FW<-OLlist$Venn_List$BW_FW
length(overlap_BW_FW)
overlap_transfer_FW<-OLlist$Venn_List$transfer_FW
length(overlap_transfer_FW)
overlap_transfer_BW<-OLlist$Venn_List$transfer_BW
length(overlap_transfer_BW)
combined_BW_FW<-union(overlap_BW_FW,combined_BW_FW)
length(combined_BW_FW)
combined_transfer_FW<-union(overlap_transfer_FW,combined_transfer_FW)
length(combined_transfer_FW)
combined_transfer_BW<-union(overlap_transfer_BW,combined_transfer_BW)
length(combined_transfer_BW)

#######
##F_rathbuni
#######

F_rathbuni<-data.1[,c(86:94)]
colnames(F_rathbuni)
col.names<-colnames(F_rathbuni)
head(F_rathbuni)
conditions = sapply(strsplit(col.names,"_"),`[`,4)
genus = sapply(strsplit(col.names,"_"),`[`,1)
species = sapply(strsplit(col.names,"_"),`[`,2)
genus_species = paste(genus,species,sep="_")
pop = sapply(strsplit(col.names,"_"),`[`,3)
genus_species_pop = paste(genus_species,pop,sep=".")
genus_species = gsub(".NA", "", genus_species_pop)
ExpDesign <- data.frame(row.names=colnames(F_rathbuni), condition = conditions,genus_species = genus_species)
ExpDesign
cds<-DESeqDataSetFromMatrix(countData=F_rathbuni, colData=ExpDesign,design= ~ condition)
cds<-DESeq(cds, betaPrior=FALSE)
log_cds<-rlog(cds)
plotPCAWithSampleNames(log_cds,intgroup="condition",ntop=40000)
res.1<-results(cds,contrast=c("condition","BW","FW"))
dim(res.1)
res.2<-results(cds,contrast=c("condition","transfer","FW"))
res.3<-results(cds,contrast=c("condition","transfer","BW"))
resultsNames(cds)
res1_ordered <-as.data.frame(res.1[order(res.1$padj),])
dim(res1_ordered)
res1_filtered <-subset(res1_ordered,res1_ordered$padj<0.05)
res1_filtered <-subset(res1_filtered,res1_filtered$log2FoldChange>1 | res1_filtered$log2FoldChange< -1)
id<-rownames(res1_filtered)
res1_filtered<-cbind(res1_filtered,id)
dim(res1_filtered)
res2_ordered <-as.data.frame(res.2[order(res.2$padj),])
res2_filtered<-subset(res2_ordered,res2_ordered$padj<0.05)
res2_filtered <-subset(res2_filtered,res2_filtered$log2FoldChange>1 | res2_filtered$log2FoldChange< -1)
id<-rownames(res2_filtered)
res2_filtered<-cbind(res2_filtered,id)
res3_ordered<-as.data.frame(res.3[order(res.3$padj),])
res3_filtered<-subset(res3_ordered,res3_ordered$padj<0.05)
res3_filtered <-subset(res3_filtered,res3_filtered$log2FoldChange>1 | res3_filtered$log2FoldChange< -1)
id<-rownames(res3_filtered)
res3_filtered<-cbind(res3_filtered,id)
F_rathbuni_norm_counts<-counts(cds,normalized=TRUE)
plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (BW vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="F. rathbuni (transfer vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.3$baseMean), res.3$log2FoldChange, 
     col=ifelse(res.3$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (transfer vs. BW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
m<-res1_filtered$id
length(m)
n<-res2_filtered$id
length(n)
o<-res3_filtered$id
length(o)
setlist <- list(BW_FW=as.vector(m),transfer_FW=as.vector(n),transfer_BW=as.vector(o))
OLlist <- overLapper(setlist=setlist, sep="", type="vennsets")
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts)
# extract intersections:
names(OLlist$Venn_List)
overlap_BW_FWtransfer_FW<-OLlist$Venn_List$BW_FWtransfer_FW
length(overlap_BW_FWtransfer_FW)
overlap_BW_FWtransfer_BW<-OLlist$Venn_List$BW_FWtransfer_BW
length(overlap_BW_FWtransfer_BW)
overlap_transfer_FWtransfer_BW<-OLlist$Venn_List$transfer_FWtransfer_BW
length(overlap_transfer_FWtransfer_BW)
overlap_BW_FW<-OLlist$Venn_List$BW_FW
length(overlap_BW_FW)
overlap_transfer_FW<-OLlist$Venn_List$transfer_FW
length(overlap_transfer_FW)
overlap_transfer_BW<-OLlist$Venn_List$transfer_BW
length(overlap_transfer_BW)
combined_BW_FW<-union(overlap_BW_FW,combined_BW_FW)
length(combined_BW_FW)
combined_transfer_FW<-union(overlap_transfer_FW,combined_transfer_FW)
length(combined_transfer_FW)
combined_transfer_BW<-union(overlap_transfer_BW,combined_transfer_BW)
length(combined_transfer_BW)

#######
##F_goodei
#######

L_goodei<-data.1[,c(112:120)]
colnames(L_goodei)
col.names<-colnames(L_goodei)
head(L_goodei)
conditions = sapply(strsplit(col.names,"_"),`[`,4)
genus = sapply(strsplit(col.names,"_"),`[`,1)
species = sapply(strsplit(col.names,"_"),`[`,2)
genus_species = paste(genus,species,sep="_")
pop = sapply(strsplit(col.names,"_"),`[`,3)
genus_species_pop = paste(genus_species,pop,sep=".")
genus_species = gsub(".NA", "", genus_species_pop)
ExpDesign <- data.frame(row.names=colnames(L_goodei), condition = conditions,genus_species = genus_species)
ExpDesign
cds<-DESeqDataSetFromMatrix(countData=L_goodei, colData=ExpDesign,design= ~ condition)
cds<-DESeq(cds, betaPrior=FALSE)
log_cds<-rlog(cds)
plotPCAWithSampleNames(log_cds,intgroup="condition",ntop=40000)
res.1<-results(cds,contrast=c("condition","BW","FW"))
dim(res.1)
res.2<-results(cds,contrast=c("condition","transfer","FW"))
res.3<-results(cds,contrast=c("condition","transfer","BW"))
resultsNames(cds)
res1_ordered <-as.data.frame(res.1[order(res.1$padj),])
dim(res1_ordered)
res1_filtered <-subset(res1_ordered,res1_ordered$padj<0.05)
res1_filtered <-subset(res1_filtered,res1_filtered$log2FoldChange>1 | res1_filtered$log2FoldChange< -1)
id<-rownames(res1_filtered)
res1_filtered<-cbind(res1_filtered,id)
dim(res1_filtered)
res2_ordered <-as.data.frame(res.2[order(res.2$padj),])
res2_filtered<-subset(res2_ordered,res2_ordered$padj<0.05)
res2_filtered <-subset(res2_filtered,res2_filtered$log2FoldChange>1 | res2_filtered$log2FoldChange< -1)
id<-rownames(res2_filtered)
res2_filtered<-cbind(res2_filtered,id)
res3_ordered<-as.data.frame(res.3[order(res.3$padj),])
res3_filtered<-subset(res3_ordered,res3_ordered$padj<0.05)
res3_filtered <-subset(res3_filtered,res3_filtered$log2FoldChange>1 | res3_filtered$log2FoldChange< -1)
id<-rownames(res3_filtered)
res3_filtered<-cbind(res3_filtered,id)
L_goodei_norm_counts<-counts(cds,normalized=TRUE)
plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="L. goodei (BW vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="L_goodei (transfer vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.3$baseMean), res.3$log2FoldChange, 
     col=ifelse(res.3$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPL (transfer vs. BW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
m<-res1_filtered$id
length(m)
n<-res2_filtered$id
length(n)
o<-res3_filtered$id
length(o)
setlist <- list(BW_FW=as.vector(m),transfer_FW=as.vector(n),transfer_BW=as.vector(o))
OLlist <- overLapper(setlist=setlist, sep="", type="vennsets")
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts)
# extract intersections:
names(OLlist$Venn_List)
overlap_BW_FWtransfer_FW<-OLlist$Venn_List$BW_FWtransfer_FW
length(overlap_BW_FWtransfer_FW)
overlap_BW_FWtransfer_BW<-OLlist$Venn_List$BW_FWtransfer_BW
length(overlap_BW_FWtransfer_BW)
overlap_transfer_FWtransfer_BW<-OLlist$Venn_List$transfer_FWtransfer_BW
length(overlap_transfer_FWtransfer_BW)
overlap_BW_FW<-OLlist$Venn_List$BW_FW
length(overlap_BW_FW)
overlap_transfer_FW<-OLlist$Venn_List$transfer_FW
length(overlap_transfer_FW)
overlap_transfer_BW<-OLlist$Venn_List$transfer_BW
length(overlap_transfer_BW)
combined_BW_FW<-union(overlap_BW_FW,combined_BW_FW)
length(combined_BW_FW)
combined_transfer_FW<-union(overlap_transfer_FW,combined_transfer_FW)
length(combined_transfer_FW)
combined_transfer_BW<-union(overlap_transfer_BW,combined_transfer_BW)
length(combined_transfer_BW)

#######
##L_parva
#######

L_parva<-data.1[,c(121:129)]
colnames(L_parva)
col.names<-colnames(L_parva)
head(L_parva)
conditions = sapply(strsplit(col.names,"_"),`[`,4)
genus = sapply(strsplit(col.names,"_"),`[`,1)
species = sapply(strsplit(col.names,"_"),`[`,2)
genus_species = paste(genus,species,sep="_")
pop = sapply(strsplit(col.names,"_"),`[`,3)
genus_species_pop = paste(genus_species,pop,sep=".")
genus_species = gsub(".NA", "", genus_species_pop)
ExpDesign <- data.frame(row.names=colnames(L_parva), condition = conditions,genus_species = genus_species)
ExpDesign
cds<-DESeqDataSetFromMatrix(countData=L_parva, colData=ExpDesign,design= ~ condition)
cds<-DESeq(cds, betaPrior=FALSE)
log_cds<-rlog(cds)
plotPCAWithSampleNames(log_cds,intgroup="condition",ntop=40000)
res.1<-results(cds,contrast=c("condition","BW","FW"))
dim(res.1)
res.2<-results(cds,contrast=c("condition","transfer","FW"))
res.3<-results(cds,contrast=c("condition","transfer","BW"))
resultsNames(cds)
res1_ordered <-as.data.frame(res.1[order(res.1$padj),])
dim(res1_ordered)
res1_filtered <-subset(res1_ordered,res1_ordered$padj<0.05)
res1_filtered <-subset(res1_filtered,res1_filtered$log2FoldChange>1 | res1_filtered$log2FoldChange< -1)
id<-rownames(res1_filtered)
res1_filtered<-cbind(res1_filtered,id)
dim(res1_filtered)
res2_ordered <-as.data.frame(res.2[order(res.2$padj),])
res2_filtered<-subset(res2_ordered,res2_ordered$padj<0.05)
res2_filtered <-subset(res2_filtered,res2_filtered$log2FoldChange>1 | res2_filtered$log2FoldChange< -1)
id<-rownames(res2_filtered)
res2_filtered<-cbind(res2_filtered,id)
res3_ordered<-as.data.frame(res.3[order(res.3$padj),])
res3_filtered<-subset(res3_ordered,res3_ordered$padj<0.05)
res3_filtered <-subset(res3_filtered,res3_filtered$log2FoldChange>1 | res3_filtered$log2FoldChange< -1)
id<-rownames(res3_filtered)
res3_filtered<-cbind(res3_filtered,id)
L_parva_norm_counts<-counts(cds,normalized=TRUE)
plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="L_parva (BW vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="L_goodei (transfer vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.3$baseMean), res.3$log2FoldChange, 
     col=ifelse(res.3$padj < 0.05, "red","gray67"),
     main="L. parva (transfer vs. BW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
m<-res1_filtered$id
length(m)
n<-res2_filtered$id
length(n)
o<-res3_filtered$id
length(o)
setlist <- list(BW_FW=as.vector(m),transfer_FW=as.vector(n),transfer_BW=as.vector(o))
OLlist <- overLapper(setlist=setlist, sep="", type="vennsets")
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts)
# extract intersections:
names(OLlist$Venn_List)
overlap_BW_FWtransfer_FW<-OLlist$Venn_List$BW_FWtransfer_FW
length(overlap_BW_FWtransfer_FW)
overlap_BW_FWtransfer_BW<-OLlist$Venn_List$BW_FWtransfer_BW
length(overlap_BW_FWtransfer_BW)
overlap_transfer_FWtransfer_BW<-OLlist$Venn_List$transfer_FWtransfer_BW
length(overlap_transfer_FWtransfer_BW)
overlap_BW_FW<-OLlist$Venn_List$BW_FW
length(overlap_BW_FW)
overlap_transfer_FW<-OLlist$Venn_List$transfer_FW
length(overlap_transfer_FW)
overlap_transfer_BW<-OLlist$Venn_List$transfer_BW
length(overlap_transfer_BW)
combined_BW_FW<-union(overlap_BW_FW,combined_BW_FW)
length(combined_BW_FW)
combined_transfer_FW<-union(overlap_transfer_FW,combined_transfer_FW)
length(combined_transfer_FW)
combined_transfer_BW<-union(overlap_transfer_BW,combined_transfer_BW)
length(combined_transfer_BW)


m<-combined_BW_FW
length(m)
n<-combined_transfer_FW
length(n)
o<-combined_transfer_BW
length(o)
setlist <- list(BW_FW=as.vector(m),transfer_FW=as.vector(n),transfer_BW=as.vector(o))
OLlist <- overLapper(setlist=setlist, sep="", type="vennsets")
counts <- sapply(OLlist$Venn_List, length)

png(filename="vann_all.png",width=3.25,height=3.25,units="in",res=1200,pointsize = 4)
vennPlot(counts=counts)
dev.off()

overlap_BW_FW_transfer <- OLlist$Venn_List$BW_FWtransfer_FW
overlap_BW_FW_transfer <- union(overlap_BW_FW_transfer,OLlist$Venn_List$BW_FWtransfer_BW)
overlap_BW_FW_transfer <- union(overlap_BW_FW_transfer,OLlist$Venn_List$transfer_FWtransfer_BW)
overlap_BW_FW_transfer <- union(overlap_BW_FW_transfer,OLlist$Venn_List$BW_FWtransfer_FWtransfer_BW)
length(overlap_BW_FW_transfer)


############
# Marine
############
# F_heteroclitus.MDPL
sig_counts <- as.data.frame(F_heteroclitus.MDPL_norm_counts)
row.names <- rownames(sig_counts)
sig_counts <- cbind(sig_counts,id=row.names)
dim(sig_counts)
# F_heteroclitus.MDPP
F_heteroclitus.MDPP_counts <- as.data.frame(F_heteroclitus.MDPP_norm_counts)
row.names <- rownames(F_heteroclitus.MDPP_counts)
F_heteroclitus.MDPP_counts <- cbind(F_heteroclitus.MDPP_counts,id=row.names)
sig_counts<-merge(sig_counts,F_heteroclitus.MDPP_counts,by="id")
dim(sig_counts)
# F_grandis
F_grandis_counts <- as.data.frame(F_grandis_norm_counts)
row.names <- rownames(F_grandis_norm_counts)
F_grandis_counts <- cbind(F_grandis_counts,id=row.names)
sig_counts<-merge(sig_counts,F_grandis_counts,by="id")
dim(sig_counts)
# L_parva
L_parva_counts <- as.data.frame(L_parva_norm_counts)
row.names <- rownames(L_parva_counts)
L_parva_counts <- cbind(L_parva_counts,id=row.names)
sig_counts<-merge(sig_counts,L_parva_counts,by="id")
dim(sig_counts)
############
# Brackish
############
# F_chrysotus
F_chrysotus_counts <- as.data.frame(F_chrysotus_norm_counts)
row.names <- rownames(F_chrysotus_counts)
F_chrysotus_counts <- cbind(F_chrysotus_counts,id=row.names)
sig_counts<-merge(sig_counts,F_chrysotus_counts,by="id")
dim(sig_counts)
# F_diaphanus
F_diaphanus_counts <- as.data.frame(F_diaphanus_norm_counts)
row.names <- rownames(F_diaphanus_counts)
F_diaphanus_counts <- cbind(F_diaphanus_counts,id=row.names)
sig_counts<-merge(sig_counts,F_diaphanus_counts,by="id")
dim(sig_counts)
# F_parvapinis
F_parvapinis_counts <- as.data.frame(F_parvapinis_norm_counts)
row.names <- rownames(F_parvapinis_counts)
F_parvapinis_counts <- cbind(F_parvapinis_counts,id=row.names)
sig_counts<-merge(sig_counts,F_parvapinis_counts,by="id")
dim(sig_counts)
############
# Fresh
###########
# F_olivaceous
F_olivaceous_counts<-as.data.frame(F_olivaceous_norm_counts)
row.names<-rownames(F_olivaceous_counts)
F_olivaceous_counts<-cbind(F_olivaceous_counts,id=row.names)
sig_counts<-merge(sig_counts,F_olivaceous_counts,by="id")
dim(sig_counts)
# F_rathbuni
F_rathbuni_counts<-as.data.frame(F_rathbuni_norm_counts)
row.names<-rownames(F_rathbuni_counts)
F_rathbuni_counts<-cbind(F_rathbuni_counts,id=row.names)
sig_counts<-merge(sig_counts,F_rathbuni_counts,by="id")
dim(sig_counts)
# L_goodei
L_goodei_counts<-as.data.frame(L_goodei_norm_counts)
row.names<-rownames(L_goodei_counts)
L_goodei_counts<-cbind(L_goodei,id=row.names)
sig_counts<-merge(sig_counts,L_goodei_counts,by="id")
dim(sig_counts)
head(sig_counts)


rownames(sig_counts)<-sig_counts$id
sig_counts<-sig_counts[,c(2:85)]
sig_counts_small<-sig_counts[rownames(sig_counts) %in% overlap_BW_FW_transfer,]



#
d <- as.matrix(sig_counts_small)
d<-na.omit(d)
hr <- hclust(as.dist(1-cor(t(d), method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- greenred(75)
png(filename="heatmap_separate.png",width=3.25,height=3.25,units="in",res=1200,pointsize = 4)
heatmap.2(d, main="Killifish (M,B,F), union of padj<0.05, log2FC +-1", 
          Rowv=as.dendrogram(hr),
          cexRow=0.45,cexCol=0.45,srtCol= 90,
          adjCol = c(NA,0),offsetCol=2,offsetRow=0.1,
          Colv=NA, dendrogram="row", 
          scale="row", col=myheatcol, 
          density.info="none", 
          trace="none", RowSideColors= myClusterSideBar)
dev.off()
