library(DESeq2)
library(RColorBrewer)
library(gplots)
source('~/Documents/scripts/plotPCAWithSampleNames.R')
source('~/Documents/scripts/overLapper_original.R')
setwd("~/Documents/UCDavis/osmotic")

data.1<-read.csv("killifish_allcounts.csv")
annotation<-read.table("kfish2rae5g.annotation.transcript.name.id", fill=TRUE,header=FALSE)
colnames(annotation)<-c("id","gene")
head(data.1)
colnames(data.1)
id<-data.1$GeneID
rownames(data.1)<-id
head(data.1)

######
##F_similis
######
F_similis<-data.1[,c(99:101,103:107)]
colnames(F_similis)
col.names<-colnames(F_similis)
head(F_similis)
conditions = sapply(strsplit(col.names,"_"),`[`,4)
genus = sapply(strsplit(col.names,"_"),`[`,1)
species = sapply(strsplit(col.names,"_"),`[`,2)
genus_species = paste(genus,species,sep="_")
pop = sapply(strsplit(col.names,"_"),`[`,3)
genus_species_pop = paste(genus_species,pop,sep=".")
genus_species = gsub(".NA", "", genus_species_pop)
ExpDesign <- data.frame(row.names=colnames(F_similis), condition = conditions,genus_species = genus_species)
ExpDesign
cds<-DESeqDataSetFromMatrix(countData=F_similis, 
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

# get normalized counts
# add id column
F_similis_norm_counts<-counts(cds,normalized=TRUE)
id<-rownames(F_similis_norm_counts)
F_similis_norm_counts<-cbind(F_similis_norm_counts,id)

# merge res1, res2, res3 with counts
# "BW","FW"
res1_df<-as.data.frame(res.1)
colnames(res1_df)<-paste(colnames(res1_df),"BW_FW", sep='.')
id<-rownames(res1_df)
res1_df<-cbind(res1_df,id)
dim(res1_df)
# "transfer","FW"
res2_df<-as.data.frame(res.2)
colnames(res2_df)<-paste(colnames(res2_df),"transfer_FW", sep='.')
id<-rownames(res2_df)
res2_df<-cbind(res2_df,id)
dim(res2_df)
# "transfer","BW"
res3_df<-as.data.frame(res.3)
colnames(res3_df)<-paste(colnames(res3_df),"transfer_BW", sep='.')
id<-rownames(res3_df)
res3_df<-cbind(res3_df,id)
dim(res3_df)
F_similis_res<-merge(F_similis_norm_counts,res1_df,by="id")
dim(F_similis_res)
colnames(F_similis_res)
F_similis_res<-merge(F_similis_res,res2_df,by="id")
dim(F_similis_res)
colnames(F_similis_res)
F_similis_res<-merge(F_similis_res,res3_df,by="id")
dim(F_similis_res)
colnames(F_similis_res)
head(F_similis_res)
# remove rows with NA
F_similis_res<-F_similis_res[complete.cases(F_similis_res),]
dim(F_similis_res)
F_similis_annotated<-merge(F_similis_res,annotation,by="id")
#F_similis_annotated<-F_similis_annotated[,c(ncol(F_similis_annotated),1:(ncol(F_similis_annotated)-1))]
write.csv(F_similis_annotated,"F_similis_results_all.csv")

plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="F. similis (BW vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="F. similis (transfer vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.3$baseMean), res.3$log2FoldChange, 
     col=ifelse(res.3$padj < 0.05, "red","gray67"),
     main="F. similis (transfer vs. BW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
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

# get lists of unique genes for each comparison
F_similis_BW_FW<-OLlist$Venn_List$BW_FW
length(F_similis_BW_FW)
F_similis_transfer_FW<-OLlist$Venn_List$transfer_FW
length(F_similis_transfer_FW)
F_similis_transfer_BW<-OLlist$Venn_List$transfer_BW
length(F_similis_transfer_BW)



######
##F_heteroclitus.MDPL
######
F_heteroclitus.MDPL<-data.1[,c(41:49)]
colnames(F_heteroclitus.MDPL)
col.names<-colnames(F_heteroclitus.MDPL)
head(F_heteroclitus.MDPL)
conditions = sapply(strsplit(col.names,"_"),`[`,4)
genus = sapply(strsplit(col.names,"_"),`[`,1)
species = sapply(strsplit(col.names,"_"),`[`,2)
genus_species = paste(genus,species,sep="_")
pop = sapply(strsplit(col.names,"_"),`[`,3)
genus_species_pop = paste(genus_species,pop,sep=".")
genus_species = gsub(".NA", "", genus_species_pop)
ExpDesign <- data.frame(row.names=colnames(F_heteroclitus.MDPL), condition = conditions,genus_species = genus_species)
ExpDesign
cds<-DESeqDataSetFromMatrix(countData=F_heteroclitus.MDPL, 
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

# get normalized counts
# add id column
F_heteroclitus.MDPL_norm_counts<-counts(cds,normalized=TRUE)
id<-rownames(F_heteroclitus.MDPL_norm_counts)
F_heteroclitus.MDPL_norm_counts<-cbind(F_heteroclitus.MDPL_norm_counts,id)

# merge res1, res2, res3 with counts
# "BW","FW"
res1_df<-as.data.frame(res.1)
colnames(res1_df)<-paste(colnames(res1_df),"BW_FW", sep='.')
id<-rownames(res1_df)
res1_df<-cbind(res1_df,id)
dim(res1_df)
# "transfer","FW"
res2_df<-as.data.frame(res.2)
colnames(res2_df)<-paste(colnames(res2_df),"transfer_FW", sep='.')
id<-rownames(res2_df)
res2_df<-cbind(res2_df,id)
dim(res2_df)
# "transfer","BW"
res3_df<-as.data.frame(res.3)
colnames(res3_df)<-paste(colnames(res3_df),"transfer_BW", sep='.')
id<-rownames(res3_df)
res3_df<-cbind(res3_df,id)
dim(res3_df)
F_heteroclitus.MDPL_res<-merge(F_heteroclitus.MDPL_norm_counts,res1_df,by="id")
dim(F_heteroclitus.MDPL_res)
colnames(F_heteroclitus.MDPL_res)
F_heteroclitus.MDPL_res<-merge(F_heteroclitus.MDPL_res,res2_df,by="id")
dim(F_heteroclitus.MDPL_res)
colnames(F_heteroclitus.MDPL_res)
F_heteroclitus.MDPL_res<-merge(F_heteroclitus.MDPL_res,res3_df,by="id")
dim(F_heteroclitus.MDPL_res)
colnames(F_heteroclitus.MDPL_res)
F_heteroclitus.MDPL_res<-F_heteroclitus.MDPL_res[complete.cases(F_heteroclitus.MDPL_res),]
dim(F_heteroclitus.MDPL_res)
F_heteroclitus.MDPL_annotated<-merge(F_heteroclitus.MDPL_res,annotation,by="id")
F_heteroclitus.MDPL_annotated<-F_heteroclitus.MDPL_annotated[,c(ncol(F_heteroclitus.MDPL_annotated),1:(ncol(F_heteroclitus.MDPL_annotated)-1))]
#write.csv(F_heteroclitus.MDPL_annotated,"F_heteroclitus.MDPL_results_all.csv")

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

# get lists of unique genes for each comparison
F_heteroclitus.MDPL_BW_FW<-OLlist$Venn_List$BW_FW
length(F_heteroclitus.MDPL_BW_FW)
F_heteroclitus.MDPL_transfer_FW<-OLlist$Venn_List$transfer_FW
length(F_heteroclitus.MDPL_transfer_FW)
F_heteroclitus.MDPL_transfer_BW<-OLlist$Venn_List$transfer_BW
length(F_heteroclitus.MDPL_transfer_BW)



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


# get normalized counts
# add id column
F_heteroclitus.MDPP_norm_counts<-counts(cds,normalized=TRUE)
id<-rownames(F_heteroclitus.MDPP_norm_counts)
F_heteroclitus.MDPP_norm_counts<-cbind(F_heteroclitus.MDPP_norm_counts,id)

# merge res1, res2, res3 with counts
# "BW","FW"
res1_df<-as.data.frame(res.1)
colnames(res1_df)<-paste(colnames(res1_df),"BW_FW", sep='.')
id<-rownames(res1_df)
res1_df<-cbind(res1_df,id)
# "transfer","FW"
res2_df<-as.data.frame(res.2)
colnames(res2_df)<-paste(colnames(res2_df),"transfer_FW", sep='.')
id<-rownames(res2_df)
res2_df<-cbind(res2_df,id)
# "transfer","BW"
res3_df<-as.data.frame(res.3)
colnames(res3_df)<-paste(colnames(res3_df),"transfer_BW", sep='.')
id<-rownames(res3_df)
res3_df<-cbind(res3_df,id)
F_heteroclitus.MDPP_res<-merge(F_heteroclitus.MDPP_norm_counts,res1_df,by="id")
F_heteroclitus.MDPP_res<-merge(F_heteroclitus.MDPP_res,res2_df,by="id")
F_heteroclitus.MDPP_res<-merge(F_heteroclitus.MDPP_res,res3_df,by="id")
dim(F_heteroclitus.MDPP_res)
F_heteroclitus.MDPP_res<-F_heteroclitus.MDPP_res[complete.cases(F_heteroclitus.MDPP_res),]
dim(F_heteroclitus.MDPP_res)
F_heteroclitus.MDPP_annotated<-merge(F_heteroclitus.MDPP_res,annotation,by="id")
F_heteroclitus.MDPP_annotated<-F_heteroclitus.MDPP_annotated[,c(ncol(F_heteroclitus.MDPP_annotated),1:(ncol(F_heteroclitus.MDPP_annotated)-1))]
write.csv(F_heteroclitus.MDPP_annotated,"F_heteroclitus.MDPP_results_all.csv")

plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPP (BW vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPP (transfer vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.3$baseMean), res.3$log2FoldChange, 
     col=ifelse(res.3$padj < 0.05, "red","gray67"),
     main="F. heteroclitus.MDPP (transfer vs. BW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
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

# get lists of unique genes for each comparison
F_heteroclitus.MDPP_BW_FW<-OLlist$Venn_List$BW_FW
length(F_heteroclitus.MDPP_BW_FW)
F_heteroclitus.MDPP_transfer_FW<-OLlist$Venn_List$transfer_FW
length(F_heteroclitus.MDPP_transfer_FW)
F_heteroclitus.MDPP_transfer_BW<-OLlist$Venn_List$transfer_BW
length(F_heteroclitus.MDPP_transfer_BW)

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

# get normalized counts
# add id column
F_chrysotus_norm_counts<-counts(cds,normalized=TRUE)
id<-rownames(F_chrysotus_norm_counts)
F_chrysotus_norm_counts<-cbind(F_chrysotus_norm_counts,id)

# merge res1, res2, res3 with counts
# "BW","FW"
res1_df<-as.data.frame(res.1)
colnames(res1_df)<-paste(colnames(res1_df),"BW_FW", sep='.')
id<-rownames(res1_df)
res1_df<-cbind(res1_df,id)
# "transfer","FW"
res2_df<-as.data.frame(res.2)
colnames(res2_df)<-paste(colnames(res2_df),"transfer_FW", sep='.')
id<-rownames(res2_df)
res2_df<-cbind(res2_df,id)
# "transfer","BW"
res3_df<-as.data.frame(res.3)
colnames(res3_df)<-paste(colnames(res3_df),"transfer_BW", sep='.')
id<-rownames(res3_df)
res3_df<-cbind(res3_df,id)
F_chrysotus_res<-merge(F_chrysotus_norm_counts,res1_df,by="id")
F_chrysotus_res<-merge(F_chrysotus_res,res2_df,by="id")
F_chrysotus_res<-merge(F_chrysotus_res,res3_df,by="id")
dim(F_chrysotus_res)
F_chrysotus_res<-F_chrysotus_res[complete.cases(F_chrysotus_res),]
dim(F_chrysotus_res)
F_chrysotus_annotated<-merge(F_chrysotus_res,annotation,by="id")
F_chrysotus_annotated<-F_chrysotus_annotated[,c(ncol(F_chrysotus_annotated),1:(ncol(F_chrysotus_annotated)-1))]
write.csv(F_chrysotus_annotated,"F_chrysotus_results_all.csv")

plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="F_chrysotus (BW vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="F_chrysotus (transfer vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
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

# get lists of unique genes for each comparison
F_chrysotus_BW_FW<-OLlist$Venn_List$BW_FW
length(F_chrysotus_BW_FW)
F_chrysotus_transfer_FW<-OLlist$Venn_List$transfer_FW
length(F_chrysotus_transfer_FW)
F_chrysotus_transfer_BW<-OLlist$Venn_List$transfer_BW
length(F_chrysotus_transfer_BW)

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

# get normalized counts
# add id column
F_diaphanus_norm_counts<-counts(cds,normalized=TRUE)
id<-rownames(F_diaphanus_norm_counts)
F_diaphanus_norm_counts<-cbind(F_diaphanus_norm_counts,id)

# merge res1, res2, res3 with counts
# "BW","FW"
res1_df<-as.data.frame(res.1)
colnames(res1_df)<-paste(colnames(res1_df),"BW_FW", sep='.')
id<-rownames(res1_df)
res1_df<-cbind(res1_df,id)
# "transfer","FW"
res2_df<-as.data.frame(res.2)
colnames(res2_df)<-paste(colnames(res2_df),"transfer_FW", sep='.')
id<-rownames(res2_df)
res2_df<-cbind(res2_df,id)
# "transfer","BW"
res3_df<-as.data.frame(res.3)
colnames(res3_df)<-paste(colnames(res3_df),"transfer_BW", sep='.')
id<-rownames(res3_df)
res3_df<-cbind(res3_df,id)
F_diaphanus_res<-merge(F_diaphanus_norm_counts,res1_df,by="id")
F_diaphanus_res<-merge(F_diaphanus_res,res2_df,by="id")
F_diaphanus_res<-merge(F_diaphanus_res,res3_df,by="id")
dim(F_diaphanus_res)
F_diaphanus_res<-F_diaphanus_res[complete.cases(F_diaphanus_res),]
dim(F_diaphanus_res)
F_diaphanus_annotated<-merge(F_diaphanus_res,annotation,by="id")
F_diaphanus_annotated<-F_diaphanus_annotated[,c(ncol(F_diaphanus_annotated),1:(ncol(F_diaphanus_annotated)-1))]
write.csv(F_diaphanus_annotated,"F_diaphanus_results_all.csv")

plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="F. diaphanus (BW vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="F. diaphanus (transfer vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.3$baseMean), res.3$log2FoldChange, 
     col=ifelse(res.3$padj < 0.05, "red","gray67"),
     main="F. diaphanus (transfer vs. BW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
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

# get lists of unique genes for each comparison
F_diaphanus_BW_FW<-OLlist$Venn_List$BW_FW
length(F_diaphanus_BW_FW)
F_diaphanus_transfer_FW<-OLlist$Venn_List$transfer_FW
length(F_diaphanus_transfer_FW)
F_diaphanus_transfer_BW<-OLlist$Venn_List$transfer_BW
length(F_diaphanus_transfer_BW)

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

# get normalized counts
# add id column
F_grandis_norm_counts<-counts(cds,normalized=TRUE)
id<-rownames(F_grandis_norm_counts)
F_grandis_norm_counts<-cbind(F_grandis_norm_counts,id)

# merge res1, res2, res3 with counts
# "BW","FW"
res1_df<-as.data.frame(res.1)
colnames(res1_df)<-paste(colnames(res1_df),"BW_FW", sep='.')
id<-rownames(res1_df)
res1_df<-cbind(res1_df,id)
# "transfer","FW"
res2_df<-as.data.frame(res.2)
colnames(res2_df)<-paste(colnames(res2_df),"transfer_FW", sep='.')
id<-rownames(res2_df)
res2_df<-cbind(res2_df,id)
# "transfer","BW"
res3_df<-as.data.frame(res.3)
colnames(res3_df)<-paste(colnames(res3_df),"transfer_BW", sep='.')
id<-rownames(res3_df)
res3_df<-cbind(res3_df,id)
F_grandis_res<-merge(F_grandis_norm_counts,res1_df,by="id")
F_grandis_res<-merge(F_grandis_res,res2_df,by="id")
F_grandis_res<-merge(F_grandis_res,res3_df,by="id")
dim(F_grandis_res)
F_grandis_res<-F_grandis_res[complete.cases(F_grandis_res),]
dim(F_grandis_res)
F_grandis_annotated<-merge(F_grandis_res,annotation,by="id")
F_grandis_annotated<-F_grandis_annotated[,c(ncol(F_grandis_annotated),1:(ncol(F_grandis_annotated)-1))]
write.csv(F_grandis_annotated,"F_grandis_results_all.csv")

plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="F. grandis (BW vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="F. grandis (transfer vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
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

# get lists of unique genes for each comparison
F_grandis_BW_FW<-OLlist$Venn_List$BW_FW
length(F_grandis_BW_FW)
F_grandis_transfer_FW<-OLlist$Venn_List$transfer_FW
length(F_grandis_transfer_FW)
F_grandis_transfer_BW<-OLlist$Venn_List$transfer_BW
length(F_grandis_transfer_BW)

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

# get normalized counts
# add id column
F_olivaceous_norm_counts<-counts(cds,normalized=TRUE)
id<-rownames(F_olivaceous_norm_counts)
F_olivaceous_norm_counts<-cbind(F_olivaceous_norm_counts,id)

# merge res1, res2, res3 with counts
# "BW","FW"
res1_df<-as.data.frame(res.1)
colnames(res1_df)<-paste(colnames(res1_df),"BW_FW", sep='.')
id<-rownames(res1_df)
res1_df<-cbind(res1_df,id)
# "transfer","FW"
res2_df<-as.data.frame(res.2)
colnames(res2_df)<-paste(colnames(res2_df),"transfer_FW", sep='.')
id<-rownames(res2_df)
res2_df<-cbind(res2_df,id)
# "transfer","BW"
res3_df<-as.data.frame(res.3)
colnames(res3_df)<-paste(colnames(res3_df),"transfer_BW", sep='.')
id<-rownames(res3_df)
res3_df<-cbind(res3_df,id)
F_olivaceous_res<-merge(F_olivaceous_norm_counts,res1_df,by="id")
F_olivaceous_res<-merge(F_olivaceous_res,res2_df,by="id")
F_olivaceous_res<-merge(F_olivaceous_res,res3_df,by="id")
dim(F_olivaceous_res)
F_olivaceous_res<-F_olivaceous_res[complete.cases(F_olivaceous_res),]
dim(F_olivaceous_res)
F_olivaceous_annotated<-merge(F_olivaceous_res,annotation,by="id")
F_olivaceous_annotated<-F_olivaceous_annotated[,c(ncol(F_olivaceous_annotated),1:(ncol(F_olivaceous_annotated)-1))]
write.csv(F_olivaceous_annotated,"F_olivaceous_results_all.csv")

plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="F. olivaceous (BW vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="F. olivaceous (transfer vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.3$baseMean), res.3$log2FoldChange, 
     col=ifelse(res.3$padj < 0.05, "red","gray67"),
     main="F. olivaceous (transfer vs. BW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
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

# get lists of unique genes for each comparison
F_olivaceous_BW_FW<-OLlist$Venn_List$BW_FW
length(F_olivaceous_BW_FW)
F_olivaceous_transfer_FW<-OLlist$Venn_List$transfer_FW
length(F_olivaceous_transfer_FW)
F_olivaceous_transfer_BW<-OLlist$Venn_List$transfer_BW
length(F_olivaceous_transfer_BW)

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

# get normalized counts
# add id column
F_parvapinis_norm_counts<-counts(cds,normalized=TRUE)
id<-rownames(F_parvapinis_norm_counts)
F_parvapinis_norm_counts<-cbind(F_parvapinis_norm_counts,id)

# merge res1, res2, res3 with counts
# "BW","FW"
res1_df<-as.data.frame(res.1)
colnames(res1_df)<-paste(colnames(res1_df),"BW_FW", sep='.')
id<-rownames(res1_df)
res1_df<-cbind(res1_df,id)
# "transfer","FW"
res2_df<-as.data.frame(res.2)
colnames(res2_df)<-paste(colnames(res2_df),"transfer_FW", sep='.')
id<-rownames(res2_df)
res2_df<-cbind(res2_df,id)
# "transfer","BW"
res3_df<-as.data.frame(res.3)
colnames(res3_df)<-paste(colnames(res3_df),"transfer_BW", sep='.')
id<-rownames(res3_df)
res3_df<-cbind(res3_df,id)
F_parvapinis_res<-merge(F_parvapinis_norm_counts,res1_df,by="id")
F_parvapinis_res<-merge(F_parvapinis_res,res2_df,by="id")
F_parvapinis_res<-merge(F_parvapinis_res,res3_df,by="id")
dim(F_parvapinis_res)
F_parvapinis_res<-F_parvapinis_res[complete.cases(F_parvapinis_res),]
dim(F_parvapinis_res)
F_parvapinis_annotated<-merge(F_parvapinis_res,annotation,by="id")
F_parvapinis_annotated<-F_parvapinis_annotated[,c(ncol(F_parvapinis_annotated),1:(ncol(F_parvapinis_annotated)-1))]
write.csv(F_parvapinis_annotated,"F_parvapinis_results_all.csv")

plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="F. parvapinis (BW vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="F. parvapinis (transfer vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.3$baseMean), res.3$log2FoldChange, 
     col=ifelse(res.3$padj < 0.05, "red","gray67"),
     main="F. parvapinis (transfer vs. BW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
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

# get lists of unique genes for each comparison
F_parvapinis_BW_FW<-OLlist$Venn_List$BW_FW
length(F_parvapinis_BW_FW)
F_parvapinis_transfer_FW<-OLlist$Venn_List$transfer_FW
length(F_parvapinis_transfer_FW)
F_parvapinis_transfer_BW<-OLlist$Venn_List$transfer_BW
length(F_parvapinis_transfer_BW)

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


# get normalized counts
# add id column
F_rathbuni_norm_counts<-counts(cds,normalized=TRUE)
id<-rownames(F_rathbuni_norm_counts)
F_rathbuni_norm_counts<-cbind(F_rathbuni_norm_counts,id)

# merge res1, res2, res3 with counts
# "BW","FW"
res1_df<-as.data.frame(res.1)
colnames(res1_df)<-paste(colnames(res1_df),"BW_FW", sep='.')
id<-rownames(res1_df)
res1_df<-cbind(res1_df,id)
# "transfer","FW"
res2_df<-as.data.frame(res.2)
colnames(res2_df)<-paste(colnames(res2_df),"transfer_FW", sep='.')
id<-rownames(res2_df)
res2_df<-cbind(res2_df,id)
# "transfer","BW"
res3_df<-as.data.frame(res.3)
colnames(res3_df)<-paste(colnames(res3_df),"transfer_BW", sep='.')
id<-rownames(res3_df)
res3_df<-cbind(res3_df,id)
F_rathbuni_res<-merge(F_rathbuni_norm_counts,res1_df,by="id")
F_rathbuni_res<-merge(F_rathbuni_res,res2_df,by="id")
F_rathbuni_res<-merge(F_rathbuni_res,res3_df,by="id")
dim(F_rathbuni_res)
F_rathbuni_res<-F_rathbuni_res[complete.cases(F_rathbuni_res),]
dim(F_rathbuni_res)
F_rathbuni_annotated<-merge(F_rathbuni_res,annotation,by="id")
F_rathbuni_annotated<-F_rathbuni_annotated[,c(ncol(F_rathbuni_annotated),1:(ncol(F_rathbuni_annotated)-1))]
write.csv(F_rathbuni_annotated,"F_rathbuni_results_all.csv")

plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="F. rathbuni (BW vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="F. rathbuni (transfer vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.3$baseMean), res.3$log2FoldChange, 
     col=ifelse(res.3$padj < 0.05, "red","gray67"),
     main="F. rathbuni (transfer vs. BW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
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

# get lists of unique genes for each comparison
F_rathbuni_BW_FW<-OLlist$Venn_List$BW_FW
length(F_rathbuni_BW_FW)
F_rathbuni_transfer_FW<-OLlist$Venn_List$transfer_FW
length(F_rathbuni_transfer_FW)
F_rathbuni_transfer_BW<-OLlist$Venn_List$transfer_BW
length(F_rathbuni_transfer_BW)

#######
##L_goodei
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

# get normalized counts
# add id column
L_goodei_norm_counts<-counts(cds,normalized=TRUE)
id<-rownames(L_goodei_norm_counts)
L_goodei_norm_counts<-cbind(L_goodei_norm_counts,id)

# merge res1, res2, res3 with counts
# "BW","FW"
res1_df<-as.data.frame(res.1)
colnames(res1_df)<-paste(colnames(res1_df),"BW_FW", sep='.')
id<-rownames(res1_df)
res1_df<-cbind(res1_df,id)
# "transfer","FW"
res2_df<-as.data.frame(res.2)
colnames(res2_df)<-paste(colnames(res2_df),"transfer_FW", sep='.')
id<-rownames(res2_df)
res2_df<-cbind(res2_df,id)
# "transfer","BW"
res3_df<-as.data.frame(res.3)
colnames(res3_df)<-paste(colnames(res3_df),"transfer_BW", sep='.')
id<-rownames(res3_df)
res3_df<-cbind(res3_df,id)
L_goodei_res<-merge(L_goodei_norm_counts,res1_df,by="id")
L_goodei_res<-merge(L_goodei_res,res2_df,by="id")
L_goodei_res<-merge(L_goodei_res,res3_df,by="id")
dim(L_goodei_res)
L_goodei_res<-L_goodei_res[complete.cases(L_goodei_res),]
dim(L_goodei_res)
L_goodei_annotated<-merge(L_goodei_res,annotation,by="id")
L_goodei_annotated<-L_goodei_annotated[,c(ncol(L_goodei_annotated),1:(ncol(L_goodei_annotated)-1))]
write.csv(L_goodei_annotated,"L_goodei_results_all.csv")


plot(log2(res.1$baseMean), res.1$log2FoldChange, 
     col=ifelse(res.1$padj < 0.05, "red","gray67"),
     main="L. goodei (BW vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.2$baseMean), res.2$log2FoldChange, 
     col=ifelse(res.2$padj < 0.05, "red","gray67"),
     main="L. goodei (transfer vs. FW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
abline(h=c(-1,1), col="blue")
plot(log2(res.3$baseMean), res.3$log2FoldChange, 
     col=ifelse(res.3$padj < 0.05, "red","gray67"),
     main="L. goodei (transfer vs. BW) (padj<0.05)",xlim=c(1,15),pch=20,cex=1)
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

# get lists of unique genes for each comparison
L_goodei_BW_FW<-OLlist$Venn_List$BW_FW
length(L_goodei_BW_FW)
L_goodei_transfer_FW<-OLlist$Venn_List$transfer_FW
length(L_goodei_transfer_FW)
L_goodei_transfer_BW<-OLlist$Venn_List$transfer_BW
length(L_goodei_transfer_BW)

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
# get normalized counts
# add id column
L_parva_norm_counts<-counts(cds,normalized=TRUE)
id<-rownames(L_parva_norm_counts)
L_parva_norm_counts<-cbind(L_parva_norm_counts,id)

# merge res1, res2, res3 with counts
# "BW","FW"
res1_df<-as.data.frame(res.1)
colnames(res1_df)<-paste(colnames(res1_df),"BW_FW", sep='.')
id<-rownames(res1_df)
res1_df<-cbind(res1_df,id)
# "transfer","FW"
res2_df<-as.data.frame(res.2)
colnames(res2_df)<-paste(colnames(res2_df),"transfer_FW", sep='.')
id<-rownames(res2_df)
res2_df<-cbind(res2_df,id)
# "transfer","BW"
res3_df<-as.data.frame(res.3)
colnames(res3_df)<-paste(colnames(res3_df),"transfer_BW", sep='.')
id<-rownames(res3_df)
res3_df<-cbind(res3_df,id)
L_parva_res<-merge(L_parva_norm_counts,res1_df,by="id")
L_parva_res<-merge(L_parva_res,res2_df,by="id")
L_parva_res<-merge(L_parva_res,res3_df,by="id")
dim(L_parva_res)
L_parva_res<-L_parva_res[complete.cases(L_parva_res),]
dim(L_parva_res)
L_parva_annotated<-merge(L_parva_res,annotation,by="id")
L_parva_annotated<-L_parva_annotated[,c(ncol(L_parva_annotated),1:(ncol(L_parva_annotated)-1))]
write.csv(L_parva_annotated,"L_parva_results_all.csv")

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

# get lists of unique genes for each comparison
L_parva_BW_FW<-OLlist$Venn_List$BW_FW
length(L_parva_BW_FW)
L_parva_transfer_FW<-OLlist$Venn_List$transfer_FW
length(L_parva_transfer_FW)
L_parva_transfer_BW<-OLlist$Venn_List$transfer_BW
length(L_parva_transfer_BW)


# get list of unique results from each species
# BW_FW
length(F_heteroclitus.MDPL_BW_FW)
length(F_heteroclitus.MDPP_BW_FW)
length(F_chrysotus_BW_FW)
length(F_diaphanus_BW_FW)
length(F_grandis_BW_FW)
length(F_olivaceous_BW_FW)
length(F_parvapinis_BW_FW)
length(F_rathbuni_BW_FW)
length(L_goodei_BW_FW)
length(L_parva_BW_FW)
# transfer_FW
length(F_heteroclitus.MDPL_transfer_FW)
length(F_heteroclitus.MDPP_transfer_FW)
length(F_chrysotus_transfer_FW)
length(F_diaphanus_transfer_FW)
length(F_grandis_transfer_FW)
length(F_olivaceous_transfer_FW)
length(F_parvapinis_transfer_FW)
length(F_rathbuni_transfer_FW)
length(L_goodei_transfer_FW)
length(L_parva_transfer_FW)
# transfer_BW
length(F_heteroclitus.MDPL_transfer_BW)
length(F_heteroclitus.MDPP_transfer_BW)
length(F_chrysotus_transfer_BW)
length(F_diaphanus_transfer_BW)
length(F_grandis_transfer_BW)
length(F_olivaceous_transfer_BW)
length(F_parvapinis_transfer_BW)
length(F_rathbuni_transfer_BW)
length(L_goodei_transfer_BW)
length(L_parva_transfer_BW)

all_transfer_BW <- list(F_heteroclitus.MDPL_transfer_BW,F_heteroclitus.MDPP_transfer_BW,F_chrysotus_transfer_BW,
                        F_diaphanus_transfer_BW, F_grandis_transfer_BW,F_olivaceous_transfer_BW,F_parvapinis_transfer_BW,
                        F_rathbuni_transfer_BW,L_goodei_transfer_BW,L_parva_transfer_BW)

all_transfer_FW <- list(F_heteroclitus.MDPL_transfer_FW,F_heteroclitus.MDPP_transfer_FW,F_chrysotus_transfer_FW,
                        F_diaphanus_transfer_FW,F_grandis_transfer_FW,F_olivaceous_transfer_FW,F_parvapinis_transfer_FW,
                        F_rathbuni_transfer_FW,L_goodei_transfer_FW,L_parva_transfer_FW)
all_BW_FW <- list(F_heteroclitus.MDPL_BW_FW,F_heteroclitus.MDPP_BW_FW,F_chrysotus_BW_FW,F_diaphanus_BW_FW,
                  F_grandis_BW_FW,F_olivaceous_BW_FW,F_parvapinis_BW_FW,F_rathbuni_BW_FW,L_goodei_BW_FW,L_parva_BW_FW)


all_overlap_BW_FW <- Reduce(intersect, all_BW_FW)
all_overlap_transfer_FW <- Reduce(intersect, all_transfer_FW)
all_overlap_transfer_BW <- Reduce(intersect, all_transfer_BW)

m<-all_overlap_BW_FW
length(m)
n<-all_overlap_transfer_FW
length(n)
o<-all_overlap_transfer_BW
length(o)
setlist <- list(BW_FW=as.vector(m),transfer_FW=as.vector(n),transfer_BW=as.vector(o))
OLlist <- overLapper(setlist=setlist, sep="", type="vennsets")
counts <- sapply(OLlist$Venn_List, length)
vennPlot(counts=counts)




overlap_BW_FW_transfer <- OLlist$Venn_List$BW_FWtransfer_FW
overlap_BW_FW_transfer <- intersect(overlap_BW_FW_transfer,OLlist$Venn_List$BW_FWtransfer_BW)
overlap_BW_FW_transfer <- intersect(overlap_BW_FW_transfer,OLlist$Venn_List$transfer_FWtransfer_BW)
overlap_BW_FW_transfer <- intersect(overlap_BW_FW_transfer,OLlist$Venn_List$BW_FWtransfer_FWtransfer_BW)



length(overlap_BW_FW_transfer)
all_overlap<-OLlist$Venn_List$BW_FWtransfer_FWtransfer_BW
length(all_overlap)

names(OLlist$Venn_List)

BW_FW <- OLlist$Venn_List$BW_FW
transfer_FW <- OLlist$Venn_List$transfer_FW
transfer_BW <- OLlist$Venn_List$transfer_BW



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


#BW_FW <- OLlist$Venn_List$BW_FW
#transfer_FW <- OLlist$Venn_List$transfer_FW
#transfer_BW <- OLlist$Venn_List$transfer_BW
length(BW_FW)
length(transfer_FW)
length(transfer_BW)


rownames(sig_counts)<-sig_counts$id
sig_counts<-sig_counts[,c(2:85)]
sig_counts_small<-sig_counts[rownames(sig_counts) %in% overlap_BW_FW_transfer,]

sig_counts_BW_FW<-sig_counts[rownames(sig_counts) %in% BW_FW, ]
sig_counts_transfer_FW<-sig_counts[rownames(sig_counts) %in% transfer_FW,]
sig_counts_transfer_BW<-sig_counts[rownames(sig_counts) %in% transfer_BW,]



# heatmap for sig_counts_BW_FW
d <- as.matrix(sig_counts_BW_FW)
d<-na.omit(d)
hr <- hclust(as.dist(1-cor(t(d), method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- greenred(75)
png(filename="heatmap_BW_FW.png",width=3.25,height=3.25,units="in",res=1200,pointsize = 4)
heatmap.2(d, main="Killifish (M,B,F), intersect of padj<0.05, log2FC +-1", 
          Rowv=as.dendrogram(hr),
          cexRow=0.45,cexCol=0.45,srtCol= 90,
          adjCol = c(NA,0),offsetCol=2,offsetRow=0.1,
          Colv=NA, dendrogram="row", 
          scale="row", col=myheatcol, 
          density.info="none", 
          trace="none", RowSideColors= myClusterSideBar)
dev.off()


## heatmap for sig_counts_transfer_FW
d <- as.matrix(sig_counts_transfer_FW)
d<-na.omit(d)
hr <- hclust(as.dist(1-cor(t(d), method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- greenred(75)
png(filename="heatmap_sig_counts_transfer_FW.png",width=3.25,height=3.25,units="in",res=1200,pointsize = 4)
heatmap.2(d, main="Killifish (M,B,F), intersect of padj<0.05, log2FC +-1", 
          Rowv=as.dendrogram(hr),
          cexRow=0.45,cexCol=0.45,srtCol= 90,
          adjCol = c(NA,0),offsetCol=2,offsetRow=0.1,
          Colv=NA, dendrogram="row", 
          scale="row", col=myheatcol, 
          density.info="none", 
          trace="none", RowSideColors= myClusterSideBar)
dev.off()

# heatmap for sig_counts_transfer_BW
d <- as.matrix(sig_counts_transfer_BW)
d<-na.omit(d)
hr <- hclust(as.dist(1-cor(t(d), method="pearson")), method="complete")
mycl <- cutree(hr, h=max(hr$height/1.5))
clusterCols <- rainbow(length(unique(mycl)))
myClusterSideBar <- clusterCols[mycl]
myheatcol <- greenred(75)
png(filename="heatmap_transfer_BW.png",width=3.25,height=3.25,units="in",res=1200,pointsize = 4)
heatmap.2(d, main="Killifish (M,B,F), intersect of padj<0.05, log2FC +-1", 
          Rowv=as.dendrogram(hr),
          cexRow=0.45,cexCol=0.45,srtCol= 90,
          adjCol = c(NA,0),offsetCol=2,offsetRow=0.1,
          Colv=NA, dendrogram="row", 
          scale="row", col=myheatcol, 
          density.info="none", 
          trace="none", RowSideColors= myClusterSideBar)
dev.off()




