# do this to save workspace
#save(list=ls(),fil="mydata.RData")

library(limma)
library(edgeR)
library(DESeq2)
library(magrittr)
library(stringr)
library(RColorBrewer)
library(dplyr)
source('~/Documents/scripts/plotPCAWithSampleNames.R')
setwd("~/Documents/UCDavis/Whitehead/osmotic/osmotic")
data<-read.csv("allcounts.csv")
head(data)
data[is.na(data)] <- 0
colnames(data)
data.1<-data[,c(3:124)]
#write.csv(data,"killifish_allcounts.csv")
col.names<-colnames(data.1)
col.names
conditions = sapply(strsplit(col.names,"_"),`[`,4)
genus = sapply(strsplit(col.names,"_"),`[`,1)
species = sapply(strsplit(col.names,"_"),`[`,2)
genus_species = paste(genus,species,sep="_")
pop = sapply(strsplit(col.names,"_"),`[`,3)
genus_species_pop = paste(genus_species,pop,sep=".")
genus_species = gsub(".NA", "", genus_species_pop)
marine<-c("A_xenica","F_grandis","F_heteroclitus.MDPP","F_heteroclitus.MDPL","F_similis","L_parva","F_zebrinus")
brackish<-c("F_chrysotus","F_diaphanus","F_parvapinis")
fresh<-c("F_catanatus","F_notatus","F_olivaceous","F_rathbuni","F_sciadicus","L_goodei","F_notti")
# species as a factor
ExpDesign <- data.frame(row.names=colnames(data.1), condition = conditions,genus_species = genus_species)
ExpDesign$physiology = ifelse(ExpDesign$genus_species %in% marine, "M",
  ifelse(ExpDesign$genus_species %in% brackish, "B",
  ifelse(ExpDesign$genus_species %in% fresh,"F",NA)))
cds<-DESeqDataSetFromMatrix(countData=data.1, colData=ExpDesign,design= ~ condition)
mm1 = model.matrix(~ 0 + condition + genus_species + condition:genus_species, ExpDesign)
idx <- which(colSums(mm == 0) == nrow(mm))
mm1 <- mm[,-idx] # removing the columns with all 0's
mm0 <- model.matrix(~ 0 + condition + genus_species, ExpDesign)
design(cds) <- ~ 0 + condition + genus_species + condition:genus_species
cds <- estimateSizeFactors(cds)
cds <- estimateDispersionsGeneEst(cds, modelMatrix=mm1)
#cds$condition <- relevel(cds$condition, "BW")


cds<-DESeq(cds)
log_cds<-rlog(cds)
plotPCAWithSampleNames(log_cds,intgroup="condition",ntop=40000)
norm_counts<-counts(cds,normalized=TRUE)
#norm_counts_log<-counts(log_cds,normalized=TRUE)
res.1<-results(cds,contrast=c("condition","BW","FW"))
res.2<-results(cds,contrast=c("condition","transfer","FW"))
res.3<-results(cds,contrast=c("condition","transfer","BW"))
resultsNames(cds)