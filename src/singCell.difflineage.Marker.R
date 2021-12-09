#' get top markers for Epi, Endoderm, TE, Meso, Amnion
# R4.0
rm(list=ls())
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(scran))
suppressMessages(library(batchelor))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))
suppressMessages(library(readxl))

# working directory
DIR <- "~/My_project/CheckBlastoids"
setwd(DIR)

#' Loading R functions
#source("~/PC/R_code/functions.R")
#source("~/PC/SnkM/SgCell.R")
feature.plot.cor <- c("yellow","red","black")


options(digits = 4)
options(future.globals.maxSize= 3001289600)

suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)

TD="Nov14_2021"
MM="Pub"
zs.limit <- 2.5
rename <- dplyr::rename
#' loading dataset
load("tmp_data/gene.meta.Rdata",verbose=T)
meta.filter<- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds"))
counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/",MM,".lognormExp.mBN.rds"))
data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/Pub.D2.umap.cord.rds"))



# meta data
temp.M <-  data.ob.umap %>% filter(cluster_EML %in% c("Epiblast","ELC","TE","TLC","Endoderm","HLC","Mesoderm","MeLC","Amnion","AMLC"))
temp.M.smt2 <- temp.M %>% filter(seqType=="smt2") 
temp.M.10X <- temp.M %>% filter(seqType=="10X")


#' select expressed genes
expG.set <- list()
expG.set$smt2 <- rownames(counts.filter )[rowSums(counts.filter[,temp.M.smt2$cell] >=1) >=5] # the most frequent expressed gene
expG.set$X10 <- rownames(counts.filter )[rowSums(counts.filter[,temp.M.10X$cell] >=1) >=5]# the most frequent expressed gene
temp.expG <- intersect(expG.set$smt2,expG.set$X10) %>% intersect(rownames(lognormExp.mBN))


#' for smart-seq dataset only
data.temp <- CreateSeuratObject(counts.filter[temp.expG ,temp.M.smt2 %>% filter(cluster_EML %in% c("Epiblast","TE","Endoderm","Mesoderm","Amnion") ) %>%pull(cell)], meta.data = (temp.M.smt2 %>% filter(cluster_EML %in% c("Epiblast","TE","Endoderm","Mesoderm","Amnion")) %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) 
data.temp@assays$RNA@data <- as.matrix(lognormExp.mBN[rownames(data.temp),colnames(data.temp)])

Idents(data.temp) <- as.factor(data.temp@meta.data$cluster_EML)

#de.out <- FindAllMarkers(data.temp,test.use="roc",verbose=F)  %>% tbl_df() %>%rename(Gene=gene)
lineage.mk <- FunRF_FindAllMarkers_para(data.temp)
#save(de.out,temp.DM, file=paste0("tmp_data/",TD,"/",MM,".DE.allL.Rdata"))
saveRDS(lineage.mk, file=paste0("tmp_data/",TD,"/",MM,".mk.allL.rds"))



