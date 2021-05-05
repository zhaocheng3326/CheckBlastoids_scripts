#' ---
#' title: "select markers for each lineages"
#' output: 
#'  html_document:
#'    code_folding: hide
#'    
#' ---


#  R3.6

#'
# loading R library
#R3.6 

#' select the Marker genes to distinct TE and Amnion
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
DIR <- "~/My_project/JP_project"
setwd(DIR)

#' Loading R functions
source("~/PC/R_code/functions.R")
source("~/PC/SnkM/SgCell.R")
feature.plot.cor <- c("yellow","red","black")


options(digits = 4)
options(future.globals.maxSize= 3001289600)

suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)

TD="D2_pub"
zs.limit <- 2.5
rename <- dplyr::rename

#' loading dataset
load("tmp_data/gene.meta.Rdata",verbose=T)
meta.filter<- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds"))
counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/lognormExp.mBN.rds"))
data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/Pub.D2.umap.cord.rds"))


#' detect markers
temp.M <-  data.ob.umap %>% filter(new_EML %in% c("Epiblast","ELC","TE","TLC","Endoderm","HLC","Mesoderm","MeLC","Amnion","AMLC"))
temp.M.smt2 <- temp.M %>% filter(seqType=="smt2")
temp.M.10X <- temp.M %>% filter(seqType=="10X")

expG.set <- list()

expG.set$smt2 <- rownames(counts.filter )[rowSums(counts.filter[,temp.M.smt2$cell] >=1) >=5] # the most frequent expressed gene
expG.set$X10 <- rownames(counts.filter )[rowSums(counts.filter[,temp.M.10X$cell] >=1) >=5]# the most frequent expressed gene

temp.expG <- intersect(expG.set$smt2,expG.set$X10) %>% intersect(rownames(lognormExp.mBN))

#' for smart-seq dataset only
data.temp <- CreateSeuratObject(counts.filter[temp.expG ,temp.M.smt2$cell], meta.data = (temp.M.smt2 %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) 
data.temp@assays$RNA@data <- as.matrix(lognormExp.mBN[rownames(data.temp),colnames(data.temp)])

Idents(data.temp) <- as.factor(data.temp@meta.data$new_EML)

de.out <- FindAllMarkers(data.temp,test.use="roc",verbose=F)  %>% tbl_df() %>%rename(Gene=gene)
temp.DM <- FunRF_FindAllMarkers_para(data.temp)
save(de.out,temp.DM, file=paste0("tmp_data/",TD,"/DE.allL.Rdata"))



