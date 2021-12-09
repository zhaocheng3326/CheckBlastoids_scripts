#' ---
#' title: "Integration of all human dataset except for AMLC dataset"
#' output: 
#'  html_document:
#'    code_folding: hide
#'    
#' ---


#R3.6
rm(list=ls())
rewrite=FALSE

condaENV <- "/home/chenzh/miniconda3/envs/R3.6"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratWrappers))


# working directory
DIR <- "~/My_project/CheckBlastoids"
setwd(DIR)
knitr::opts_knit$set(root.dir=DIR)

#' Loading R functions
#source("~/PC/R_code/functions.R")
#source("~/PC/SnkM/SgCell.R")
source("src/local.quick.fun.R")
feature.plot.cor <- c("yellow","red","black")

suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)


options(digits = 4)
options(future.globals.maxSize= 3001289600)
TD="Nov14_2021"
MM="Pub"

nGene <- 2000;pc=30

rename <- dplyr::rename
select<- dplyr::select
filter <- dplyr::filter

#' loading data
load("tmp_data/gene.meta.Rdata",verbose=T)
counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/",MM,".lognormExp.mBN.rds"))

meta.filter<- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% filter(pj %in% c("EBD2","IBD2","D3post","CS7","SPH2016","nBGuo","nicolBla") ) %>% mutate(EML=ifelse(EML %in% c("Ectoderm"),"Amnion",EML))

#temp.M <- meta.filter %>% filter(pj %in% c("IBD2","D3post","CS7","SPH2016","JPF2019"))
temp.M <- meta.filter 
temp.sel.expG <-rownames(lognormExp.mBN )

data.merge <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) 
data.merge@assays$RNA@data <- as.matrix(lognormExp.mBN[temp.sel.expG,colnames(data.merge)])
data.spt <- SplitObject(data.merge, split.by = "pj")%>% lapply(function(x){x=FindVariableFeatures(x,verbose=F,nfeatures=2000)})

# release memory
rm(lognormExp.mBN,data.merge,counts.filter)


data.spt <- data.spt[c("IBD2","EBD2","nicolBla","D3post","CS7","SPH2016","nBGuo")]
set.seed(123)
data.ob <- RunFastMNN(data.spt,verbose=F,features=nGene) %>% RunUMAP( reduction = "mnn", dims = 1:pc,verbose=F) %>% FindNeighbors( reduction = "mnn", dims = 1:pc)
rm(data.spt)

data.temp <- data.ob
data.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(c(cell,SID:mt.perc)) %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") 

#' update sph2016 annotation
cell_lineage_data <- read.delim("~/My_project/ScRNA_analysis/data/Published/E-MTAB-3929.short.txt",stringsAsFactors = F,row.names = 1)%>% tibble::rownames_to_column("cell") %>% tbl_df()  %>% setNames(c("cell","Embryo","Days","Lineage1","Lineage2")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="primitive endoderm", "PrE")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="epiblast", "Epiblast")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="trophectoderm", "TE")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="not applicable", "Prelineage")) %>% mutate(new_lineage=Lineage1)%>% select(cell,new_lineage)

data.ob.umap <- data.ob.umap %>% left_join(cell_lineage_data,by="cell") %>% mutate(EML=ifelse(is.na(new_lineage),EML,new_lineage)) %>% select(-new_lineage)

#'  unified EML annotation
data.ob.umap <- data.ob.umap %>% mutate(EML=ifelse(EML=="EM_NA","Prelineage",EML)) %>% mutate(EML=ifelse(EML=="EPI","Epiblast",EML))%>% mutate(EML=ifelse(EML=="PE","PrE",EML))  %>% mutate(EML=ifelse(pj=="D3post" & EML=="ICM","3D_ICM",EML))

#' rename EML annotation
data.ob.umap <- data.ob.umap %>% mutate(rename_EML=EML) %>% mutate(rename_EML=ifelse(rename_EML%in% c("AdvMes","AxMes","EmMes","NasMes"),"Mesoderm",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML%in% c("CTB","EVT","STB","TE","EarlyTE"),"TE",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML%in% c("MeLC2","MeLC1"),"MeLC",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML%in% c("Tsw-AMLC","AMLC"),"AMLC",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML%in% c("PrE"),"Endoderm",rename_EML))  %>% mutate(rename_EML=ifelse(rename_EML %in% c("EB_ELC","IB_Epi","nicolBla_ELC"),"ELC",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML %in% c("EB_HLC","IB_PE","nicolBla_HLC"),"HLC",rename_EML))%>% mutate(rename_EML=ifelse(rename_EML %in% c("EB_TLC","IB_TE","nicolBla_TLC"),"TLC",rename_EML))  %>% mutate(rename_EML=ifelse(rename_EML%in% c("EB_U10","EB_U5","EB_U6","EB_UE7","EB_UI13","IB_IM1","IB_IM2","IB_IM3","IB_uc","IB_NR","Trans","nicolBla_nc"),"Undef",rename_EML))


#' saving dataset
if (!file.exists(paste0("tmp_data/",TD,"/Pub.D2.withoutAMLC.mnn.rds")) | rewrite) {
  print("save output")
  saveRDS(data.ob ,file=paste0("tmp_data/",TD,"/Pub.D2.withoutAMLC.mnn.rds"))
  saveRDS(data.ob.umap ,file=paste0("tmp_data/",TD,"/Pub.D2.withoutAMLC.umap.cord.rds"))
}


#' visulization check
AP(data.ob)
AP4(data.ob)

