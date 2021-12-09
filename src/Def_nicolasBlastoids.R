#' ---
#' title: "define the cell type of nicolas-Blastoids "
#' output: 
#'  html_document:
#'    code_folding: hide
#'    
#' ---


#  R4.0

#' reference:[https://github.com/RivronLab/Human_Blastoid_Kagawa_et_al-)



rm(list=ls())
rewrite=FALSE

condaENV <- "/home/chenzh/miniconda3/envs/R4.0"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(edgeR))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(pheatmap))
suppressMessages(library(Seurat))
suppressMessages(library(readxl))

DIR="/home/chenzh/My_project/CheckBlastoids"
setwd(DIR)
knitr::opts_knit$set(root.dir=DIR)

#source("~/PC/SnkM/SgCell.R")
source("src/local.quick.fun.R")
load("tmp_data/gene.meta.Rdata",verbose=T)

TD="Nov14_2021"

#' QC has been down in merging dataset ( the cells used in original analysis)
savefile=paste0("tmp_data/",TD,"/nicolBla.Rdata")

if (file.exists(savefile)) {
  load(savefile,verbose = T)
}else{
  
  # loading all counts
  load(paste0("tmp_data/",TD,"/all.counts.meta.Rdata"),verbose=T)
  load("tmp_data/gene.meta.Rdata",verbose=T)

  #' get the detailed lineages
  temp.M <- meta.all %>% filter(pj=="nicolBla") %>% inner_join((read.delim("data/GSE177689/GSE177689_series_matrix.trans.txt",stringsAsFactors = F,head=F,sep="\t")) %>% tbl_df() %>% setNames(c("Sample_title","cell","GSM")),by="cell")
  temp.cell <-temp.M$cell
  temp.sel.expG <- rownames(counts.all)[rowSums(counts.all[,temp.cell] >=1) >=3] %>% setdiff(mt.gene)
  
  data.temp <- CreateSeuratObject(counts.all [temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)  %>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>% ScaleData(verbose=F)%>% RunPCA(verbose=F,npcs=30) %>% RunUMAP(dims=1:20,verbose=F) %>% FindNeighbors( dims = 1:20,verbose = FALSE)  %>% FindClusters(reso=0.5,verbose=F)
  save(data.temp,temp.M,file=savefile)
}

cowplot::plot_grid(
  DimPlot(data.temp,group.by="devTime")+NoAxes(),
 DimPlot(data.temp,label=T)+NoAxes()+NoLegend(),
 nrow=2,ncol=2
)
DimPlot(data.temp,group.by="Sample_title")+NoAxes()
FeaturePlot(data.temp,c("GATA2","GATA3","POU5F1","KLF17","GATA4","SOX17"))

#' check the AdMes/ExE Mes marker genes expression
FeaturePlot(data.temp,c("LIX1","PMP22","COL3A1","ANXA1","COL1A1","VIM"))


data.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(c(cell,SID:mt.perc)) %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% mutate(seurat_cluster=paste0('C',as.vector(Idents(data.temp)))) %>% mutate(EML=ifelse(EML %in% c("naive_H9","primed_H9","okae_bts5"),EML,"nicolBla_nc")) %>% mutate(EML=ifelse(EML=="nicolBla_nc" & seurat_cluster %in% c("C1","C4"), "nicolBla_TLC",EML ))%>% mutate(EML=ifelse(EML=="nicolBla_nc" & seurat_cluster %in% c("C6"), "nicolBla_HLC",EML )) %>% mutate(EML=ifelse(EML=="nicolBla_nc" & seurat_cluster %in% c("C0","C2","C3","C5","C9"), "nicolBla_ELC",EML ))

data.mod <- data.temp
Idents(data.mod) <- as.factor((data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.mod@meta.data),"EML"])
DimPlot(data.mod,label=T)+NoAxes()+NoLegend()

nicolBla.meta.filter <- data.ob.umap %>% select(cell,SID:mt.perc)
nicolBla.meta.filter.ds <- nicolBla.meta.filter  %>% split(nicolBla.meta.filter$EML) %>% lapply(function(x){return(FunMaSF(x,100))}) %>% do.call("bind_rows",.)

#' save object 
if (!file.exists(paste0("tmp_data/",TD,"/meta.nicolBla.further.rds")) | rewrite) {
  saveRDS(nicolBla.meta.filter,file=paste0("tmp_data/",TD,"/meta.nicolBla.further.rds"))
  saveRDS(nicolBla.meta.filter.ds,file=paste0("tmp_data/",TD,"/meta.nicolBla.further.ds.rds"))
}
