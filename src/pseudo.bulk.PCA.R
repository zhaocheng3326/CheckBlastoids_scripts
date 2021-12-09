#' ---
#' title: "PCA analysis based on Pseudo-bulk dataset"
#' output: 
#'  html_document:
#'    code_folding: hide
#'    
#' ---
#' 
#R3.6

rm(list=ls())
rewrite=FALSE
condaENV <- "/home/chenzh/miniconda3/envs/R3.6"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(DESeq))
suppressMessages(library(Seurat))
suppressMessages(library(VennDiagram))
suppressMessages(library(grDevices))
#source("~/PC/SnkM/SgCell.R")
rename <- dplyr::rename


# working directory
DIR <- "~/My_project/CheckBlastoids"
setwd(DIR)
knitr::opts_knit$set(root.dir=DIR)
source("src/local.quick.fun.R")

options(digits = 4)
options(future.globals.maxSize= 3001289600)
TD="Nov14_2021"
MM="Pub"


rename <- dplyr::rename
select<- dplyr::select
filter <- dplyr::filter

pv.cutoff <- 0.05

#' loading data
load("tmp_data/gene.meta.Rdata",verbose=T)

if (file.exists(paste0("tmp_data/",TD,"/Psd.norm.exp.Rdata"))) {
  load(paste0("tmp_data/",TD,"/Psd.norm.exp.Rdata"),verbose = T)
}else{
  #' loading pseudo counts object
  load(paste0("tmp_data/",TD,"/counts.merge.Rdata"),verbose=T)
  
  #' get universal normalization values (DESeq, median-of-ratios)
  
  #human
  temp.counts <- counts.merge[,colnames(counts.merge)%>% setdiff(c("D3post_TE","EBD2_TLC.mod","nicolBla_TLC.mod"))]#,"JPF2019_Tsw_AMLC"
  temp.sel.gene <- rownames(temp.counts)[rowSums(temp.counts >= 3) >0] ## at least expressed with 3 counts in one dataset
  temp.counts <- temp.counts[temp.sel.gene,]
  temp.cds=newCountDataSet(temp.counts,colnames(temp.counts) ) %>% estimateSizeFactors()
  counts.merge.norm <- counts(temp.cds, normalized=TRUE)
  
  # monkey
  temp.counts <- NHP.counts.merge[ ,colnames(NHP.counts.merge) %>% setdiff("D14_Amnion1")]
  temp.sel.gene <- rownames(temp.counts )[rowSums(temp.counts >= 3) >0] ## at least expressed with 3 counts in one dataset
  temp.counts <- temp.counts[temp.sel.gene ,]
  NHP.cds=newCountDataSet(temp.counts,colnames(temp.counts)) %>% estimateSizeFactors() %>% estimateDispersions(method="blind",sharingMode="fit-only" )
  NHP.counts.merge.norm <- counts(NHP.cds, normalized=TRUE)
  
  # human-monkey cds
  # loading NHP dataset
  load(paste0("tmp_data/",TD,"/NHP.Rdata"),verbose=T)
  temp.sel.gene <- rownames(counts.merge.norm) %>% intersect(rownames(NHP.counts.merge.norm))%>% intersect(GI.list$ov)
  print(length(temp.sel.gene))
  
  temp.counts <- NHP.counts.merge[temp.sel.gene,colnames(NHP.counts.merge) ]#%>% setdiff("D14_Amnion1")
  colnames(temp.counts) <- paste("NHP",colnames(temp.counts),sep="_")
  temp.counts <- temp.counts %>% cbind(counts.merge[temp.sel.gene,])
  
  temp.counts <- temp.counts[,colnames(temp.counts)%>% setdiff(c("EBD2_TLC.mod","D3post_TE","nicolBla_TLC.mod")) ]#""NHP_D14_Amnion1",,JPF2019_Tsw_AMLC"
  IT.merge.norm <- newCountDataSet(temp.counts,colnames(temp.counts) ) %>% estimateSizeFactors() %>%counts(normalized=TRUE)
  
  save(counts.merge.norm,NHP.counts.merge.norm,IT.merge.norm ,file=paste0("tmp_data/",TD,"/Psd.norm.exp.Rdata"))

}

#' loading DE analysis based on single-cell dataset
AmnVSTE.mk <- readRDS(paste0("tmp_data/",TD,"/human.AmnionVsTE.sigMarker.rds"))

pca.out <- list()

#' for human only (based on typical marker)
temp.sel.gene <- AmnVSTE.mk$gene %>% intersect(rownames(counts.merge.norm))
temp.sel.sample <- c("Blakeley_TE","SPH2016_TE","D3post_TE_pre","D3post_TE_post","nBGuo_TE","CS7_Amnion","IBD2_TLC","EBD2_TLC","nBGuo_TLC","JPF2019_AMLC","JPF2019_Tsw_AMLC","nicolBla_TLC","nicolBla_Okae")
temp.norm <- counts.merge.norm[temp.sel.gene ,temp.sel.sample]
set.seed(123)
pca.out$human <- prcomp(log1p(t(temp.norm )), retx=TRUE)

#' human+monkey (based on typical marker)
temp.sel.gene <- AmnVSTE.mk$gene %>% intersect(rownames(IT.merge.norm))
temp.sel.sample <- c("NHP_D10_TE","NHP_D12_TE","NHP_D14_TE","NHP_D14_Amnion1","NHP_D14_Amnion2","Blakeley_TE","SPH2016_TE","D3post_TE_pre","D3post_TE_post","nBGuo_TE","CS7_Amnion","IBD2_TLC","EBD2_TLC","nBGuo_TLC","JPF2019_AMLC","JPF2019_Tsw_AMLC","nicolBla_TLC","nicolBla_Okae")
temp.norm <- IT.merge.norm[temp.sel.gene ,temp.sel.sample]
set.seed(123)
pca.out$HuMonk <-prcomp(log1p(t(temp.norm )), retx=TRUE)

#' save pca object
if (!file.exists(paste0("tmp_data/",TD,"/psd.AmnVsTE.pca.out.rds")) | rewrite) {
  print("save output")
  saveRDS(pca.out,file=paste0("tmp_data/",TD,"/psd.AmnVsTE.pca.out.rds"))
}

#' check the pca results
for (n in names(pca.out)) {
  pca.it <- pca.out[[n]]
 
  pc1.imp <- round(100*summary(pca.it)$importance[2,1],digits=2) 
  pc2.imp <-round(100*summary(pca.it)$importance[2,2],digits=2) 
  
  temp <- pca.it$x[,c("PC1","PC2")] %>% as.data.frame()%>% tibble::rownames_to_column("Sample")
  
  print (
    ggplot(temp,mapping=aes(x=PC1,y=PC2))+geom_point()+geom_text(mapping=aes(label=Sample))+xlab(paste("PC1(",pc1.imp,"% Proportion of Variance)"))+ylab(paste("PC2(",pc2.imp,"% Proportion of Variance)"))
  )
}







