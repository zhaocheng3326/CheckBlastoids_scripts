#' ---
#' title: "normalization"
#' output: 
#'  html_document:
#'    code_folding: hide
#' ---




# ### Loading R library

#R3.6

rm(list=ls())
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(scran))
suppressMessages(library(batchelor))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(pheatmap))

# working directory
DIR <- "~/My_project/JP_project"
setwd(DIR)

#' Loading R functions
source("~/PC/R_code/functions.R")
source("~/PC/SnkM/SgCell.R")
feature.plot.cor <- c("yellow","red","black")


options(digits = 4)
options(future.globals.maxSize= 3001289600)

TD="D2_pub"

#' loading dataset
load("tmp_data/gene.meta.Rdata",verbose=T)
counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
meta.filter<- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds"))


meta.filter<- meta.filter %>% filter(pj %in% c("D3post","CS7","SPH2016","EBD2","IBD2","JPF2019"))

#' get the expressed genes ( expressed in at least 2 datatsets)
expG.set <- list()
for (b in unique(meta.filter$pj  %>% unique() %>% as.vector())) { 
  temp.cell <- meta.filter %>% filter(pj==b) %>% pull(cell)
  expG.set[[b]] <- rownames(counts.filter )[rowSums(counts.filter[,temp.cell] >=1) >=5]
}
sel.expG <-unlist(expG.set) %>% unique() %>% as.vector()

# run scran normalization followed by multiBatchNorm 
sce.ob <- list()

for (b in unique(meta.filter$pj  %>% unique() %>% as.vector())) { 
  print(b)
  temp.M <- meta.filter %>% filter(pj==b) 
  temp.sce <-  SingleCellExperiment(list(counts=as.matrix(counts.filter[sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors()
  sce.ob[[b]] <- temp.sce
}
mBN.sce.ob <- multiBatchNorm(sce.ob$SPH2016,sce.ob$D3post,sce.ob$CS7,sce.ob$JPF2019,sce.ob$IBD2,sce.ob$EBD2)
names(mBN.sce.ob) <- c("SPH2016","D3post","CS7","JPF2019","IBD2","EBD2")   
lognormExp.mBN<- mBN.sce.ob %>% lapply(function(x) {logcounts(x) %>% as.data.frame()  %>% return()}) %>% do.call("bind_cols",.)


#save object
saveRDS(lognormExp.mBN,file=paste0("tmp_data/",TD,"/lognormExp.mBN.rds"))

sessionInfo()








