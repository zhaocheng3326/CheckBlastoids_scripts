#' title: "QC and downsampling"
#' ### Loading R library

#R3.6
rm(list=ls())
rewrite=FALSE

condaENV <- "/home/chenzh/miniconda3/envs/R3.6"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(pheatmap))
suppressMessages(library(Seurat))


suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)

# working directory
DIR <- "~/My_project/CheckBlastoids"
setwd(DIR)
knitr::opts_knit$set(root.dir=DIR)

#' Loading R functions
#source("~/PC/R_code/functions.R")
#source("~/PC/SnkM/SgCell.R")
feature.plot.cor <- c("yellow","red","black")


options(digits = 4)
options(future.globals.maxSize= 3001289600)


#' loading data
TD="Nov14_2021"
load(paste0("tmp_data/",TD,"/all.counts.meta.Rdata"),verbose=T)
load("tmp_data/gene.meta.Rdata",verbose = T)


#' QC 
meta.filter <- meta.all %>% filter(seqType=="smt2" & pj !="Blakeley" & pj !="nBGuo" & pj!="nicolBla") %>% filter(nGene >2000 & mt.perc < 0.125)  %>% bind_rows(meta.all %>% filter(pj  %in% c("Blakeley","nBGuo")) %>% filter(nGene >2000 & mt.perc < 0.2)) %>% bind_rows(readRDS(paste0("tmp_data/",TD,"/meta.AMN.further.rds")))  %>% bind_rows(readRDS(paste0("tmp_data/",TD,"/meta.IBD2.further.rds"))) %>% bind_rows(readRDS(paste0("tmp_data/",TD,"/meta.EBD2.LW60.further.rds"))) %>% bind_rows(readRDS(paste0("tmp_data/",TD,"/meta.nicolBla.further.rds")) %>% filter(! devTime %in% c("24h","60h")))  %>% mutate(EML=gsub("Primitive_Streak","PriS",EML)) %>% mutate(EML=gsub("Emergent_Mesoderm","EmMes",EML)) %>% mutate(EML=gsub("YS_Mesoderm","YsMes",EML)) %>% mutate(EML=gsub("Nascent_Mesoderm","NasMes",EML)) %>% mutate(EML=gsub("Axial_Mesoderm","AxMes",EML)) %>% mutate(EML=gsub("Advanced_Mesoderm","AdvMes",EML)) %>% mutate(EML=gsub("Non-Neural_Ectoderm","Ectoderm",EML)) %>% mutate(EML=gsub("ExE_Mesoderm","ExE_Mes",EML)) %>% filter(! EML %in% c("Hemogenic_Endothelial_Progenitors","Erythroblasts"))



#' loading Downsampled cells
meta.AMN.ds <- readRDS(paste0("tmp_data/",TD,"/meta.AMN.further.ds.rds"))
meta.IBD2.ds <- readRDS(paste0("tmp_data/",TD,"/meta.IBD2.further.ds.rds"))
meta.EBD2.ds <- readRDS(paste0("tmp_data/",TD,"/meta.EBD2.LW60.further.ds.rds"))
meta.nicolBla.ds <- readRDS(paste0("tmp_data/",TD,"/meta.nicolBla.further.ds.rds"))  %>% filter(! devTime %in% c("24h","60h")) 



temp <- meta.filter %>% filter(pj=="SPH2016") %>% filter(cellType=="EM")
meta.sph <- temp



meta.filter.ds <- meta.filter %>% filter(pj %in% c("CS7","D3post","SPH2016","Blakeley","nBGuo")) %>% bind_rows(meta.AMN.ds) %>% bind_rows(meta.IBD2.ds) %>% bind_rows(meta.EBD2.ds) %>% bind_rows(meta.nicolBla.ds) # full sph2016



#' filtering counts
counts.filter <- counts.all[setdiff(rownames(counts.all), mt.gene),meta.filter$cell]
counts.filter.ds <- counts.all[setdiff(rownames(counts.all), mt.gene),meta.filter.ds$cell]

#' saving data

if (!file.exists(paste0("tmp_data/",TD,"/meta.filter.rds")) | rewrite) {
  print("save output")
  saveRDS(counts.filter,file=paste0("tmp_data/",TD,"/counts.filter.rds"))
  saveRDS(counts.filter.ds,file=paste0("tmp_data/",TD,"/counts.filter.ds.rds"))
  saveRDS(meta.filter,file=paste0("tmp_data/",TD,"/meta.filter.rds"))
  saveRDS(meta.filter.ds,file=paste0("tmp_data/",TD,"/meta.filter.ds.rds"))
}





sessionInfo()