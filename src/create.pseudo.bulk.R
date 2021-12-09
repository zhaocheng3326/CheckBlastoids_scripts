#' title: crete combined pseudo-bulk dataset
#R4.0

rm(list=ls())
rewrite=FALSE
condaENV <- "/home/chenzh/miniconda3/envs/R4.0"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(scran))
suppressMessages(library(batchelor))
suppressMessages(library(Seurat))
suppressMessages(library(harmony))

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

#' loading data
load("tmp_data/gene.meta.Rdata",verbose=T)
counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
meta.filter<- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% filter(pj %in% c("EBD2","IBD2","JPF2019","D3post","CS7","SPH2016","Blakeley","nBGuo","nicolBla") ) %>% mutate(EML=ifelse(EML %in% c("Ectoderm"),"Amnion",EML)) 

#' loading whole dataset integration results
data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/Pub.D2.umap.cord.rds"))

#' type of cells
cells.list <- list()

cells.list$Blakeley_TE <- meta.filter %>% filter(pj=="Blakeley" & EML=="TE")  %>% pull(cell)
cells.list$SPH2016_TE <- meta.filter %>% filter(pj=="SPH2016" & EML=="TE")  %>% pull(cell)
cells.list$D3post_TE <- meta.filter %>% filter(pj=="D3post" & EML %in% c("CTB","STB","EVT"))  %>% pull(cell)
cells.list$D3post_TE_pre <- meta.filter %>% filter(pj=="D3post" & EML %in% c("CTB","STB","EVT") & devTime %in% c("E6","E7"))  %>% pull(cell)
cells.list$D3post_TE_post <- meta.filter %>% filter(pj=="D3post" & EML %in% c("CTB","STB","EVT") & devTime %in% c("E8","E9","E10","E12","E14"))  %>% pull(cell)
cells.list$nBGuo_TE <- meta.filter %>% filter(pj=="nBGuo" & EML %in% c("EarlyTE","TE"))  %>% pull(cell)
cells.list$CS7_Amnion <- meta.filter %>% filter(pj=="CS7" & EML %in% c("Amnion")) %>% pull(cell)

cells.list$IBD2_TLC <-  meta.filter %>% filter(pj=="IBD2" & EML %in% c("IB_TE")) %>% pull(cell)
cells.list$EBD2_TLC <-  meta.filter %>% filter(pj=="EBD2" & EML %in% c("EB_TLC")) %>% pull(cell)
cells.list$nicolBla_TLC <-  meta.filter %>% filter(pj=="nicolBla" & EML %in% c("nicolBla_TLC")) %>% pull(cell)
cells.list$nicolBla_Okae <-  meta.filter %>% filter(pj=="nicolBla" & EML %in% c("okae_bts5")) %>% pull(cell)
cells.list$EBD2_TLC.mod <-  data.ob.umap %>% filter(pj=="EBD2" & cluster_EML %in% c("TLC")) %>% pull(cell)
cells.list$nicolBla_TLC.mod <-  data.ob.umap %>% filter(pj=="nicolBla" & cluster_EML %in% c("TLC")) %>% pull(cell)
cells.list$nBGuo_TLC <-  meta.filter %>% filter(pj=="nBGuo" & EML %in% c("TLC")) %>% pull(cell)
cells.list$JPF2019_AMLC <-  meta.filter %>% filter(pj=="JPF2019" & EML %in% c("AMLC")) %>% pull(cell)
cells.list$JPF2019_Tsw_AMLC <-  meta.filter %>% filter(pj=="JPF2019" & EML %in% c("Tsw-AMLC")) %>% pull(cell)

lapply(cells.list,length)

counts.merge <- list()
for (n in names(cells.list)) {
  temp.cells <- cells.list[[n]]
  counts.merge[[n]] <- as.data.frame(rowSums(counts.filter[,temp.cells]) ) %>% setNames(n)
}
counts.merge <- counts.merge %>% do.call("cbind",.)


#' loading NHP dataset
load(paste0("tmp_data/",TD,"/NHP.Rdata"),verbose=T)

NHP.cells.list <- list()
NHP.cells.list$D10_TE <- NHP.meta %>% filter(devTime=="E10" & EML=="TE" ) %>% pull(cell) 
NHP.cells.list$D12_TE <- NHP.meta %>% filter(devTime=="E12" & EML=="TE" ) %>% pull(cell) 
NHP.cells.list$D14_TE <- NHP.meta %>% filter(devTime=="E14" & EML=="TE" ) %>% pull(cell) 
NHP.cells.list$D14_Amnion1 <- NHP.meta %>% filter(devTime=="E14" & EML=="Amnion1" ) %>% pull(cell) 
NHP.cells.list$D14_Amnion2 <- NHP.meta %>% filter(devTime=="E14" & EML=="Amnion2" ) %>% pull(cell) 

NHP.counts.merge <- list()

for (n in names(NHP.cells.list)) {
  temp.cells <- NHP.cells.list[[n]]
  NHP.counts.merge[[n]] <- as.data.frame(rowSums(NHP.counts.mod[,temp.cells]) ) %>% setNames(n)
}
NHP.counts.merge <- NHP.counts.merge %>% do.call("cbind",.)

#' saving object
if (!file.exists(paste0("tmp_data/",TD,"/counts.merge.Rdata")) | rewrite) {
  print("save output")
  save(cells.list,counts.merge,NHP.cells.list, NHP.counts.merge,file=paste0("tmp_data/",TD,"/counts.merge.Rdata"))
}

