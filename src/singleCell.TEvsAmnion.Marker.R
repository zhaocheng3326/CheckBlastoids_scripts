#' title: "find the markers to distinct TE and Amnion"

#R4.0
rm(list=ls())

condaENV <- "/home/chenzh/miniconda3/envs/R4.0"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Seurat))



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

rename <- dplyr::rename
select<- dplyr::select
filter <- dplyr::filter


#' cutoffs
power.cutoff <- 0.5
log2FC.cutoff <- 0.5
sel.sets <-  c("SPH2016_TE","D3post_TE_pre","D3post_TE_post") # c("Blakeley_TE","SPH2016_TE","D3post_TE_pre","D3post_TE_post","nBGuo_TE")
nC.cutoff <- length(sel.sets) 
pct.min.cutoff <- c(0.25)
pct.max.cutoff <- c(0.5)

#' loading Gene annotation
Gene.anno <- read.delim("~/Genome/Human/RefSeq/Homo_sapiens.GRCh38.gene.description",sep="\t",stringsAsFactors=F,header = F)  %>% tbl_df() %>% rename(gene=V1,annotaton=V2)


# number of conserved markers
if (file.exists(paste0("tmp_data/",TD,"/human.AmnionVsTE.Rdata")))  {
  load(paste0("tmp_data/",TD,"/human.AmnionVsTE.Rdata"),verbose=T)
}else{
  #' loading data
  load("tmp_data/gene.meta.Rdata",verbose=T)
  counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
  lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/",MM,".lognormExp.mBN.rds"))
  #' whole dataset integration results
  data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/Pub.D2.umap.cord.rds")) 
  
  meta.filter<- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% filter(pj %in% c("EBD2","IBD2","JPF2019","D3post","CS7","SPH2016","Blakeley","nBGuo","nicolBla") ) %>% mutate(EML=ifelse(EML %in% c("Ectoderm"),"Amnion",EML))
  
  #temp.M <- meta.filter %>% filter(pj %in% c("IBD2","D3post","CS7","SPH2016","JPF2019"))
  temp.M <- meta.filter 
  temp.sel.expG <-rownames(lognormExp.mBN )
  
  data.merge <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) 
  data.merge@assays$RNA@data <- as.matrix(lognormExp.mBN[temp.sel.expG,colnames(data.merge)])
  
  
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
  cells.list$nBGuo_TLC <-  meta.filter %>% filter(pj=="nBGuo" & EML %in% c("TLC")) %>% pull(cell)
  cells.list$nicolBla_TLC <-  meta.filter %>% filter(pj=="nicolBla" & EML %in% c("nicolBla_TLC")) %>% pull(cell)
  
  cells.list$EBD2_TLC.mod <-  data.ob.umap %>% filter(pj=="EBD2" & cluster_EML %in% c("TLC")) %>% pull(cell)
  cells.list$nicolBla_TLC.mod <-  data.ob.umap %>% filter(pj=="nicolBla" & cluster_EML %in% c("TLC")) %>% pull(cell)
  
  
  lapply(cells.list,length)

  hm.DE.list <- list()
  hm.DE.list$CS7_Amnion <- list()
  hm.DE.list$CS7_Amnion$Blakeley_TE <- 1
  hm.DE.list$CS7_Amnion$SPH2016_TE <- 1
  hm.DE.list$CS7_Amnion$D3post_TE <- 1
  hm.DE.list$CS7_Amnion$D3post_TE_pre <- 1
  hm.DE.list$CS7_Amnion$D3post_TE_post <- 1
  hm.DE.list$CS7_Amnion$nBGuo_TE <- 1
  
  # hm.DE.list$IBD2_TLC$EBD2_TLC <- 1
  # hm.DE.list$IBD2_TLC$EBD2_TLC.mod <- 1
  # hm.DE.list$IBD2_TLC$nBGuo_TLC <- 1
  # hm.DE.list$EBD2_TLC$nBGuo_TLC <- 1
  # hm.DE.list$EBD2_TLC.mod$nBGuo_TLC <- 1

  hm.AMvsTE.out <- list()
  for (n in names(hm.DE.list)) {
    for (m in names(hm.DE.list[[n]])) {
      print(paste(m,n))
      hm.AMvsTE.out[[n]][[m]] <- FunDEG(data.merge,cells.list,n,m,psd.counts=1) %>% mutate(SampleA=n,SampleB=m)
    }
  }
  
  save(hm.AMvsTE.out,file=paste0("tmp_data/",TD,"/human.AmnionVsTE.Rdata"))
}


# loading NHP Amnion vs TE resutls
#load(paste0("tmp_data/",TD,"/NHP.AmnionVsTE.Rdata"),verbose=T)

# loading psed-bulk DEG results
#load(paste0("tmp_data/",TD,"/merge.DEseq.Rdata"),verbose = T)

AmnVSTE.mk <-  hm.AMvsTE.out$CS7_Amnion %>% do.call("bind_rows",.) %>% filter(SampleB %in% c(sel.sets))%>% mutate(Type=ifelse(power > power.cutoff & avg_log2FC > log2FC.cutoff & pct.1 > pct.max.cutoff & pct.2 < pct.min.cutoff, "UpDEG",ifelse(power > power.cutoff & avg_log2FC < -log2FC.cutoff & pct.1 < pct.min.cutoff & pct.2 > pct.max.cutoff, "DownDEG","NotDEG"))) %>% filter(Type !="NotDEG")%>% group_by(gene,Type) %>% summarise(nC=n()) %>% filter(nC >=nC.cutoff) %>% mutate(compair="em") %>% select(-nC) %>% ungroup()%>% select(-compair) %>% mutate(Type=gsub("DownDEG","TE",Type)) %>% mutate(Type=gsub("UpDEG","Amnion",Type)) %>% arrange(Type)%>% left_join(Gene.anno,by="gene")

if (!file.exists(paste0("tmp_data/",TD,"/human.AmnionVsTE.sigMarker.rds")) | rewrite) {
  print("save output")
  saveRDS(AmnVSTE.mk,file=paste0("tmp_data/",TD,"/human.AmnionVsTE.sigMarker.rds"))
}

#' ### some check
#' froom keyGene paper https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4472038/
keyGene.set <- c("TNMD","GABRP","ITGB6","HOXD9","VTCN1","KRT24","MUC16")
keyGene.set[keyGene.set %in% AmnVSTE.mk$gene]
#' from Io et al., Cell Stem Cell 2021,
c("TNC", "VIM", "VTCN1" ) %in%  AmnVSTE.mk$gene
hm.AMvsTE.out$CS7_Amnion %>% do.call("bind_rows",.) %>% filter(SampleB %in% c(sel.sets)) %>% filter(gene %in% c("TNC", "VIM", "VTCN1" )) %>% arrange(gene)

CSC.gene <- c("FABP3","GCM1","S100A14","DNMT3L","DPPA3","HAVCR1","CGA","GTSF1","ISL1","HEY1","GABRP","MIF","PLA2G2A","IGFBP7","IGFBP3")
CSC.gene[CSC.gene %in% AmnVSTE.mk$gene]

hm.AMvsTE.out$CS7_Amnion %>% do.call("bind_rows",.) %>% filter(SampleB %in% c(sel.sets)) %>% filter(gene %in% CSC.gene[!CSC.gene %in% AmnVSTE.mk$gene]) %>% arrange(gene)
#AmnVSTE.mk=readRDS(paste0("tmp_data/",TD,"/human.AmnionVsTE.sigMarker.rds"))
#AmnVSTE.mk %>% write.table("~/AmnVSTE.mk.tsv",quote=F,col.names = T,row.names = F,sep="\t")

########################### not used
#key.marker <- c("ISL1","GABRP","IGFBP5","HEY1","GCM1","GATA2","GATA3","KRT19")
# key.marker2 <- c("GATA2","GATA3","CDX2","TBX3","KRT7","KRT19","GCM1","ITGA6","EGFR","TP63","NR2F2","VGLL1","TFAP2A","TNC", "CDH10","PKDCC", "IGFBP3","PGFBP5","GABRP", "VIM","HEY1", "ISL1", "POSTN", "WNT6")

# check the cutoff
#hm.AMvsTE.out$CS7_Amnion %>% do.call("bind_rows",.)  %>% filter(SampleB!="D3post_TE") %>% filter(gene %in% key.marker) %>% pull(power) %>% min()
#hm.AMvsTE.out$CS7_Amnion %>% do.call("bind_rows",.)  %>% filter(SampleB!="D3post_TE") %>% filter(gene %in% key.marker) %>% pull(avg_log2FC) %>% abs() %>%  min()

#hm.AMvsTE.out$CS7_Amnion %>% do.call("bind_rows",.)  %>% filter(SampleB!="D3post_TE") %>% filter(gene %in% c("ISL1","GABRP","IGFBP5","HEY1"))

#hm.AMvsTE.out$CS7_Amnion %>% do.call("bind_rows",.)  %>% filter(SampleB!="D3post_TE") %>% filter(gene %in% c("GCM1","GATA2","GATA3","KRT19"))

