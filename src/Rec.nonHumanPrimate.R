#' title: "get non-human primate expression"

#' ### Loading R library
#R4.0
rm(list=ls())
rewrite=FALSE

condaENV <- "/home/chenzh/miniconda3/envs/R4.0"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(pheatmap))
suppressMessages(library(Seurat))
suppressMessages(library(scran))
suppressMessages(library(batchelor))


suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)

# working directory
DIR="/home/chenzh/My_project/CheckBlastoids"
setwd(DIR)
knitr::opts_knit$set(root.dir=DIR)

#' Loading R functions
#source("~/PC/R_code/functions.R")
#source("~/PC/SnkM/SgCell.R")
source("src/local.quick.fun.R")

options(digits = 4)
options(future.globals.maxSize= 3001289600)


#' loading data
TD="Nov14_2021"



if (file.exists(paste0("tmp_data/",TD,"/NHP.Rdata")) & !rewrite) {
  print("Done")
}else{
  hm.counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
  human.gi <- rownames(hm.counts.filter)
  
  Funid_translate <- function(id) {
    annotations <- read.table('~/My_project/CheckBlastoids/data/Primate_10X_Amn/raw_data/mart_export_idtranslate_fun.txt',header = T, sep = ",", stringsAsFactors = F)
    translated <- vector(mode = 'character',length = length(id)) 
    for(i in 1:length(id)){
      if(!is.na(match(id[i], annotations[,1]))){
        if(annotations[match(id[i], annotations[,1]),2] != ''){
          translated[i] <- annotations[match(id[i], annotations[,1]),2]
        }
        else{
          translated[i] <- id[i]
        }
        if(length(grep('LOC', translated[i])!=0)){
          if(annotations[match(id[i], annotations[,1]),3] != ''){
            translated[i] <- annotations[match(id[i], annotations[,1]),3]
          }
        }
      }
      else if (!is.na(match(id[i], annotations[,2]))){
        translated[i] <- annotations[match(id[i], annotations[,2]),1]
      }
      else if (!is.na(match(id[i], annotations[,3]))){
        translated[i] <- annotations[match(id[i], annotations[,3]),1]
      }
      else{
        translated[i] <- id[i]
        warning(paste0('Gene ',id[i], ' not found'))
      }
    }
    return(translated)
  } # thanks Alex
  
  
  NHP.counts <- read.delim("data/Primate_10X_Amn/GSE148683_CM_filtered_uncorrected.tsv",stringsAsFactors = F,head=T,row.names = 1) 
  
  temp.meta <- readRDS("data/Primate_10X_Amn/metadata_full.Rds")
  NHP.meta <- temp.meta$d10  %>% mutate( devTime="E10") %>% bind_rows(temp.meta$d12  %>% mutate( devTime="E12"))%>% bind_rows(temp.meta$d14  %>% mutate( devTime="E14"))%>% mutate(EML=lineage)   %>% mutate(batch=ifelse(grepl("b1",cell),"NHP_b1",ifelse(grepl("b2",cell),"NHP_b2","NHP_b"))) %>% mutate(cell=paste("NHP",cell,sep="_"),seqType="10X",cellType="EM",subCT=NA,mt.perc=NA,pj="NHP") %>% mutate(EML=ifelse(EML=="amnion","Amnion",EML)) %>% mutate(EML=ifelse(EML=="amnion 1","Amnion1",EML))%>% mutate(EML=ifelse(EML=="amnion 2","Amnion2",EML))%>% mutate(EML=ifelse(EML=="endoderm","Endoderm",EML))%>% mutate(EML=ifelse(EML=="epiblast","EPI",EML))%>% mutate(EML=ifelse(EML=="extraembryonic mesenchyme","ExE_Mech",EML))%>% mutate(EML=ifelse(EML=="extraembryonic mesoderm","ExE_Mes",EML))%>% mutate(EML=ifelse(EML=="mesoderm 1","Mes1",EML)) %>% mutate(EML=ifelse(EML=="mesoderm 2","Mes2",EML))%>% mutate(EML=ifelse(EML=="transition","Trans",EML))%>% mutate(EML=ifelse(EML=="trophoblast","TE",EML)) %>% select(-UMAP_1,-UMAP_2,-lineage)
  
  
  
  colnames(NHP.counts) <- paste("NHP",colnames(NHP.counts),sep="_")
  NHP.counts <- NHP.counts[,NHP.meta$cell]
  
  NHP.meta <- data.frame(cell=colnames(NHP.counts),nGene=colSums(NHP.counts > 0),libsize=colSums(NHP.counts)) %>% tbl_df() %>% inner_join(NHP.meta,by="cell")
  
  table(duplicated(NHP.meta$cell))
  
  #' modify the gene name
  NHP.HM.ID.anno <- read.table('~/My_project/CheckBlastoids/data/Primate_10X_Amn/raw_data/mart_export_idtranslate_fun.txt',header = T, sep = ",", stringsAsFactors = F) %>% tbl_df() %>% mutate(human=ifelse(WikiGene.name %in% human.gi,WikiGene.name,Gene.stable.ID)) %>% mutate(human=ifelse(Gene.name %in% human.gi,Gene.name,Gene.stable.ID)) %>% mutate(NHP=Gene.stable.ID) %>% select(NHP,human) %>% unique() %>% filter(NHP %in% rownames(NHP.counts))
  
  table(rownames(NHP.counts) %in% NHP.HM.ID.anno$NHP)
  
  human.gi.need.merge <- NHP.HM.ID.anno %>% filter(human %in%  c(NHP.HM.ID.anno$human[duplicated(NHP.HM.ID.anno$human)])) %>% pull(NHP) %>% unique()
  
  temp.counts.s1 <- NHP.counts[setdiff(rownames(NHP.counts),human.gi.need.merge),]
  rownames(temp.counts.s1 ) <- (NHP.HM.ID.anno  %>% tibble::column_to_rownames("NHP"))[rownames(temp.counts.s1 ),]
  temp.counts.s1 <- temp.counts.s1 %>% tibble::rownames_to_column("human") %>% tbl_df()
  
  temp.counts.s2 <- NHP.counts[intersect(rownames(NHP.counts),human.gi.need.merge),] %>% tibble::rownames_to_column("NHP") %>%  inner_join(NHP.HM.ID.anno,by="NHP") %>% select(-NHP) %>% gather(cell,counts,-human)%>% group_by(human,cell) %>% summarise(counts=sum(counts)) %>% spread(cell,counts) 
  
  NHP.counts.mod <- temp.counts.s1 %>% bind_rows(temp.counts.s2) %>% tibble::column_to_rownames("human")
  
  GI.list <- list()
  GI.list$ov <- rownames(NHP.counts.mod ) %>% intersect(human.gi)
  GI.list$NHP <- rownames(NHP.counts.mod)
  GI.list$human <- human.gi
  GI.list$all <- c(rownames(NHP.counts.mod ),human.gi) %>% unique()
  
  save(GI.list,NHP.counts.mod,NHP.meta,file=paste0("tmp_data/",TD,"/NHP.Rdata"))
  
  
  #' scran normalization for 2 batches
  expG.set <- list()
  for (b in unique(NHP.meta$batch  %>% unique() %>% as.vector())) { 
    temp.cell <- NHP.meta %>% filter(batch==b) %>% pull(cell)
    expG.set[[b]] <- rownames(NHP.counts.mod )[rowSums(NHP.counts.mod[,temp.cell] >=1) >=5]
  }
  sel.expG <-unlist(expG.set) %>% unique() %>% as.vector()
  
  # run scran normalization followed by multiBatchNorm 
  sce.ob <- list()
  for (b in unique(NHP.meta$batch   %>% unique() %>% as.vector())) {  ### 2 batch for the monkey
    print(b)
    temp.M <- NHP.meta %>% filter(batch==b)
    temp.sce <-  SingleCellExperiment(list(counts=as.matrix(NHP.counts.mod[sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors()
    sce.ob[[b]] <- temp.sce
  }
  mBN.sce.ob <- multiBatchNorm(sce.ob$NHP_b1,sce.ob$NHP_b2)
  lognormExp.mBN<- mBN.sce.ob %>% lapply(function(x) {logcounts(x) %>% as.data.frame()  %>% return()}) %>% do.call("bind_cols",.)
  saveRDS(lognormExp.mBN,file=paste0("tmp_data/",TD,"/NHP.scran.norm.rds"))
  
  
  #' create DEG seurat object
  temp.sel.expG <- rownames(lognormExp.mBN)
  temp.M <- NHP.meta
  data.merge <- CreateSeuratObject(NHP.counts.mod[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) 
  data.merge@assays$RNA@data <- as.matrix(lognormExp.mBN[temp.sel.expG,colnames(data.merge)])
  
  #' monkey Amn vs TE
  cells.list <- list()
  cells.list$D14_Amnion2 <-  NHP.meta %>% filter(devTime=="E14" & EML=="Amnion2") %>% pull(cell)
  cells.list$D14_TE <-  NHP.meta %>% filter(devTime=="E14" & EML=="TE") %>% pull(cell)
  cells.list$D10_TE <-  NHP.meta %>% filter(devTime=="E10" & EML=="TE") %>% pull(cell)
  cells.list$D12_TE <-  NHP.meta %>% filter(devTime=="E12" & EML=="TE") %>% pull(cell)
  
  NHP.AMvsTE.out <- list()
  NHP.AMvsTE.out$D14_Amnion2 <- list()
  NHP.AMvsTE.out$D14_Amnion2$D14_TE <- 1
  NHP.AMvsTE.out$D14_Amnion2$D10_TE <- 1
  NHP.AMvsTE.out$D14_Amnion2$D12_TE <- 1
  
  for (n in names(NHP.AMvsTE.out)) {
    for (m in names(NHP.AMvsTE.out[[n]])) {
      NHP.AMvsTE.out[[n]][[m]] <- FunDEG(data.merge,cells.list,n,m) %>% mutate(SampleA=n,SampleB=m) %>% left_join(NHP.HM.ID.anno %>% dplyr::rename(gene=NHP,human.gene=human),by="gene")
    }
  }
  save(NHP.AMvsTE.out,file=paste0("tmp_data/",TD,"/NHP.AmnionVsTE.Rdata"))
}

