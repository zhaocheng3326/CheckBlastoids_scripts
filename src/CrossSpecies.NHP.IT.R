#' ---
#' title: "human dataset integrated with NHP dataset"
#' output: 
#'  html_document:
#'    code_folding: hide
#'    
#' ---

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

nGene <- 2000;pc <- 25

rename <- dplyr::rename
select<- dplyr::select
filter <- dplyr::filter

#' loading data
load("tmp_data/gene.meta.Rdata",verbose=T)
counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))

meta.filter<- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% filter(pj %in% c("EBD2","IBD2","JPF2019","D3post","CS7","SPH2016","Blakeley","nBGuo","nicolBla") ) %>% mutate(EML=ifelse(EML %in% c("Ectoderm"),"Amnion",EML)) %>% filter(pj %in% c("CS7","EBD2","IBD2","nBGuo","nicolBla"))

#' loading NHP dataset
load(paste0("tmp_data/",TD,"/NHP.Rdata"),verbose=T)


#' check the overlapped genes
GI.list$ov %in% rownames(counts.filter) %>% table()
GI.list$ov %in% rownames(NHP.counts.mod) %>% table()
# 14846 are overlapped

#' combine human+monkey, for NHP dataset only choose the D14
meta <- meta.filter %>% bind_rows(NHP.meta %>% mutate(seqType="10X", cellType="EM") %>% filter(devTime=="E14"))
counts <- (counts.filter[GI.list$ov,] %>% cbind(NHP.counts.mod[GI.list$ov,]))[,meta$cell]

#' rename EML annotation
meta <- meta %>% mutate(rename_EML=EML) %>% mutate(rename_EML=ifelse(rename_EML %in% "Trans" & pj=="NHP","EPI-AM-trans",rename_EML))%>% mutate(rename_EML=ifelse(rename_EML%in% c("AdvMes","AxMes","EmMes","NasMes","Mes1","Mes2"),"Mesoderm",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML%in% c("CTB","EVT","STB","TE","EarlyTE"),"TE",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML%in% c("MeLC2","MeLC1"),"MeLC",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML%in% c("Tsw-AMLC","AMLC"),"AMLC",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML%in% c("PrE"),"Endoderm",rename_EML))  %>% mutate(rename_EML=ifelse(rename_EML %in% c("EB_ELC","IB_Epi","nicolBla_ELC"),"ELC",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML %in% c("EB_HLC","IB_PE","nicolBla_HLC"),"HLC",rename_EML))%>% mutate(rename_EML=ifelse(rename_EML %in% c("EB_TLC","IB_TE","nicolBla_TLC"),"TLC",rename_EML))  %>% mutate(rename_EML=ifelse(rename_EML%in% c("EB_U10","EB_U5","EB_U6","EB_UE7","EB_UI13","IB_IM1","IB_IM2","IB_IM3","IB_uc","IB_NR","Trans","nicolBla_nc"),"Undef",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML %in% "Amnion1","E-Amnion",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML %in% "Amnion2","Amnion",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML %in% "EPI","Epiblast",rename_EML))



expG.set <- list()
for (b in unique(meta$pj  %>% unique() %>% as.vector())) { 
  temp.cell <- meta %>% filter(pj==b) %>% pull(cell)
  expG.set[[b]] <- rownames(counts )[rowSums(counts[,temp.cell] >=1) >=5]
}
sel.expG <-unlist(expG.set) %>% unique() %>% as.vector()

if (file.exists(file=paste0("tmp_data/",TD,"/NHP.IT.umap.rds")) & !rewrite) {
  data.umap.list <- readRDS(paste0("tmp_data/",TD,"/NHP.IT.umap.rds"))
}else{
  
  data.umap.list <- list()
  for (s in c("CS7","EBD2","IBD2","nBGuo","nicolBla")) {
    print(paste0("now it is ",s ))
    
    # run scran normalization followed by multiBatchNorm  ( it doesn't effect the integration results can be removed)
    # sce.ob <- list()
    # for (b in c("NHP",s)) {
    #   print(b)
    #   temp.M <- meta %>% filter(pj==b)
    #   temp.sce <-  SingleCellExperiment(list(counts=as.matrix(counts[sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors()
    #   sce.ob[[b]] <- temp.sce
    # }
    # mBN.sce.ob <- multiBatchNorm(sce.ob[[1]],sce.ob[[2]])
    # lognormExp.mBN<- mBN.sce.ob %>% lapply(function(x) {logcounts(x) %>% as.data.frame()  %>% return()}) %>% do.call("bind_cols",.)

    #' prepare for integration
    temp.M <-  meta %>% filter(pj=="NHP") %>% bind_rows(meta %>% filter(pj==s))
    temp.sel.cell <- temp.M$cell
    data.merge <- CreateSeuratObject(counts[sel.expG,temp.sel.cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) 
    #data.merge@assays$RNA@data <- as.matrix(lognormExp.mBN[rownames(lognormExp.mBN),colnames(data.merge)]) 
    data.spt <- SplitObject(data.merge, split.by = "pj")%>% lapply(function(x){x=FindVariableFeatures(x,verbose=F,nfeatures=2000)})
    
    #' fastmnn  
    data.spt <- data.spt;set.seed(123)
    data.temp <- SeuratWrappers::RunFastMNN(data.spt,verbose=F,features=nGene) %>% RunUMAP( reduction = "mnn", dims = 1:pc,verbose=F) 
    saveRDS(data.temp,paste0("tmp_data/",TD,"/NHP.",s,".mnn.rds"))
    data.umap.list[[paste0("NHP.",s,".mnn")]]<- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(c(cell,SID:rename_EML)) %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") 
    
    #' seurat CCA
    sel.features <- SelectIntegrationFeatures(object.list = data.spt, nfeatures = nGene,verbose=F)
    data.temp <- FindIntegrationAnchors(object.list = data.spt, dims = 1:pc,anchor.features=sel.features,verbose=F) %>% IntegrateData(dims = 1:pc,verbose=F) %>% ScaleData(verbose = FALSE) %>% RunPCA(npcs = pc, verbose = FALSE) %>%  RunUMAP( reduction = "pca", dims = 1:pc,verbose=F) ### 
    saveRDS(data.temp,paste0("tmp_data/",TD,"/NHP.",s,".cca.rds"))
    data.umap.list[[paste0("NHP.",s,".cca")]]<- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(c(cell,SID:rename_EML)) %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") 
    
    #' harmony
    data.temp <- data.merge  %>% FindVariableFeatures(selection.method = "vst", nfeatures = nGene, verbose = FALSE) %>% ScaleData(verbose = FALSE) %>% RunPCA( npcs = pc, verbose = FALSE) %>%  RunHarmony("pj", plot_convergence = TRUE, verbose = FALSE) %>%  RunUMAP(reduction = "harmony", dims = 1:pc, verbose = FALSE)
    saveRDS(data.temp,paste0("tmp_data/",TD,"/NHP.",s,".harmony.rds"))
    data.umap.list[[paste0("NHP.",s,".harmony")]]<- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(c(cell,SID:rename_EML)) %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") 
    
  }
 
  saveRDS(data.umap.list,file=paste0("tmp_data/",TD,"/NHP.IT.umap.rds"))
}

#' the integration results
#+ fig.width=15, fig.height=15
for (m in c("mnn","cca","harmony")) {
  temp.plot <- list()
  for (s in c("CS7","EBD2","IBD2","nBGuo","nicolBla")) {
    data.temp <- readRDS(paste0("tmp_data/",TD,"/NHP.",s,".",m,".rds"))
    temp.plot[[paste(s,"all")]] <- DimPlot(data.temp,group.by="EML",label=T)+NoAxes()+NoLegend()+ggtitle(paste0("NHP + ", s))+theme(plot.title = element_text(hjust=0.5,face="bold"))
    temp.plot[[paste(s,"split")]] <- DimPlot(data.temp,group.by="EML",split.by="pj",label=T)+NoAxes()+NoLegend()+ggtitle("")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  }
  print(
    cowplot::plot_grid(plotlist=temp.plot,nrow=5,rel_widths=c(1:2))
  )
}






  