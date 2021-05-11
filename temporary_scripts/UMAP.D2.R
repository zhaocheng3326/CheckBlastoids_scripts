#' ---
#' title: "UMAP projection"
#' output: 
#'  html_document:
#'    code_folding: hide
#' ---




# ### Loading R library

#R3.6
rm(list=ls())

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratWrappers))


# working directory
DIR <- "~/My_project/JP_project"
setwd(DIR)

#' Loading R functions
source("~/PC/R_code/functions.R")
source("~/PC/SnkM/SgCell.R")
source("loca.quick.fun.R")
feature.plot.cor <- c("yellow","red","black")

suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)


options(digits = 4)
options(future.globals.maxSize= 3001289600)
TD="D2_pub"
nGene <- 2000
pc=25

rename <- dplyr::rename
select<- dplyr::select
filter <- dplyr::filter

#' loading data
load("tmp_data/gene.meta.Rdata",verbose=T)
counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/lognormExp.mBN.rds"))

meta.filter<- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% filter(pj %in% c("EBD2","IBD2","JPF2019","D3post","CS7","SPH2016") ) %>% mutate(EML=ifelse(EML %in% c("Ectoderm"),"Amnion",EML))


#temp.M <- meta.filter %>% filter(pj %in% c("IBD2","D3post","CS7","SPH2016","JPF2019"))
temp.M <- meta.filter 
temp.sel.expG <-rownames(lognormExp.mBN )

data.merge <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) 
data.merge@assays$RNA@data <- as.matrix(lognormExp.mBN[temp.sel.expG,colnames(data.merge)])
data.spt <- SplitObject(data.merge, split.by = "pj")%>% lapply(function(x){x=FindVariableFeatures(x,verbose=F,nfeatures=2000)})

#' release memory
rm(lognormExp.mBN,data.merge,counts.filter)
data.spt <- data.spt[c("JPF2019","IBD2","EBD2","D3post","CS7","SPH2016")]


#' #### mnn correction
set.seed(123)
data.ob <- RunFastMNN(data.spt,verbose=F,features=nGene) %>% RunUMAP( reduction = "mnn", dims = 1:pc,verbose=F) %>% FindNeighbors( reduction = "mnn", dims = 1:pc)
rm(data.spt)

data.temp <- data.ob %>% FindClusters(reso=0.4,verbose=F)

#' #### check the UMAP
#+ fig.width=9,fig.height=9
temp.plot <- list()
temp.plot$raw <- DimPlot(data.temp,label=T)+theme_void()+NoLegend()
temp.plot$EML <- DimPlot(data.temp,group.by="EML",label=T)+theme_void()+NoLegend()
print(cowplot::plot_grid(plotlist=temp.plot,nrow=2,ncol=2))


data.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(c(cell,SID:mt.perc)) %>% mutate(seurat_clusters=paste0("C",as.vector(Idents(data.temp))))  %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$mnn@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:mnn_10),by="cell")

#' using old sph2016 annotation
cell_lineage_data <- read.delim("~/My_project/ScRNA_analysis/data/Published/E-MTAB-3929.short.txt",stringsAsFactors = F,row.names = 1)%>% tibble::rownames_to_column("cell") %>% tbl_df()  %>% setNames(c("cell","Embryo","Days","Lineage1","Lineage2")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="primitive endoderm", "PrE")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="epiblast", "Epiblast")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="trophectoderm", "TE")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="not applicable", "Prelineage")) %>% mutate(new_lineage=Lineage1)%>% select(cell,new_lineage)

data.ob.umap <- data.ob.umap %>% left_join(cell_lineage_data,by="cell") %>% mutate(EML=ifelse(is.na(new_lineage),EML,new_lineage)) %>% select(-new_lineage)

#' rename EML annotation
data.ob.umap <- data.ob.umap %>% mutate(EML=ifelse(EML=="EM_NA","Prelineage",EML)) %>% mutate(EML=ifelse(EML=="EPI","Epiblast",EML))%>% mutate(EML=ifelse(EML=="PE","PrE",EML))  %>% mutate(EML=ifelse(pj=="D3post" & EML=="ICM","3D_ICM",EML))

#' create modified annotation
data.ob.umap <- data.ob.umap %>% mutate(new_EML=EML) %>% mutate(new_EML=ifelse(new_EML%in% c("AdvMes","AxMes","EmMes","NasMes"),"Mesoderm",new_EML)) %>% mutate(new_EML=ifelse(new_EML%in% c("CTB","EVT","STB","TE"),"TE",new_EML)) %>% mutate(new_EML=ifelse(new_EML%in% c("MeLC2","MeLC1"),"MeLC",new_EML)) %>% mutate(new_EML=ifelse(new_EML%in% c("Tsw-AMLC","AMLC"),"AMLC",new_EML)) %>% mutate(new_EML=ifelse(new_EML%in% c("PrE"),"Endoderm",new_EML))  %>% mutate(new_EML=ifelse(new_EML=="EB_ELC","ELC",new_EML)) %>% mutate(new_EML=ifelse(new_EML=="EB_HLC","HLC",new_EML))%>% mutate(new_EML=ifelse(new_EML=="EB_TLC","TLC",new_EML)) %>% mutate(new_EML=ifelse(new_EML=="IB_Epi","ELC",new_EML))%>% mutate(new_EML=ifelse(new_EML=="IB_PE","HLC",new_EML))%>% mutate(new_EML=ifelse(new_EML=="IB_TE","TLC",new_EML))  %>% mutate(new_EML=ifelse(new_EML%in% c("EB_U10","EB_U5","EB_U6","EB_UE7","EB_UI13","IB_IM1","IB_IM2","IB_IM3","IB_uc","IB_NR"),"Undef",new_EML))


#' create cluster annotation
data.ob.umap <- data.ob.umap %>% filter(pj %in% c("EBD2","IBD2")) %>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C4","C8"),"ELC","Undef")) %>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C5"),"HLC",cluster_EML))%>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C3"),"TLC",cluster_EML))%>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C0","C2","C11"),"MeLC",cluster_EML)) %>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C1","C6"),"AMLC",cluster_EML)) %>% bind_rows(data.ob.umap %>% filter(!pj %in% c("EBD2","IBD2")) %>% mutate(cluster_EML=new_EML))
#' saving dataset


saveRDS(data.ob ,file=paste0("tmp_data/",TD,"/Pub.D2.mnn.rds"))
saveRDS(data.ob.umap ,file=paste0("tmp_data/",TD,"/Pub.D2.umap.cord.rds"))

sessionInfo()
