#' ---
#' title: "Integration of all human dataset"
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
suppressMessages(library(scran))
suppressMessages(library(batchelor))

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

meta.filter<- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% filter(pj %in% c("EBD2","IBD2","JPF2019","D3post","CS7","SPH2016","Blakeley","nBGuo","nicolBla") ) %>% mutate(EML=ifelse(EML %in% c("Ectoderm"),"Amnion",EML))

temp.M <- meta.filter 
temp.sel.expG <-rownames(lognormExp.mBN )

data.merge <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) 
data.merge@assays$RNA@data <- as.matrix(lognormExp.mBN[temp.sel.expG,colnames(data.merge)])
data.spt <- SplitObject(data.merge, split.by = "pj")%>% lapply(function(x){x=FindVariableFeatures(x,verbose=F,nfeatures=2000)})
#' release memory
rm(lognormExp.mBN,data.merge,counts.filter)

data.spt <- data.spt[c("JPF2019","IBD2","EBD2","nicolBla","D3post","CS7","SPH2016","nBGuo")]
set.seed(123)

 # pdf("temp1.pdf")
 # for (nGene in c(2000,2500,3000,3500,4000,4500,5000)) {
 #   data.ob <- RunFastMNN(data.spt,verbose=F,features=nGene) %>% RunUMAP( reduction = "mnn", dims = 1:pc,verbose=F)
 #    print(nGene)
 #   AP(data.ob)
 #   APPJ(data.ob,"IBD2")
 # }
# dev.off()
data.ob <- RunFastMNN(data.spt,verbose=F,features=nGene) %>% RunUMAP( reduction = "mnn", dims = 1:pc,verbose=F) %>% FindNeighbors( reduction = "mnn", dims = 1:pc)



data.temp <- data.ob %>% FindClusters(reso=0.6,verbose=F)
data.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(c(cell,SID:mt.perc)) %>% mutate(seurat_clusters=paste0("C",as.vector(Idents(data.temp))))%>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$mnn@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:mnn_10),by="cell")


#' update sph2016 annotation
cell_lineage_data <- read.delim("~/My_project/ScRNA_analysis/data/Published/E-MTAB-3929.short.txt",stringsAsFactors = F,row.names = 1)%>% tibble::rownames_to_column("cell") %>% tbl_df()  %>% setNames(c("cell","Embryo","Days","Lineage1","Lineage2")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="primitive endoderm", "PrE")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="epiblast", "Epiblast")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="trophectoderm", "TE")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="not applicable", "Prelineage")) %>% mutate(new_lineage=Lineage1)%>% select(cell,new_lineage)

data.ob.umap <- data.ob.umap %>% left_join(cell_lineage_data,by="cell") %>% mutate(EML=ifelse(is.na(new_lineage),EML,new_lineage)) %>% select(-new_lineage)

#'  unified EML annotation
data.ob.umap <- data.ob.umap %>% mutate(EML=ifelse(EML=="EM_NA","Prelineage",EML)) %>% mutate(EML=ifelse(EML=="EPI","Epiblast",EML))%>% mutate(EML=ifelse(EML=="PE","PrE",EML))  %>% mutate(EML=ifelse(pj=="D3post" & EML=="ICM","3D_ICM",EML))

#' rename EML annotation
data.ob.umap <- data.ob.umap %>% mutate(rename_EML=EML) %>% mutate(rename_EML=ifelse(rename_EML%in% c("AdvMes","AxMes","EmMes","NasMes"),"Mesoderm",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML%in% c("CTB","EVT","STB","TE","EarlyTE"),"TE",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML%in% c("MeLC2","MeLC1"),"MeLC",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML%in% c("Tsw-AMLC","AMLC"),"AMLC",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML%in% c("PrE"),"Endoderm",rename_EML))  %>% mutate(rename_EML=ifelse(rename_EML %in% c("EB_ELC","IB_Epi","nicolBla_ELC"),"ELC",rename_EML)) %>% mutate(rename_EML=ifelse(rename_EML %in% c("EB_HLC","IB_PE","nicolBla_HLC"),"HLC",rename_EML))%>% mutate(rename_EML=ifelse(rename_EML %in% c("EB_TLC","IB_TE","nicolBla_TLC"),"TLC",rename_EML))  %>% mutate(rename_EML=ifelse(rename_EML%in% c("EB_U10","EB_U5","EB_U6","EB_UE7","EB_UI13","IB_IM1","IB_IM2","IB_IM3","IB_uc","IB_NR","Trans","nicolBla_nc"),"Undef",rename_EML))
#' create cluster EML annotation
data.ob.umap <- data.ob.umap %>% filter(pj %in% c("EBD2","IBD2","nBGuo")) %>% filter(cellType!="EM")  %>% bind_rows(data.ob.umap %>% filter(pj %in% c("nicolBla")) %>% filter(!EML %in% c("okae_bts5","naive_H9","primed_H9")))%>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C0","C7","C9"),"ELC","Undef")) %>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C5"),"HLC",cluster_EML))%>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C3","C12","C14"),"TLC",cluster_EML))%>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C4","C6"),"MeLC",cluster_EML)) %>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C1","C8","C10"),"AMLC",cluster_EML)) %>% bind_rows(data.ob.umap %>% filter(!pj %in% c("EBD2","IBD2","nicolBla")) %>% filter(!(pj=="nBGuo" & cellType!="EM")) %>% mutate(cluster_EML=rename_EML)) %>% bind_rows(data.ob.umap %>% filter(pj %in% "nicolBla") %>% filter(EML %in% c("okae_bts5","naive_H9","primed_H9")) %>% mutate(cluster_EML=rename_EML))

#' saving dataset
if (!file.exists(paste0("tmp_data/",TD,"/Pub.D2.mnn.rds")) | rewrite) {
  print("save output")
  saveRDS(data.ob ,file=paste0("tmp_data/",TD,"/Pub.D2.mnn.rds"))
  saveRDS(data.ob.umap ,file=paste0("tmp_data/",TD,"/Pub.D2.umap.cord.rds"))
}  



#' visualization check
temp.plot <- list()
temp.plot$raw <- DimPlot(data.temp,label=T)+theme_void()+NoLegend()
temp.plot$EML <- DimPlot(data.temp,group.by="EML",label=T)+theme_void()+NoLegend()
print(cowplot::plot_grid(plotlist=temp.plot,nrow=2,ncol=2))

AP4(data.ob)
APnGuo(data.ob)

temp.plot <- list()
for ( n in unique(data.ob.umap$rename_EML)) {
  if ( ! n %in% c("3D_ICM","PSA-EPI","Undef")) {
    temp <- data.ob.umap %>% filter(rename_EML==n) %>% group_by(pj, seurat_clusters,rename_EML) %>% summarise(nCell=n_distinct(cell))%>% group_by(pj) %>% mutate(Perc=nCell/sum(nCell))
    temp.plot[[n]] <-  temp %>% ggplot()+geom_bar(mapping=aes(x=seurat_clusters,fill=pj,y=Perc),stat = "identity")+theme_classic()+ggtitle(n)+FunTitle()
    print(temp.plot[[n]] )
  }else if (n %in% "Undef"){
    temp <- data.ob.umap %>% filter(rename_EML==n) %>% group_by(pj, seurat_clusters,EML) %>% summarise(nCell=n_distinct(cell))%>% group_by(pj) %>% mutate(Perc=nCell/sum(nCell))
    temp.plot[[n]] <-  temp %>% ggplot()+geom_bar(mapping=aes(x=seurat_clusters,fill=EML,y=Perc),stat = "identity")+theme_classic()+ggtitle(n)+FunTitle() 
  }
}

#' check the cluster annotation
temp.plot <- list()
temp.plot$raw <- DimPlot(data.temp,label=T)+theme_void()+NoLegend()
temp.plot$EML <- DimPlot(data.temp,group.by="EML",label=T)+theme_void()+NoLegend()
temp.plot$EM <- DimPlot(subset(data.ob,cells=(data.ob.umap %>% filter(cellType=="EM") %>% filter(pj!="JPF2019") %>% pull(cell))),group.by="EML",label=T)+theme_void()+NoLegend()
temp.plot$ELC <- DimPlot(data.ob,cells.highlight=(data.ob.umap %>% filter(cluster_EML=="ELC") %>% pull(cell)) )+theme_void()+NoLegend()+ggtitle("ELC")+theme(plot.title = element_text(hjust=0.5,face="bold"))
temp.plot$HLC <- DimPlot(data.ob,cells.highlight=(data.ob.umap %>% filter(cluster_EML=="HLC") %>% pull(cell)) )+theme_void()+NoLegend()+ggtitle("HLC")+theme(plot.title = element_text(hjust=0.5,face="bold"))
temp.plot$TLC <- DimPlot(data.ob,cells.highlight=(data.ob.umap %>% filter(cluster_EML=="TLC") %>% pull(cell)) )+theme_void()+NoLegend()+ggtitle("TLC")+theme(plot.title = element_text(hjust=0.5,face="bold"))
temp.plot$AMLC <- DimPlot(data.ob,cells.highlight=(data.ob.umap %>% filter(cluster_EML=="AMLC") %>% pull(cell)) )+theme_void()+NoLegend()+ggtitle("AMLC")+theme(plot.title = element_text(hjust=0.5,face="bold"))
temp.plot$MeLC <- DimPlot(data.ob,cells.highlight=(data.ob.umap %>% filter(cluster_EML=="MeLC") %>% pull(cell)) )+theme_void()+NoLegend()+ggtitle("MeLC")+theme(plot.title = element_text(hjust=0.5,face="bold"))
temp.plot$Und <- DimPlot(data.ob,cells.highlight=(data.ob.umap %>% filter(cluster_EML=="Undef") %>% pull(cell)) )+theme_void()+NoLegend()+ggtitle("Undef")+theme(plot.title = element_text(hjust=0.5,face="bold"))
print(cowplot::plot_grid(plotlist=temp.plot))

