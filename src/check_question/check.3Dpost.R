#' ---
#' title: "check the Xiang et al dataset"
#' output: 
#'  html_document:
#'    code_folding: hide
#'    
#' ---


#R3.6
rm(list=ls())

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratWrappers))


# working directory
DIR <- "~/My_project/CheckBlastoids"
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
TD="Nov14_2021"
MM="Pub"

nGene <- 2000;pc=30

rename <- dplyr::rename
select<- dplyr::select
filter <- dplyr::filter

#' loading data
load("tmp_data/gene.meta.Rdata",verbose=T)
counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))

meta.filter<- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% filter(pj %in% c("EBD2","IBD2","JPF2019","D3post","CS7","SPH2016","Blakeley","nBGuo") ) %>% mutate(EML=ifelse(EML %in% c("Ectoderm"),"Amnion",EML))


temp.M <- meta.filter %>% filter(pj %in% c("D3post"))
temp.cell <- temp.M$cell
temp.sel.expG <- rownames(counts.filter )[rowSums(counts.filter[,temp.cell] >=1) >=5]

data.temp <-  CreateSeuratObject(counts.filter[temp.sel.expG,temp.M$cell], meta.data = (temp.M%>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) %>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000, verbose = FALSE)  %>% ScaleData(verbose=F)%>% RunPCA(verbose=F,npcs=20) %>% RunUMAP(dims=1:20,verbose=F) 
Idents(data.temp) <- data.temp@meta.data$EML

DimPlot(data.temp,group.by="EML",label=T)+theme_void()+NoLegend()
DimPlot(data.temp,group.by="devTime")+theme_void()
FeaturePlot(data.temp,c("GATA2","GATA3","SOX2","NANOG"))
FeaturePlot(data.temp,c("HES1","SOX11"))

VlnPlot(data.temp,c("HES1","SOX11"))


temp.DM <- FunRF_FindAllMarkers_para(data.temp)

temp.DM.sel <- temp.DM$sig %>% filter(NP=="pos") %>% group_by(set)  %>% top_n(15,power)%>% ungroup() 
Lineage.markers <- split(temp.DM.sel,temp.DM.sel$set)
Lineage.markers <- Lineage.markers[c("EPI","PrE","PSA-EPI","ICM","EVT","CTB","STB")]
temp.exp <- data.temp@assays$RNA@data[unlist(lapply(Lineage.markers,function(x){return(x$gene)})),] %>% as.data.frame()
temp.sel.exp <- t(apply(temp.exp,1,scale))
colnames(temp.sel.exp) <- colnames(temp.exp)
rownames(temp.sel.exp) <- rownames(temp.exp)
temp.sel.exp[temp.sel.exp>3.5] <- 3.5
temp.sel.exp[temp.sel.exp<  -3.5] <- -3.5
temp.anno <- temp.M %>%select(cell,pj,devTime,EML) %>% mutate(EML=factor(EML,levels=c("EPI","PrE","PSA-EPI","ICM","EVT","CTB","STB"),ordered = T)) %>% arrange(EML,devTime) %>% tibble::column_to_rownames("cell")

#+ fig.width=18,fig.height=18
pheatmap::pheatmap(temp.sel.exp[,rownames(temp.anno)],cluster_rows=F,cluster_cols=F,scale="none",annotation_col=temp.anno,show_colnames=F,color=(colorRampPalette(c("royalblue3","white","firebrick4"))(50)), fontsize_row=7,width=7,height=3,border_color="NA",main="Top lineage markers",gaps_row = unlist(lapply(Lineage.markers,function(x){return(nrow(x))})) %>% cumsum())
