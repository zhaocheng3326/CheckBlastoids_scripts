#' ---
#' title: "check the CS7 annotation of Epi,PriS and Amnion"
#' output: 
#'  html_document:
#'    code_folding: hide
#'    
#' ---
#' 
#' 
#' check the subType of CS7 dataset
# R4.0
rm(list=ls())

rewrite=FALSE
condaENV <- "/home/chenzh/miniconda3/envs/R4.0"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Seurat))
suppressMessages(library(SeuratWrappers))
suppressMessages(library(patchwork))

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
#lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/",MM,".lognormExp.mBN.rds"))
meta.filter<- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% filter(pj %in% c("EBD2","IBD2","JPF2019","D3post","CS7","SPH2016","Blakeley","nBGuo") ) %>% mutate(EML=ifelse(EML %in% c("Ectoderm"),"Amnion",EML))

#' loading the multi-data integration results
data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/Pub.D2.umap.cord.rds"))
data.ob <- readRDS(paste0("tmp_data/",TD,"/Pub.D2.mnn.rds"))
#data.ob@meta.data$seurat_clusters=(data.ob.umap %>% tibble::column_to_rownames("cell"))[colnames(data.ob),"seurat_clusters"]
#DimPlot(data.ob,group.by="seurat_clusters",label=T)+NoAxes()+NoLegend()


#' loading CS7+NHP cross-species integration dataset
data.CS7.NHP.ob <- readRDS(paste0("tmp_data/",TD,"/NHP.","CS7",".mnn.rds"))

sus.pris.amnion.cells <- c("sc7788294","sc7786321","sc7785994","sc7786393","sc7788314","sc7786320","sc7785681","sc7788101","sc7788122","sc7786111") 

PGC.marker <- c("NANOS3", "SOX17", "DND1","LAMA4", "DPPA5")
#' check the location of sus cells in cross-speciese integration
temp.plot <- list()
temp.plot$EML <- DimPlot(data.CS7.NHP.ob,group.by="EML",label=T)+NoAxes()+NoLegend()
temp.plot$Amnion1 <- DimPlot(data.CS7.NHP.ob,cells.highlight=colnames(data.CS7.NHP.ob)[data.CS7.NHP.ob@meta.data$EML=="Amnion1"])+NoAxes()+NoLegend()+ggtitle("Amnion1(Early,NHP)")
temp.plot$PriS <- DimPlot(data.CS7.NHP.ob,cells.highlight=colnames(data.CS7.NHP.ob)[data.CS7.NHP.ob@meta.data$EML=="PriS" & data.CS7.NHP.ob@meta.data$pj=="CS7"])+NoAxes()+NoLegend()+ggtitle("PriS(human)")
temp.plot$sus.pris.amnion.cells  <- DimPlot(data.CS7.NHP.ob,cells.highlight=sus.pris.amnion.cells )+NoAxes()+NoLegend()+ggtitle("sus.pris.amnion.cells")
cowplot::plot_grid(plotlist=temp.plot,ncol=2)



#' check the whole CS7 data integration
temp.M <- meta.filter %>% filter(pj %in% c("CS7"))

temp.cell <- temp.M %>% pull(cell)
temp.sel.expG <- rownames(counts.filter )[rowSums(counts.filter[,temp.cell] >=1) >=5]

verbose <- F
data.temp <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)  %>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>% ScaleData(verbose=F)%>% RunPCA(verbose=F,npcs=20) %>% RunUMAP(dims=1:20,verbose=F) %>% FindNeighbors( dims = 1:20,verbose = FALSE) %>% FindClusters(reso=0.3,verbose=F)
data.temp@meta.data$big_seurat_clusters <- (data.ob.umap %>% tibble::column_to_rownames("cell"))[rownames(data.temp@meta.data),"seurat_clusters"]


temp.plot <- list()
temp.plot$raw <- DimPlot(data.temp,label=T)+theme_void()+NoLegend()
temp.plot$EML <- DimPlot(data.temp,group.by="EML",label=T)+theme_void()+NoLegend()
#temp.plot$new_EML <- DimPlot(data.temp,group.by="subCT",label=T)+theme_void()+NoLegend()
#temp.plot$st <- DimPlot(data.temp,group.by="seqType",label=T)+theme_void()
temp.plot$bs <- DimPlot(data.temp,group.by="big_seurat_clusters",label=T)+theme_void()
temp.plot$PGC <- DimPlot(data.temp,cells.highlight=(data.ob.umap %>% filter(subCT=="PGC") %>% pull(cell)))+NoAxes()+NoLegend()+ggtitle("PGC")+FunTitle()
temp.plot$CS7_C11 <- DimPlot(data.temp,cells.highlight=(data.ob.umap %>% filter(seurat_clusters=="C11" & pj=="CS7") %>% pull(cell)))+NoAxes()+NoLegend()+ggtitle("CS7-C11")+FunTitle()
temp.plot$sus.pris.amnion <- DimPlot(data.temp,cells.highlight=sus.pris.amnion.cells)+NoAxes()+NoLegend()+ggtitle("sus.pris.amnion")+FunTitle()
cowplot::plot_grid(plotlist=temp.plot,ncol=2) # not used
#rm(data.temp)

#' check the Amnion , PriS, And Epi cells only 
temp.M <-  meta.filter %>% filter(pj %in% c("CS7")) %>% filter(EML %in% c("Epiblast","PriS","Amnion")) %>% filter(subCT !="PGC") %>% inner_join(data.ob.umap %>% select(cell,seurat_clusters) %>% rename(big_seurat_clusters=seurat_clusters),by="cell") 
temp.cell <-temp.M$cell
temp.sel.expG <- rownames(counts.filter )[rowSums(counts.filter[,temp.cell] >=1) >=3]
verbose <- F
data.AEP <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE)  %>% FindVariableFeatures( selection.method = "vst", nfeatures = 500, verbose = FALSE) %>% ScaleData(verbose=F)%>% RunPCA(verbose=F,npcs=30) %>% RunUMAP(dims=1:15,verbose=F) %>% FindNeighbors( dims = 1:15,verbose = FALSE)  %>% FindClusters(reso=0.3,verbose=F)

temp.M$subCluster <- paste0("SC",as.vector(Idents(data.AEP)[temp.M$cell]))
temp.M %>% select(subCluster,EML) %>% table()
temp.M <- temp.M %>% mutate(new_EML=EML) %>% mutate(new_EML=ifelse(EML %in% c("PriS") & subCluster=="SC0", "EPI_PriS",new_EML)) %>% mutate(new_EML=ifelse(EML %in% c("PriS") & subCluster=="SC3", "PriS(Amnion-like)",new_EML))%>% mutate(new_EML=ifelse(EML %in% c("Epiblast") & subCluster=="SC3", "Epiblast(Amnion-like)",new_EML))%>% mutate(new_EML=ifelse(EML %in% c("Epiblast","PriS") & subCluster=="SC2", "EPI_PriS",new_EML))

data.AEP@meta.data$new_EML <- (temp.M %>% tibble::column_to_rownames("cell"))[rownames(data.AEP@meta.data),"new_EML"]
#' check the distribution
temp.plot <- list()
temp.plot$raw <- DimPlot(data.AEP,label=T)+theme_void()+NoLegend()
temp.plot$EML <- DimPlot(data.AEP,group.by="EML",label=T)+theme_void()+NoLegend()+FunTitle()+ggtitle("old annotation")
temp.plot$new_EML <- DimPlot(data.AEP,group.by="new_EML",label=T)+theme_void()+NoLegend()+FunTitle()+ggtitle("new annotation")
temp.plot$bs <- DimPlot(data.AEP,group.by="big_seurat_clusters",label=T)+theme_void()+FunTitle()
cowplot::plot_grid(plotlist=temp.plot,nrow=2,ncol=2) 


#' get the pure lineage cell
data.deg <- subset(data.AEP,cells=(temp.M %>% filter(EML==new_EML) %>% pull(cell)))

#' get top DE among pure Epi,Amnion and PriS 
Idents(data.deg) <- factor(data.deg@meta.data$EML)
temp.DM <-  FunRF_FindAllMarkers_para(data.deg)

#' create module score based on top mk genes
temp.mk.list <-  temp.DM$sig %>% group_by(set) %>% filter(power >0.7) %>% top_n(25,power) 
temp.mk.list$set %>% table()
temp.mk.list
temp.mk.list <- temp.mk.list$gene %>% split(temp.mk.list$set)
#' add gene module score
names(temp.mk.list) <- paste(names(temp.mk.list),"Score",sep="_")
#' check gene module score for original Amnion, Epi and PriS
data.AEP <- data.AEP %>% AddModuleScore(temp.mk.list,name=names(temp.mk.list))
colnames(data.AEP@meta.data)[grepl("Amnion_Score",colnames(data.AEP@meta.data))] <- "Amnion_Score"
colnames(data.AEP@meta.data)[grepl("Epiblast_Score",colnames(data.AEP@meta.data))] <- "Epiblast_Score"
colnames(data.AEP@meta.data)[grepl("PriS_Score",colnames(data.AEP@meta.data))] <- "PriS_Score"
#' check the module score
VlnPlot(data.AEP,c("Amnion_Score","Epiblast_Score","PriS_Score"),group.by="new_EML")
#data.AEP@meta.data[data.AEP@meta.data$new_EML=="Epi(PriS-like)",c("Amnion_Score","Epiblast_Score","PriS_Score")]

# detail check and modify the identities
temp.meta <- data.AEP@meta.data %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:big_seurat_clusters,Amnion_Score:PriS_Score)%>% rows_update(temp.M %>% mutate(EML=new_EML)%>%  select(cell,EML),by="cell")

#+ fig.width=12,fig.height=6
print(
  (temp.meta %>% ggplot+geom_point(mapping=aes(x=Amnion_Score,y=PriS_Score,col=EML))) + (temp.meta %>% ggplot+geom_point(mapping=aes(x=Amnion_Score,y=Epiblast_Score,col=EML)))+(temp.meta %>% ggplot+geom_point(mapping=aes(x=PriS_Score,y=Epiblast_Score,col=EML)))+ plot_layout(guides = 'collect')
)

#' update meta data information
temp.M <- temp.M  %>% rows_update(temp.meta %>% filter(Amnion_Score >1 & PriS_Score <0) %>% mutate(new_EML="Amnion")%>% select(cell,new_EML) ,by="cell")
temp.M <- temp.M  %>% rows_update(temp.meta %>% filter(EML %in% c("EPI_PriS"))%>% filter(Epiblast_Score >0 & PriS_Score <0) %>% mutate(new_EML="Epiblast")%>% select(cell,new_EML) ,by="cell")
temp.M <- temp.M  %>% rows_update(temp.meta %>% filter(EML %in% c("EPI_PriS"))%>% filter(Epiblast_Score <0 & PriS_Score >0) %>% mutate(new_EML="PriS")%>% select(cell,new_EML) ,by="cell")

#' check the heatmap of top marker genes
zs.limit=2.5
temp.exp <- data.AEP@assays$RNA@counts[unlist(temp.mk.list),]
temp.sel.exp <- t(apply(temp.exp,1,scale))
colnames(temp.sel.exp) <- colnames(temp.exp)
rownames(temp.sel.exp) <- rownames(temp.exp)
temp.sel.exp[temp.sel.exp>zs.limit] <- zs.limit
temp.sel.exp[temp.sel.exp<  (-1*zs.limit)] <- -1*zs.limit

temp.anno <- (temp.M %>% select(cell,EML,new_EML)%>% mutate(od=factor(new_EML,c("Amnion","Epiblast(Amnion-like)","PriS(Amnion-like)","Epiblast","EPI_PriS","PriS"),ordered=T))%>% arrange(od,EML) %>% mutate(old.EML=EML,EML=new_EML)%>% mutate(EML=ifelse(EML=="Epiblast(Amnion-like)","EPI_Amnion",EML))%>% mutate(EML=ifelse(EML=="PriS(Amnion-like)","PriS_Amnion",EML))  %>% select(cell,old.EML,EML) %>% tibble::column_to_rownames("cell")) # %>% mutate(EML=ifelse(EML=="Epi(PriS-like)","PriS",EML))

#+ fig.width=12,fig.height=9
pheatmap::pheatmap( temp.sel.exp[,rownames(temp.anno )],scale="row",show_rownames=T,show_colnames=F,cluster_cols=F,cluster_rows=F,main="",annotation_col = temp.anno )

#' update CS7 annotation
CS7.meta <- meta.filter %>% filter(pj=="CS7")
APE.meta <- temp.M %>% mutate(EML=new_EML) %>% select(cell,EML) %>% mutate(EML=ifelse(EML=="Epiblast(Amnion-like)","EPI_Amnion",EML))%>% mutate(EML=ifelse(EML=="PriS(Amnion-like)","PriS_Amnion",EML)) 
CS7.meta <- CS7.meta %>% rows_update(APE.meta %>% select(cell,EML) ,by="cell")

data.AEP.UMAP <- data.AEP@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% inner_join( meta.filter,by="cell") %>% inner_join(CS7.meta %>% mutate(new_EML=EML) %>% select(cell,new_EML),by="cell" )
data.AEP.pure.DE <- temp.DM

if (!file.exists(paste0("tmp_data/",TD,"/CS7.updata.anno.rds")) | rewrite) {
  print("save output")
  saveRDS(CS7.meta, file=paste0("tmp_data/",TD,"/CS7.updata.anno.rds"))
  save(data.AEP.UMAP,data.AEP.pure.DE,file=paste0("tmp_data/",TD,"/check.CS7.sub.Rdata"))
}
