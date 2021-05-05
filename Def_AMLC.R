#' ---
#' title: "Restore the cell type of AMLC dataset (Zheng et al., 2019)"
#' output: 
#'  html_document:
#'    code_folding: hide
#'    
#' ---


#  R3.6

#' reference:[Controlled modelling of human epiblast and amnion development using stem cells](https://www.nature.com/articles/s41586-019-1535-2)

# loading R library
#R3.6 
rm(list=ls())
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(pheatmap))
suppressMessages(library(Seurat))



# working directory
DIR <- "~/My_project/JP_project"
setwd(DIR)

#' Loading R functions
source("~/PC/R_code/functions.R")
source("~/PC/SnkM/SgCell.R")
feature.plot.cor <- c("yellow","red","black")

suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)

options(digits = 4)
options(future.globals.maxSize= 3001289600)
TD="D2_pub"



#' loading gene metadata
load("tmp_data/gene.meta.Rdata",verbose=T)

savefile=paste0("tmp_data/",TD,"/EPSC_Amnion.Rdata")

if (file.exists(savefile)) {
  load(savefile,verbose = T)
}else{
  load(paste0("tmp_data/",TD,"/all.counts.meta.Rdata"),verbose=T)
  
 
  #' loading published seurat object
  AMLC.ob <- readRDS("data/GSE134571/GSE134571_Posterior48h_H9_Amnion_Merged.rds")
  AMLC.ob <- AMLC.ob %>% RenameIdents('0'="Tsw-AMLC",'1'='MeLC2','2'='hES','3'='MeLC1','4'='hPGCLC','5'='AMLC')
  
  AMLC.DM <- FunRF_FindAllMarkers_para(AMLC.ob)
  temp.plot1 <- list()
  temp.plot1$id <- DimPlot(AMLC.ob,label=T) +NoAxes()+NoLegend()
  temp.plot1$oid <- DimPlot(AMLC.ob,label=T,group.by="orig.ident") +NoAxes()+NoLegend()
  temp.mk.gene <- (AMLC.DM$sig  %>% filter(NP=="pos") %>% filter(gene %in% rownames(counts.all)) %>% group_by(set) %>% top_n(2,power) %>% pull(gene))
  for (g in temp.mk.gene) {
    temp.plot1[[g]] <- FeaturePlot(AMLC.ob,g)+NoAxes()+NoLegend()+ggtitle(g)+FunTitle()
  }
  
  #same QC from the paper
  temp.raw <-  meta.all %>% filter(pj=="JPF2019")
  temp.M <- meta.all %>% filter(pj=="JPF2019" & EML=="Pos48ELS") %>% filter(nGene >=3200 & nGene <= 6200 & mt.perc < 0.06) %>% bind_rows(meta.all %>% filter(pj=="JPF2019" & EML=="H9AMN") %>% filter(nGene >=3600 & nGene <= 6400 & mt.perc < 0.06))
 
  counts.filter <- counts.all[setdiff(rownames(counts.all), mt.gene),temp.M$cell]
  DR.TPM.filter <- apply(counts.filter,2,function(x){x/sum(x)*1000000})
  
  # expressed genes
  temp.sel.expG <- rownames(counts.filter)[rowSums(counts.filter[,c(temp.M$cell)] >=1) >5]
  
  data.temp <-  CreateSeuratObject(counts.filter[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) %>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000, verbose = FALSE)  %>% ScaleData(verbose=F,vars.to.regress=c("mt.perc"))%>% RunPCA(verbose=F,npcs=20) %>% RunUMAP(dims=1:20,verbose=F,min.dist=0.4) %>% FindNeighbors( dims = 1:20,verbose = FALSE) %>%  FindClusters(resolution = 0.45,verbose = FALSE)  #"mt.perc",
  data.mod <- data.temp
  id <- Idents(data.mod) %>% as.vector()
  names(id) <- colnames(data.mod)
  id[colnames(data.mod)[data.mod@meta.data$EML=="H9AMN" & (id %in% c(0,5))]] <- "Tsw-AMLC"
  id[colnames(data.mod)[data.mod@meta.data$EML=="H9AMN" & (id %in% c(1))]] <- "hES"
  id[colnames(data.mod)[data.mod@meta.data$EML=="Pos48ELS" & (id %in% c(3))]] <- "hPGCLC"
  id[colnames(data.mod)[data.mod@meta.data$EML=="Pos48ELS" & (id %in% c(2))]] <- "MeLC1"
  id[colnames(data.mod)[data.mod@meta.data$EML=="Pos48ELS" & (id %in% c(4,7))]] <- "MeLC2"
  id[colnames(data.mod)[data.mod@meta.data$EML=="Pos48ELS" & (id %in% c(6))]] <- "AMLC"
  id[! id %in% c("Tsw-AMLC","hES","hPGCLC","MeLC1","MeLC2","AMLC")] <- "mix"
  Idents(data.mod) <- as.factor(id)
  data.mod <- subset(data.mod,cells=names(id)[id!="mix"])
  
  table(data.mod %>% Idents())

  save(data.mod,temp.sel.expG,temp.M,temp.raw,temp.plot1,AMLC.DM,temp.mk.gene,file=savefile)
}


#' ### repeat paper's figure
#+ fig.width=12,fig.height=12
cowplot::plot_grid(plotlist=temp.plot1)

#' ### my new annotation
#+ fig.width=12,fig.height=12
temp.plot2 <- list()
temp.plot2$id <- DimPlot(data.mod,label=T) +NoAxes()+NoLegend()
temp.plot2$oid <- DimPlot(data.mod,label=T,group.by="EML") +NoAxes()+NoLegend()
for (g in temp.mk.gene) {
  temp.plot2[[g]] <- FeaturePlot(data.mod,g)+NoAxes()+NoLegend()+ggtitle(g)+FunTitle()
}
cowplot::plot_grid(plotlist=temp.plot2)
#VlnPlot(data.mod,temp.mk.gene)

#' ### General distribution of nGene and mt.perc
#+ fig.width=9,fig.height=9
temp.plot <- list()
temp.plot$AMN1 <- temp.raw %>% filter(EML=="Pos48ELS")%>% ggplot +geom_point(mapping=aes(x=mt.perc,y=nGene,col=EML))+geom_vline(xintercept = 0.06,linetype = "dashed" )+geom_hline(yintercept = 6200,linetype = "dashed")+geom_hline(yintercept = 3200,linetype = "dashed")+theme_classic()+ggtitle("QC for 10X(Pos48ELS)")+FunTitle()
temp.plot$AMN2 <- temp.raw %>% filter(EML=="H9AMN")%>% ggplot +geom_point(mapping=aes(x=mt.perc,y=nGene,col=EML))+geom_vline(xintercept = 0.06,linetype = "dashed" )+geom_hline(yintercept = 6400,linetype = "dashed")+geom_hline(yintercept = 3600,linetype = "dashed")+theme_classic()+ggtitle("QC for 10X(H9AMN)")+FunTitle()
cowplot::plot_grid(plotlist=temp.plot,nrow=2,ncol=2)


#' #### saving the new meta data
AMN.meta <- temp.M %>% filter(cell %in% colnames(data.mod)) 
AMN.meta$EML <- Idents(data.mod)[AMN.meta$cell] %>% as.vector() 
saveRDS(AMN.meta,paste0("tmp_data/",TD,"/meta.AMN.further.rds"))

#' #### downsample dataset
FunMaSF <- function(temp,n) {
  set.seed(123)
  return(head(temp[sample(nrow(temp)),],n))
}
temp.out <- NULL

meta.AMN.further.ds <- split(AMN.meta,AMN.meta$EML) %>% lapply(function(x){return(FunMaSF(x,100))}) %>% do.call("bind_rows",.)
saveRDS(meta.AMN.further.ds,file=paste0("tmp_data/",TD,"/meta.AMN.further.ds.rds"))

sessionInfo()                    
                        


