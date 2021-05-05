#' ---
#' title: "Restore the cell type of Iblastoids dataset (Liu et al., 2021) "
#' output: 
#'  html_document:
#'    code_folding: hide
#'    
#' ---


#  R3.6

#' reference:[Modelling human blastocysts by reprogramming fibroblasts into iBlastoids](https://www.nature.com/articles/s41586-021-03372-y) 
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

suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)

options(digits = 4)
options(future.globals.maxSize= 3001289600)

TD="D2_pub"

savefile1=paste0("tmp_data/",TD,"/IBD2.Rdata")
savefile2=paste0("tmp_data/",TD,"/IBD2.counts.meta.Rdata")
if (file.exists(savefile1)) {
  load(savefile1,verbose = T)
  load(savefile2,verbose=T)
}else{
  #' loading raw data
  load(paste0("tmp_data/",TD,"/all.counts.meta.Rdata"),verbose=T)
  load("tmp_data/gene.meta.Rdata",verbose=T)
  
  
  temp.mk.gene <- list()
  temp.mk.gene$IB_Epi <- c("DPPA5","ALPPL2","POU5F1","MT1G","UTF1")
  temp.mk.gene$IB_TE <- c("WFDC2","FABP7","KRT19","EZR","TMEM54")
  temp.mk.gene$IB_PE <- c("APOA1","APOA2","CKB","APOC1","LCN15")
  temp.mk.gene$IB_IIM1 <- c("VIM","CD70","POU5F1B","TUBB2A","GSN")
  temp.mk.gene$IB_IM2 <- c("PITX1","RGS5","CTC-276P9.1","LIX1","HGF")
  temp.mk.gene$IB_IM3 <- c("IFIT1","IFIT3","CCL5","IFIT2","OASL")
  temp.mk.gene$IB_NR <- c("IFI6","HTRA1","DCN","SOD3","NBL1")
  # above marker genes are from 
  #https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-03372-y/MediaObjects/41586_2021_3372_MOESM5_ESM.xlsx
  
  #QC
  temp.raw <-  meta.all %>% filter(pj=="IBD2")
  temp.M <- meta.all %>% filter(pj=="IBD2") %>% filter(nGene >2000 & nGene <5000 & mt.perc < 0.125) 
  counts.filter <- counts.all[setdiff(rownames(counts.all), mt.gene),temp.M$cell]
  IBD2.counts <-  counts.filter 
  IBD2.meta <- temp.M
  
  
  # expressed genes
  temp.sel.expG <- rownames(counts.filter)[rowSums(counts.filter[,c(temp.M$cell)] >=1) >5]
  
  data.temp <-  CreateSeuratObject(counts.filter[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) %>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000, verbose = FALSE)  %>% ScaleData(verbose=F,vars.to.regress=c("mt.perc"))%>% RunPCA(verbose=F,npcs=20) %>% RunUMAP(dims=1:20,verbose=F,min.dist=0.4) %>% FindNeighbors( dims = 1:20,verbose = FALSE) %>%  FindClusters(resolution = 0.2,verbose = FALSE)  #"mt.perc",
  
  data.mod <- data.temp %>% RenameIdents(c('0'='IB_IM1','1'='IB_PE','2'='IB_TE','3'='IB_Epi','4'='IB_IM2','5'='IB_PE','6'='IB_IM3','7'='IB_NR'))
  #assign the top right corner cells to uncertain cells
  uc.cell <- colnames(data.mod)[Idents(data.mod)=="IB_IM3" & data.mod@reductions$umap@cell.embeddings[,1] > 5]
  data.mod <- SetIdent(data.mod,cells=uc.cell,value="IB_uc")
  
  table(data.mod %>% Idents())
  save(IBD2.counts,IBD2.meta,file=savefile2)
  save(data.mod,temp.sel.expG,temp.M,temp.raw,temp.mk.gene,uc.cell,file=savefile1)
}

#DR.TPM.filter <- apply(counts.filter,2,function(x){x/sum(x)*1000000})
#DimPlot(data.mod,cells.highlight=uc.cell)


#' repeat paper's figure
#+ fig.width=12,fig.height=12
print("https://www.nature.com/articles/s41586-021-03372-y.pdf fig2C")
print("marker genes are extrated from https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-021-03372-y/MediaObjects/41586_2021_3372_MOESM5_ESM.xlsx")
#' my new annotation
#+ fig.width=12,fig.height=12
temp.plot2 <- list()
temp.plot2$id <- DimPlot(data.mod,label=T) +NoAxes()+NoLegend()
for (g in names(temp.mk.gene)) {
  temp.plot2[[g]] <- FeaturePlot(data.mod,temp.mk.gene[[g]][1])+NoAxes()+NoLegend()+ggtitle(paste(g,temp.mk.gene[[g]][1],sep="\n"))+FunTitle()
}
cowplot::plot_grid(plotlist=temp.plot2)
#VlnPlot(data.mod,temp.mk.gene)

#' ### General distribution of nGene and mt.perc
#+ fig.width=6,fig.height=6
print(
  temp.raw %>% ggplot +geom_point(mapping=aes(x=mt.perc,y=nGene,col=EML))+geom_vline(xintercept = 0.125,linetype = "dashed" )+geom_hline(yintercept = 2000,linetype = "dashed")+geom_hline(yintercept = 5000,linetype = "dashed")+theme_classic()+ggtitle("QC for 10X")+FunTitle()
)

#temp.M <- meta.all %>% filter(pj=="JPF2019" & EML=="Pos48ELS") %>% filter(nGene >=3200 & nGene <= 6200 & mt.perc < 0.06) %>% bind_rows(meta.all %>% filter(pj=="JPF2019" & EML=="H9AMN") %>% filter(nGene >=3600 & nGene <= 6400 & mt.perc < 0.06))

#' saving the new meta data
IBD2.meta <- IBD2.meta %>% filter(cell %in% colnames(data.mod)) 
IBD2.meta$EML <- Idents(data.mod)[IBD2.meta$cell] %>% as.vector() 
saveRDS(IBD2.meta,paste0("tmp_data/",TD,"/meta.IBD2.further.rds"))

#' downsample dataset

temp.out <-split(IBD2.meta,IBD2.meta$EML)

meta.IBD2.further.ds <- temp.out[c("IB_Epi","IB_PE","IB_TE")] %>% lapply(function(x){return(FunMaSF(x,150))}) %>% do.call("bind_rows",.) %>% bind_rows( temp.out[c("IB_IM1","IB_IM2","IB_IM3")] %>% lapply(function(x){return(FunMaSF(x,100))}) %>% do.call("bind_rows",.))
saveRDS(meta.IBD2.further.ds ,file=paste0("tmp_data/",TD,"/meta.IBD2.further.ds.rds"))

                        
                        
#' ### check the expression of marker gene for Ectoderm(Amnion)
#+ fig.width=9,fig.height=9
temp.plot <- list()
for (g in c("TRIML2","ISL1","STOM","TMEM54","CTSV","GABRP","MSX2","EPAS1","ANXA3","DLX5","HEY1","KRT7","SERPINB9","ATP2B1","TFAP2A"  )) {
  temp.plot[[g]] <- FeaturePlot(data.mod,g)+NoAxes()+NoLegend()+ggtitle(paste(g))+FunTitle()
}
cowplot::plot_grid(plotlist=temp.plot[1:15])

#' ### check marker gene expression for TE lineages
#+ fig.width=9,fig.height=9
temp.plot <- list()
for (g in c("GATA3","GATA2","CDX2","KRT18","TEAD3","CGA","CCR7","PDF","GCM1"  )) {
  temp.plot[[g]] <- FeaturePlot(data.mod,g)+NoAxes()+NoLegend()+ggtitle(paste(g))+FunTitle()
}
cowplot::plot_grid(plotlist=temp.plot)


#' ### check other marker genes
#+ fig.width=9,fig.height=9
temp.plot <- list()
for (g in c("CCR7","CGA","CSF3R","CYP19A1","DLX5","GCM1","GREM2","KRT23","MLLT11","MRGPRX1","MUC15","MYO10","NRP1","OVOL1","PGF","PTPRE","S1PR2","SDC1","SERPINB9","SLCO4A1")) {
  if (g %in% rownames(data.mod@assays$RNA@counts)) {
    temp.plot[[g]] <- FeaturePlot(data.mod,g)+NoAxes()+NoLegend()+ggtitle(paste(g))+FunTitle()
  }else{
    temp.plot[[g]] <-ggplot() +  annotate("text", x = 4, y = 25, size=2, label = c("not expressed(<6)")) + theme_void()+ggtitle(paste(g))+FunTitle()
  }
 
}
cowplot::plot_grid(plotlist=temp.plot)




