#' ---
#' title: "define the cell type of Stem-Blastoids(Yu et al., 2021) "
#' output: 
#'  html_document:
#'    code_folding: hide
#'    
#' ---


#  R3.6

#' reference:[Blastocyst-like structures generated from human pluripotent stem cells](https://www.nature.com/articles/s41586-021-03356-y%C2%A0)


#R3.6 
rm(list=ls())
rewrite=FALSE

condaENV <- "/home/chenzh/miniconda3/envs/R3.6"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

#' loading R library
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(edgeR))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(pheatmap))
suppressMessages(library(Seurat))


# set working directory
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


savefile=paste0("tmp_data/",TD,"/EBD2.Rdata")

if (file.exists(savefile)) {
  load(savefile,verbose = T)
}else{
  
  # loading all counts
  load(paste0("tmp_data/",TD,"/all.counts.meta.Rdata"),verbose=T)
  load("tmp_data/gene.meta.Rdata",verbose=T)
  
  
  # QC 
  EBD2.meta <- meta.all %>% filter(pj =="EBD2") %>% filter(nGene >1000 & nGene <6500 & mt.perc < 0.125) 
  EBD2.counts <- counts.all[setdiff(rownames(counts.all), mt.gene),EBD2.meta$cell]
 
  # using previous Petropolous annotation (from https://pubmed.ncbi.nlm.nih.gov/29361568/)
  temp.anno <- readxl::read_xlsx(("data/stirparo2018_tableS4.xlsx"), sheet = 1) %>% tbl_df() %>% filter(Study == "Petropoulos et al., 2016 (ERP012552)")%>% select("Cell","Revised lineage (this study)")
  colnames(temp.anno) <- c("cell","EML")
  temp.anno <- temp.anno %>% mutate(cell=gsub("_", ".", cell)) %>% mutate(EML=gsub("undefined","Undef",EML))%>% mutate(EML=gsub("epiblast","Epi",EML))%>% mutate(EML=gsub("Inner cell mass","ICM",EML))%>% mutate(EML=gsub("intermediate","INT",EML))%>% mutate(EML=gsub("primitive_endoderm","PrE",EML))%>% mutate(EML=gsub("trophectoderm","TE",EML))
  SPH.meta <- meta.all %>% filter(pj=="SPH2016") %>% filter(nGene >2000 & mt.perc < 0.125) %>% select(-EML) %>% inner_join(temp.anno,by="cell")
  
  meta.filter <- bind_rows(EBD2.meta,SPH.meta )
  counts.filter <-   counts.all[setdiff(rownames(counts.all), mt.gene),meta.filter $cell]
  
  # release memory
  rm(list=c("counts.all","meta.all"))
  
  temp.meta <- list()
  temp.meta$LW36 <- meta.filter %>% filter(cellType=="LW36")
  temp.meta$LW60 <- meta.filter %>% filter(cellType=="LW60")
  temp.meta$LW61 <- meta.filter %>% filter(cellType=="LW61")
  temp.meta$sph <- meta.filter %>% filter(pj=="SPH2016")
  
  s.genes <- cc.genes$s.genes
  g2m.genes <- cc.genes$g2m.genes
  
  data.list <- list()
  
  n="LW60"
  temp.M <- temp.meta[[n]]
  temp.sel.expG <- rownames(counts.filter)[rowSums(counts.filter[,c(temp.M$cell)] >=1) >5]
  data.temp <- CreateSeuratObject(counts.filter[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) %>% FindVariableFeatures( selection.method = "vst", nfeatures = 2000, verbose = FALSE) %>% CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes)  %>% ScaleData(verbose=F,vars.to.regress=c("mt.perc"))%>% RunPCA(verbose=F,npcs=30) %>% RunUMAP(dims=1:20,verbose=F) %>% FindNeighbors( dims = 1:20,verbose = FALSE) %>%  FindClusters(resolution = 0.6,verbose = FALSE)
  data.list[[n]] <- RenameIdents(data.temp, '0' = 'EB_TLC','1'='EB_ELC','2'="EB_TLC",'3'='EB_U5','4'='EB_ELC','5'="EB_ELC",'6'="EB_ELC",'7'="EB_U6",'8'='EB_TLC','9'="EB_U10",'10'="EB_HLC")
  EBD2.LW60.meta <- data.frame(EML=Idents(data.list[[n]]),cell=colnames(data.list[[n]])) %>% tbl_df() %>% mutate_all(as.vector) %>% inner_join(EBD2.meta %>% select(-EML),by="cell") %>% mutate(subCT=EML ) %>% mutate(subCT=ifelse(cell %in% (colnames(data.temp)[Idents(data.temp)==0]),'EB_TLC0',subCT))%>% mutate(subCT=ifelse(cell %in% (colnames(data.temp)[Idents(data.temp)==2]),'EB_TLC49',subCT))%>% mutate(subCT=ifelse(cell %in% (colnames(data.temp)[Idents(data.temp)==8]),'EB_TLC8',subCT))
  
  data.list[[n]]@meta.data$subCT <- (EBD2.LW60.meta %>% tibble::column_to_rownames("cell"))[rownames( data.list[[n]]@meta.data),"subCT"]
                                                 
  #' do integration for LW60, LW61 and sph2016
  temp.M <- meta.filter %>% filter(cellType %in% c("LW60","LW61") | pj=="SPH2016") %>% mutate(pj=ifelse(cellType == "LW60","EBD2_LW60",pj)) %>%  mutate(pj=ifelse(cellType == "LW61","EBD2_LW61",pj))

  temp.sel.expG <- rownames(counts.filter)[rowSums(counts.filter[,c(temp.M$cell)] >=1) >5]
  data.spt <- CreateSeuratObject(counts.filter[temp.sel.expG,temp.M$cell], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) %>% CellCycleScoring(s.features = s.genes, g2m.features = g2m.genes) %>% SplitObject( split.by = "pj") %>% lapply(function(x){x=FindVariableFeatures(x,verbose=F)})
  
  nf=1000
  k.filter=200
  pc1=pc2=12
  
  data.it.features <- SelectIntegrationFeatures(object.list = data.spt, nfeatures = nf)
  
  print(paste(nf,k.filter,pc1))
  pc2=pc1
  data.it.anchors <- FindIntegrationAnchors(object.list = data.spt,  anchor.features = data.it.features, verbose = FALSE,dims=1:pc1, k.filter = k.filter)
  data.it <- IntegrateData(anchorset = data.it.anchors, verbose = FALSE)
  DefaultAssay(object = data.it) <- "integrated"
  data.it <- data.it %>% ScaleData(verbose = FALSE,vars.to.regress=c("mt.perc"))%>% RunPCA( verbose = FALSE,npcs =30) %>%  RunUMAP(dims = 1:pc2,verbose = FALSE)%>% FindNeighbors( dims = 1:pc2,verbose = FALSE,nn.method="annoy",annoy.metric="cosine") %>% FindClusters(resolution = 0.6, verbose = FALSE)
  data.ob <- RenameIdents(data.it, '0' = 'EB_TLC','1'='EB_ELC','2'="EB_ELC",'3'='EB_ELC','4'='EB_TLC','5'="EB_U5",'6'="EB_U6",'7'="EB_ELC",'8'='EB_TLC','9'="EB_TLC",'10'="EB_U10",'11'="EB_HLC",'12'='EB_UN34','13'='EB_UN34')
 
  DefaultAssay(data.ob) <- "RNA"
  data.ob <- NormalizeData(data.ob, verbose = FALSE)
  
  data.list$it <- data.ob
  data.list$it.LW60 <- subset(data.ob,cells=c(temp.meta$LW60$cell))
  
  #' assign the identity to the meta.data
  EBD2.meta <- data.frame(EML=Idents(data.list$it),cell=colnames(data.list$it)) %>% tbl_df() %>% mutate_all(as.vector) %>% inner_join(EBD2.meta %>% select(-EML),by="cell") %>% mutate(subCT=EML ) %>% mutate(subCT=ifelse(cell %in% (colnames(data.it)[Idents(data.it)==0]),'EB_TLC0',subCT))%>% mutate(subCT=ifelse(cell %in% (colnames(data.it)[Idents(data.it)==4]),'EB_TLC4',subCT)) %>% mutate(subCT=ifelse(cell %in% (colnames(data.it)[Idents(data.it)==8]),'EB_TLC8',subCT))%>% mutate(subCT=ifelse(cell %in% (colnames(data.it)[Idents(data.it)==9]),'EB_TLC9',subCT))
  
  save(data.list,meta.filter,counts.filter,EBD2.meta,EBD2.LW60.meta,file=savefile)
}

#' #### check cells' distribution
temp.plot <- list()
temp.plot$cl <- DimPlot(data.list$it,label=T)+NoAxes()+NoLegend()+ggtitle("LW61/LW60/sph")+FunTitle()
temp.plot$sph <- DimPlot(data.list$it,cells.highlight=colnames(data.list$it)[data.list$it@meta.data$pj=="SPH2016"])+NoAxes()+NoLegend()+ggtitle("SPH2016")+FunTitle()
temp.plot$stage <- DimPlot(data.list$it,cells=colnames(data.list$it)[data.list$it@meta.data$pj=="SPH2016"],group.by="devTime",label=T)+NoAxes()+NoLegend()+ggtitle("SPH2016")+FunTitle()
temp.plot$EML <- DimPlot(data.list$it,cells=colnames(data.list$it)[data.list$it@meta.data$pj=="SPH2016"],group.by="EML",label=T)+NoAxes()+NoLegend()+ggtitle("SPH2016")+FunTitle()
temp.plot$LW60 <- DimPlot(data.list$it.LW60,label=T)+NoAxes()+NoLegend()+ggtitle("it.LW60")+FunTitle()

#" #### check reported marker gene expression
for (g in c("SOX2","NODAL","NANOG","FOXA2","GATA4","COL4A1","GATA2","GATA3","CCR7","CYP19A1","MUC15")) {
  temp.plot[[g]] <- FeaturePlot(data.list$it.LW60,g)+ggtitle(paste(g))+FunTitle()+NoAxes()+NoLegend()
}
#+ fig.width=9,fig.height=9
print(cowplot::plot_grid(plotlist=temp.plot[1:4])) 
#+ fig.width=12,fig.height=12
print(cowplot::plot_grid(plotlist=temp.plot[5:16])) 


#' ### check cell proportation for each cluster
EBD2.meta %>% filter(cellType=="LW60") %>% pull(EML) %>% table()

#' check the marker gene expression from their paper
# above marker genes are from https://htmlpreview.github.io/?https://github.com/jlduan/Human_blastoid/blob/master/notebooks/analyze_blastoids.html (genes after GATA2/3 are removed)
temp.mk.gene <- list()
temp.mk.gene$EB_Epi <- c("POU5F1","GDF3","DPPA5","SOX2","NANOG","TDGF1","NODAL","AC064802.1","C19orf85","VENTX","VSNL1","CBR3","SFRP2","ETV1","PRDM14")%>% intersect(rownames(counts.filter))
temp.mk.gene$EB_PE <- c("GATA4","SPINK1","CACHD1","PITX2","NID2","SLC6A19","KLB","SOAT2","APOC3","CDH2","FLRT3","RNASE1","LGALS2","PDGFRA","SOX17","IGF1","HNF1B","APOA1","MYL4","FGG","FOXA2","RSPO3","FN1","RHOBTB3","GATA6","COL4A1")%>% intersect(rownames(counts.filter))
temp.mk.gene$EB_TE <- c("MME","PPP3CA","FYB1","ATP7B","GRAMD2A","ISM2","GJA5","ERVMER34-1","ASCL2","PLEKHA6","TMC5","ITGB4","C1orf105","SOD3","ECM1","HLA-G","LIFR","TAP1","CCBE1","DUSP13","FABP3","CADM1","COL1A2","PTK7","CLDN10","RNF43","MEX3A","HPGD","MAGI3","EGR1","MMD","PHACTR1","PRTG","PRSS12","CRABP2","HAPLN1","HAND1","PEG3","BMP7","IGDCC3","RFX1","CEP170B","ABCD1","ACKR2","CA12","FOSB","NPB","UCA1","DAB2","CCR7","DLX5","COL21A1","TULP1","INSL4","MYCNUT","SLC13A4","ADAMTS20","ESRRG","GREM2","MRGPRX1","AC010328.1","OVOL1","GCM1","CYP19A1","CPM","HOPX","LGALS16","KRT23","MUC15","LCMT1-AS2","ERVV-1","TNFRSF1B","TREM1","RHOBTB1","TNFAIP2","MICALL2","ERVFRD-1","PLAC4","KRT18","KRT8","CLIC3","GATA3","GATA2")%>% intersect(rownames(counts.filter))
 #"AL162231.1","CCDC57","BTBD2","PACS2","C2orf74","CD70","UBE2L6","TP53TG1","ACTG2","ACTA1","IFI27L2","LINC01405","FBXO2","LAMB3","PLAU","EMP3","PLAUR","MMP9","CD44","PLAT","ITGA2","NKAIN4","OLFML3","SOX9","PMEPA1","SPP1","IGFBP7","COL18A1","CSRNP2","UBTD1","ACTR1B","AC022075.3","AL392172.1","EIF5A2")

#' check the mk gene expression in LW60 & LW61
data.temp <- data.list$it
temp.sel.cell <- EBD2.meta %>% filter(cellType %in% c("LW61","LW60")) %>% filter(EML %in% c("EB_ELC","EB_HLC","EB_TLC","EB_UT9")) %>% pull(cell)
temp.exp <- data.temp@assays$RNA@data[unlist(temp.mk.gene),temp.sel.cell] %>% as.data.frame()
temp.sel.exp <- t(apply(temp.exp,1,scale))
colnames(temp.sel.exp) <- colnames(temp.exp)
rownames(temp.sel.exp) <- rownames(temp.exp)
temp.sel.exp[temp.sel.exp>3.5] <- 3.5
temp.sel.exp[temp.sel.exp<  -3.5] <- -3.5
temp.anno <- EBD2.meta %>% filter(cell %in% temp.sel.cell ) %>% select(cell,EML,cellType,subCT)%>% arrange(EML,subCT,cellType) %>% tibble::column_to_rownames("cell") 
#+ fig.width=9,fig.height=9
pheatmap(temp.sel.exp[,rownames(temp.anno)],cluster_rows=F,cluster_cols=F,scale="none",annotation_col=temp.anno,show_colnames=F,color=colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100), border_color="NA",fontsize_row = 5,main="Top markers",gaps_row = unlist(lapply(temp.mk.gene,function(x){return(length(x))})) %>% cumsum())

#' check integrated cluster on the local object LW60 
data.temp=data.list$LW60
temp.plot <- list()
temp.plot$EML <- DimPlot(data.temp ,label=T)+NoAxes()+NoLegend()
#for (n in unique(EBD2.meta$subCT)) {
for (n in c("EB_UN34","EB_UE7","EB_ELC","EB_HLC", "EB_TLC0","EB_TLC4","EB_TLC8","EB_TLC9","EB_U10","EB_U5","EB_U6")) {
  temp.plot[[n]] <- DimPlot(data.temp,cells.highlight=(EBD2.meta %>% filter(subCT==n) %>% pull(cell)))+NoAxes()+NoLegend()+ggtitle(n)+FunTitle()
}
#+ fig.width=9,fig.height=9
print(cowplot::plot_grid(plotlist=temp.plot))



#' downsample dataset (chose the cells conserved in integration and local cluster resutls)
temp.cell <- list()
for (n in c("EB_ELC","EB_HLC", "EB_TLC0","EB_TLC8","EB_U10","EB_U5","EB_U6")) {
  temp.cell[[n]] <- EBD2.LW60.meta %>% filter(subCT==n ) %>% filter(cell %in% (EBD2.meta %>% filter(subCT==n) %>% pull(cell)))
}
temp.cell$EB_TLC4 <- EBD2.LW60.meta %>% filter(subCT=="EB_TLC49" ) %>% filter(cell %in% (EBD2.meta %>% filter(subCT=="EB_TLC4") %>% pull(cell)))
temp.cell$EB_TLC9 <- EBD2.LW60.meta %>% filter(subCT=="EB_TLC49" ) %>% filter(cell %in% (EBD2.meta %>% filter(subCT=="EB_TLC9") %>% pull(cell)))

EBD2.LW60.meta.ds <- do.call("bind_rows",lapply(temp.cell[c("EB_ELC","EB_HLC")],function(x){return(FunMaSF(x,150))})) %>% bind_rows(
  do.call("bind_rows",lapply(temp.cell[c("EB_U10","EB_U5","EB_U6")],function(x){return(FunMaSF(x,100))}))
) %>% bind_rows(do.call("bind_rows",lapply(temp.cell[c("EB_TLC0","EB_TLC8","EB_TLC4","EB_TLC9")],function(x){return(FunMaSF(x,50))})) %>% FunMaSF(150) )

DimPlot(data.list$LW60,cells.highlight=EBD2.LW60.meta.ds$cell)+NoAxes()+NoLegend()+ggtitle("Downsampled cells in LW60")+FunTitle()



#' save object 
if (!file.exists(paste0("tmp_data/",TD,"/meta.EBD2.LW60.further.rds")) | rewrite) {
  print("save output")
  saveRDS(EBD2.meta,paste0("tmp_data/",TD,"/meta.EBD2.further.rds"))
  saveRDS(EBD2.LW60.meta,paste0("tmp_data/",TD,"/meta.EBD2.LW60.further.rds"))
  saveRDS(EBD2.LW60.meta.ds,file=paste0("tmp_data/",TD,"/meta.EBD2.LW60.further.ds.rds"))
}

  
