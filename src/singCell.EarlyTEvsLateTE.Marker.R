#' ---
#' title: "find the markers to distinct Early TE and Late TE"
#' output: 
#'  html_document:
#'    code_folding: hide
#'    
#' ---

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
log2FC.cutoff <- 0.25
pct.min.cutoff <- c(0.25)
pct.max.cutoff <- c(0.5)

#' loading Gene annotation
Gene.anno <- read.delim("~/Genome/Human/RefSeq/Homo_sapiens.GRCh38.gene.description",sep="\t",stringsAsFactors=F,header = F)  %>% tbl_df() %>% rename(gene=V1,annotaton=V2)

counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/",MM,".lognormExp.mBN.rds"))
#' whole dataset integration results
data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/Pub.D2.umap.cord.rds")) 

meta.filter<- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% filter(pj %in% c("EBD2","IBD2","JPF2019","D3post","CS7","SPH2016","Blakeley","nBGuo","nicolBla") ) %>% mutate(EML=ifelse(EML %in% c("Ectoderm"),"Amnion",EML))

#' loading updated reference meta data
ref.meta <- readRDS(paste0("tmp_data/",TD,"/ref.meta.rds"))

#' updatae meta data
meta.filter <- rows_update(meta.filter, ref.meta %>% filter(ref_pj=="SPH2016") %>% mutate(cell=ref_cell,EML=ref_EML) %>% select(cell,EML),by="cell")


zs.limit <- 2.5
heat.col <- colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100)
#' loading top DEG between mural and polar TE 
mvp.DE.tsv <- read.delim("doc/mural.vs.polar.txt",stringsAsFactors = F,head=F) %>% tbl_df() %>% filter(V1 %in% rownames(lognormExp.mBN)) %>% mutate(Type=ifelse(V3 >0, "polar","mural")) %>% group_by(Type) %>% top_n(100,-V2) %>% arrange(Type)

for (n in c("EBD2","nicolBla","nBGuo")) {
  temp.M <-  data.ob.umap %>% filter(pj==n & rename_EML=="TLC" & cluster_EML=="TLC")
  temp.anno <- temp.M %>% select(cell,seurat_clusters,pj) %>% tibble::column_to_rownames("cell")
  temp.sel.cell <- temp.M$cell
  temp.sel.gene <- mvp.DE.tsv$V1
  row.gaps <-   100
  temp.exp <- lognormExp.mBN[temp.sel.gene,temp.sel.cell] 
  temp.sel.exp <- t(apply(temp.exp,1,scale))
  colnames(temp.sel.exp) <- colnames(temp.exp)
  rownames(temp.sel.exp) <- rownames(temp.exp)
  temp.sel.exp[temp.sel.exp>zs.limit] <- zs.limit
  temp.sel.exp[temp.sel.exp<  (-1*zs.limit)] <- -1*zs.limit
  gaps <- 10
  pheatmap(temp.sel.exp,cluster_rows=F,,cluster_cols=T,scale="none",show_colnames=F,show_rownames=F, fontsize_row=7,border_color="NA",main=paste("pure TLC of", n),gaps_row = rep(row.gaps,each=2),col=heat.col,annotation_col = temp.anno,clustering_method = "ward.D2")
  
}

n="IBD2"
temp.M <-  data.ob.umap %>% filter(pj==n & EML=="IB_TE")
temp.anno <- temp.M %>% select(cell,seurat_clusters,pj) %>% tibble::column_to_rownames("cell")
temp.sel.cell <- temp.M$cell
temp.sel.gene <- mvp.DE.tsv$V1
row.gaps <-   100
temp.exp <- lognormExp.mBN[temp.sel.gene,temp.sel.cell] 
temp.sel.exp <- t(apply(temp.exp,1,scale))
colnames(temp.sel.exp) <- colnames(temp.exp)
rownames(temp.sel.exp) <- rownames(temp.exp)
temp.sel.exp[temp.sel.exp>zs.limit] <- zs.limit
temp.sel.exp[temp.sel.exp<  (-1*zs.limit)] <- -1*zs.limit
gaps <- 10
pheatmap(temp.sel.exp,cluster_rows=F,,cluster_cols=T,scale="none",show_colnames=F,show_rownames=F, fontsize_row=7,border_color="NA",main=paste("previous TLC of", n),gaps_row = rep(row.gaps,each=2),col=heat.col,annotation_col = temp.anno,clustering_method = "ward.D2")

for (n in "D3post","SPH2016","nBGuo") {
  temp.M <-  data.ob.umap %>% filter(pj==n & EML %in% c("TE","CTB","EVT","STB","EarlyTE" ))
  temp.anno <- temp.M %>% select(cell,EML,devTime,pj) %>% tibble::column_to_rownames("cell")
  temp.sel.cell <- temp.M$cell
  temp.sel.gene <- mvp.DE.tsv$V1
  row.gaps <-   100
  temp.exp <- lognormExp.mBN[temp.sel.gene,temp.sel.cell] 
  temp.sel.exp <- t(apply(temp.exp,1,scale))
  colnames(temp.sel.exp) <- colnames(temp.exp)
  rownames(temp.sel.exp) <- rownames(temp.exp)
  temp.sel.exp[temp.sel.exp>zs.limit] <- zs.limit
  temp.sel.exp[temp.sel.exp<  (-1*zs.limit)] <- -1*zs.limit
  gaps <- 10
  pheatmap(temp.sel.exp,cluster_rows=F,,cluster_cols=T,scale="none",show_colnames=F,show_rownames=F, fontsize_row=7,border_color="NA",main=paste("TE of ",n ),gaps_row = rep(row.gaps,each=2),col=heat.col,annotation_col = temp.anno,clustering_method = "ward.D2")
  
}

# number of conserved markers
if (file.exists(paste0("tmp_data/",TD,"/human.EarlyLateTE.Rdata")))  {
  load(paste0("tmp_data/",TD,"/human.EarlyLateTE.Rdata"),verbose=T)
}else{
  #' loading data
  load("tmp_data/gene.meta.Rdata",verbose=T)
  
  temp.M <- meta.filter 
  temp.sel.expG <-rownames(lognormExp.mBN )
  
  data.merge <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) 
  data.merge@assays$RNA@data <- as.matrix(lognormExp.mBN[temp.sel.expG,colnames(data.merge)])
  
  
  cells.list <- list()
  
  cells.list$Blakeley_TE <- meta.filter %>% filter(pj=="Blakeley" & EML=="TE")  %>% pull(cell)
  cells.list$SPH2016_TE <- meta.filter %>% filter(pj=="SPH2016" & EML=="TE" & devTime !="E5")  %>% pull(cell)
  cells.list$D3post_TE_Early <- meta.filter %>% filter(pj=="D3post" & EML %in% c("CTB","ICM") & devTime %in% c("E6","E7"))  %>% pull(cell)#"E8",
  cells.list$D3post_TE_LateCTB <- meta.filter %>% filter(pj=="D3post" & EML %in% c("CTB") & devTime %in% c("E10","E12","E14"))  %>% pull(cell) 
  cells.list$D3post_TE_LateEVT <- meta.filter %>% filter(pj=="D3post" & EML %in% c("EVT") & devTime %in% c("E10","E12","E14"))  %>% pull(cell) 
  cells.list$D3post_TE_LateSTB <- meta.filter %>% filter(pj=="D3post" & EML %in% c("STB") & devTime %in% c("E10","E12","E14"))  %>% pull(cell) 
  cells.list$nBGuo_TE <- meta.filter %>% filter(pj=="nBGuo" & EML %in% c("TE","EarlyTE"))  %>% pull(cell)
  
  # cells.list$IBD2_TLC <-  meta.filter %>% filter(pj=="IBD2" & EML %in% c("IB_TE")) %>% pull(cell)
  cells.list$EBD2_TLC <-  meta.filter %>% filter(pj=="EBD2" & EML %in% c("EB_TLC")) %>% pull(cell)
  #cells.list$EBD2_TLC.C12 <-  data.ob.umap %>% filter(pj=="EBD2" & EML %in% c("EB_TLC") & seurat_clusters %in% c("C12")) %>% pull(cell)
  # cells.list$EBD2_TLC.nC12 <-  data.ob.umap %>% filter(pj=="EBD2" & EML %in% c("EB_TLC") & (!seurat_clusters %in% c("C12"))) %>% pull(cell)
  cells.list$nBGuo_TLC <-  meta.filter %>% filter(pj=="nBGuo" & EML %in% c("TLC")) %>% pull(cell)
  cells.list$nicolBla_TLC <-  meta.filter %>% filter(pj=="nicolBla" & EML %in% c("nicolBla_TLC")) %>% pull(cell)
  
  lapply(cells.list,length)
  
  hm.DE.list <- list()
  hm.DE.list$D3post_TE_LateCTB <- list()
  hm.DE.list$D3post_TE_LateCTB$SPH2016_TE <- 1
  hm.DE.list$D3post_TE_LateCTB$D3post_TE_Early<- 1
  
  hm.DE.list$D3post_TE_LateSTB <- list()
  hm.DE.list$D3post_TE_LateSTB$SPH2016_TE<- 1
  hm.DE.list$D3post_TE_LateSTB$D3post_TE_Early<- 1
  
  hm.DE.list$D3post_TE_LateEVT <- list()
  hm.DE.list$D3post_TE_LateEVT$SPH2016_TE<- 1
  hm.DE.list$D3post_TE_LateEVT$D3post_TE_Early<- 1
  
  
  hm.EarlyvsLate.TE.out <- list()
  for (n in names(hm.DE.list)) {
    for (m in names(hm.DE.list[[n]])) {
      hm.EarlyvsLate.TE.out[[n]][[m]] <- FunDEG(data.merge,cells.list,n,m) %>% mutate(SampleA=n,SampleB=m)
    }
  }
  
  save(hm.EarlyvsLate.TE.out,file=paste0("tmp_data/",TD,"/human.EarlyLateTE.Rdata"))
}

nC.cutoff <- length(hm.DE.list$D3post_TE_LateCTB)
LateTE.mk <- list()
LateTE.mk$CTB <-   hm.EarlyvsLate.TE.out$D3post_TE_LateCTB %>% do.call("bind_rows",.) %>% mutate(Type=ifelse(power > power.cutoff & avg_log2FC > log2FC.cutoff & pct.1 > pct.max.cutoff & pct.2 < pct.min.cutoff, "UpDEG",ifelse(power > power.cutoff & avg_log2FC < -log2FC.cutoff & pct.1 < pct.min.cutoff & pct.2 > pct.max.cutoff, "DownDEG","NotDEG"))) %>% filter(Type !="NotDEG")%>% group_by(gene,Type) %>% summarise(nC=n(),power=mean(power),avg_log2FC=mean(avg_log2FC)) %>% filter(nC >=nC.cutoff)%>% select(-nC) %>% ungroup() %>% mutate(Type=gsub("DownDEG","EarlyTE",Type)) %>% mutate(Type=gsub("UpDEG","CTB",Type)) %>% arrange(Type) %>% filter(Type=="CTB") %>% arrange(desc(power)) %>% head(9)
LateTE.mk$STB <-    hm.EarlyvsLate.TE.out$D3post_TE_LateSTB %>% do.call("bind_rows",.) %>% mutate(Type=ifelse(power > power.cutoff & avg_log2FC > log2FC.cutoff & pct.1 > pct.max.cutoff & pct.2 < pct.min.cutoff, "UpDEG",ifelse(power > power.cutoff & avg_log2FC < -log2FC.cutoff & pct.1 < pct.min.cutoff & pct.2 > pct.max.cutoff, "DownDEG","NotDEG"))) %>% filter(Type !="NotDEG")%>% group_by(gene,Type) %>% summarise(nC=n(),power=mean(power),avg_log2FC=mean(avg_log2FC)) %>% filter(nC >=nC.cutoff)%>% select(-nC) %>% ungroup() %>% mutate(Type=gsub("DownDEG","EarlyTE",Type)) %>% mutate(Type=gsub("UpDEG","STB",Type)) %>% arrange(Type) %>% filter(Type=="STB")%>% arrange(desc(power))%>% head(9)

LateTE.mk$EVT <-    hm.EarlyvsLate.TE.out$D3post_TE_LateEVT %>% do.call("bind_rows",.) %>% mutate(Type=ifelse(power > power.cutoff & avg_log2FC > log2FC.cutoff & pct.1 > pct.max.cutoff & pct.2 < pct.min.cutoff, "UpDEG",ifelse(power > power.cutoff & avg_log2FC < -log2FC.cutoff & pct.1 < pct.min.cutoff & pct.2 > pct.max.cutoff, "DownDEG","NotDEG"))) %>% filter(Type !="NotDEG")%>% group_by(gene,Type) %>% summarise(nC=n(),power=mean(power),avg_log2FC=mean(avg_log2FC)) %>% filter(nC >=nC.cutoff)%>% select(-nC) %>% ungroup() %>% mutate(Type=gsub("DownDEG","EarlyTE",Type)) %>% mutate(Type=gsub("UpDEG","EVT",Type)) %>% arrange(Type) %>% filter(Type=="EVT")%>% arrange(desc(power))%>% head(9)


saveRDS(LateTE.mk,file=paste0("tmp_data/",TD,"/human.EarlyVSLate.TE.sigMarker.rds"))

#' check the genes expression
#' ### VlnPlot for key marker genes
for (s in names(LateTE.mk)) {
  temp.sel.gene <- LateTE.mk[[s]]$gene
  temp.M <- as.list(names(cells.list)) %>% lapply(function(x) {data.frame(cell=cells.list[[x]],anno=x)}) %>% do.call("bind_rows",.) %>% tbl_df()
  
  temp.sel.exp <- lognormExp.mBN[temp.sel.gene ,temp.M$cell]%>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,logExp,-Gene) %>% inner_join(temp.M,by="cell") %>% mutate(od=factor(anno,c("Blakeley_TE","nBGuo_TE" ,"SPH2016_TE","D3post_TE_Early","D3post_TE_LateCTB","D3post_TE_LateSTB","D3post_TE_LateEVT","EBD2_TLC","nBGuo_TLC","nicolBla_TLC"),ordered = T)) %>% arrange(od)  
  
  temp.plot <- list()
  for (g in temp.sel.gene ) {
    temp.plot[[paste0("MK.Vln.",g)]] <- temp.sel.exp%>% filter(Gene==g) %>% ggplot(mapping=aes(x=od,y=logExp)) + geom_violin(mapping=aes(fill=anno),color="black",scale = "width")+geom_jitter(size=0.05)+ theme_classic() + theme(axis.text.x=element_text(angle = 40,hjust = 1)) +NoLegend()+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="bold"))+ylab("Expression level")+xlab("")+ylab("")
  }
  cowplot::plot_grid(plotlist = temp.plot)
  
}


temp.sel.gene <- LateTE.mk %>% lapply(function(x){head(x,5) %>% pull(gene)}) %>% unlist() %>% unique()
temp.M <- as.list(names(cells.list)) %>% lapply(function(x) {data.frame(cell=cells.list[[x]],anno=x)}) %>% do.call("bind_rows",.) %>% tbl_df()

temp.sel.exp <- lognormExp.mBN[temp.sel.gene ,temp.M$cell]%>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,logExp,-Gene) %>% inner_join(temp.M,by="cell") %>% mutate(od=factor(anno,c("Blakeley_TE","nBGuo_TE" ,"SPH2016_TE","D3post_TE_Early","D3post_TE_LateCTB","D3post_TE_LateSTB","D3post_TE_LateEVT","EBD2_TLC","nBGuo_TLC","nicolBla_TLC"),ordered = T)) %>% arrange(od)  

temp.plot <- list()
for (g in temp.sel.gene ) {
  temp.plot[[paste0("MK.Vln.",g)]] <- temp.sel.exp%>% filter(Gene==g) %>% ggplot(mapping=aes(x=od,y=logExp)) + geom_violin(mapping=aes(fill=anno),color="black",scale = "width")+geom_jitter(size=0.05)+ theme_classic() + theme(axis.text.x=element_text(angle = 40,hjust = 1)) +NoLegend()+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="bold"))+ylab("Expression level")+xlab("")+ylab("")
}
cowplot::plot_grid(plotlist = temp.plot)
