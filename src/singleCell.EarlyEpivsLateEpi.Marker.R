#' ---
#' title: "find the markers to distinct TE and Amnion"
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


# number of conserved markers
if (file.exists(paste0("tmp_data/",TD,"/human.EarlyLateEpi.Rdata")))  {
  load(paste0("tmp_data/",TD,"/human.EarlyLateEpi.Rdata"),verbose=T)
}else{
  #' loading data
  load("tmp_data/gene.meta.Rdata",verbose=T)
  counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))
  lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/",MM,".lognormExp.mBN.rds"))
  #' whole dataset integration results
  data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/Pub.D2.umap.cord.rds")) 
  
  meta.filter<- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% filter(pj %in% c("EBD2","IBD2","JPF2019","D3post","CS7","SPH2016","Blakeley","nBGuo","nicolBla") ) %>% mutate(EML=ifelse(EML %in% c("Ectoderm"),"Amnion",EML))
  
  #' loading updated reference meta data
  ref.meta <- readRDS(paste0("tmp_data/",TD,"/ref.meta.rds"))
  
  #' updatae meta data
  meta.filter <- rows_update(meta.filter, ref.meta %>% mutate(cell=ref_cell,EML=ref_EML) %>% select(cell,EML),by="cell")
  temp.M <- meta.filter 
  temp.sel.expG <-rownames(lognormExp.mBN )
  
  data.merge <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) 
  data.merge@assays$RNA@data <- as.matrix(lognormExp.mBN[temp.sel.expG,colnames(data.merge)])
  
  
  cells.list <- list()
  
  cells.list$Blakeley_Epi <- meta.filter %>% filter(pj=="Blakeley" & EML=="EPI")  %>% pull(cell)
  cells.list$SPH2016_Epi <- meta.filter %>% filter(pj=="SPH2016" & EML=="Epi" & devTime !="E5")  %>% pull(cell)
  cells.list$D3post_Epi_Early <- meta.filter %>% filter(pj=="D3post" & EML %in% c("Epi") & devTime %in% c("E6","E7"))  %>% pull(cell)#"E8",
  cells.list$D3post_Epi_Late <- meta.filter %>% filter(pj=="D3post" & EML %in% c("Epi") & devTime %in% c("E12","E14"))  %>% pull(cell) %>% intersect(data.ob.umap %>% filter(seurat_clusters=="C9") %>% pull(cell))
  cells.list$CS7_Epi <- meta.filter %>% filter(pj=="CS7" & EML %in% c("Epi"))  %>% pull(cell)
  cells.list$nBGuo_Epi <- meta.filter %>% filter(pj=="nBGuo" & EML %in% c("Epiblast"))  %>% pull(cell)
  
  # cells.list$IBD2_ELC <-  meta.filter %>% filter(pj=="IBD2" & EML %in% c("IB_Epi")) %>% pull(cell)
  # cells.list$EBD2_ELC <-  meta.filter %>% filter(pj=="EBD2" & EML %in% c("EB_ELC")) %>% pull(cell)
  cells.list$IBD2_ELC.C0 <-  data.ob.umap %>% filter(pj=="IBD2" & EML %in% c("IB_Epi") & seurat_clusters %in% c("C0")) %>% pull(cell)
  cells.list$IBD2_ELC.nC0 <-  data.ob.umap %>% filter(pj=="IBD2" & EML %in% c("IB_Epi") & (!seurat_clusters %in% c("C0"))) %>% pull(cell)
  cells.list$EBD2_ELC.C0 <-  data.ob.umap  %>% filter(pj=="EBD2" & EML %in% c("EB_ELC")  & seurat_clusters %in% c("C0") ) %>% pull(cell)
  cells.list$EBD2_ELC.nC0 <-  data.ob.umap  %>% filter(pj=="EBD2" & EML %in% c("EB_ELC")  & (!seurat_clusters %in% c("C0") )) %>% pull(cell)
  
  cells.list$nBGuo_ELC <-  meta.filter %>% filter(pj=="nBGuo" & EML %in% c("ELC")) %>% pull(cell)
  cells.list$nicolBla_ELC <-  meta.filter %>% filter(pj=="nicolBla" & EML %in% c("nicolBla_ELC")) %>% pull(cell)
  
  lapply(cells.list,length)
  
  hm.DE.list <- list()
  hm.DE.list$SPH2016_Epi <- list()
  hm.DE.list$SPH2016_Epi$D3post_Epi_Late <- 1
  hm.DE.list$SPH2016_Epi$CS7_Epi <- 1
  
  hm.DE.list$D3post_Epi_Early <- list()
  hm.DE.list$D3post_Epi_Early$D3post_Epi_Late <- 1
  hm.DE.list$D3post_Epi_Early$CS7_Epi <- 1
  
  hm.EarlyvsLate.Epi.out <- list()
  for (n in names(hm.DE.list)) {
    for (m in names(hm.DE.list[[n]])) {
      hm.EarlyvsLate.Epi.out[[n]][[m]] <- FunDEG(data.merge,cells.list,n,m) %>% mutate(SampleA=n,SampleB=m)
    }
  }
  
  save(hm.EarlyvsLate.Epi.out,file=paste0("tmp_data/",TD,"/human.EarlyLateEpi.Rdata"))
}

nC.cutoff <- length(hm.DE.list$SPH2016_Epi)

EarlyVSLate.Epi.mk <-  hm.EarlyvsLate.Epi.out$SPH2016_Epi %>% do.call("bind_rows",.) %>% mutate(Type=ifelse(power > power.cutoff & avg_log2FC > log2FC.cutoff & pct.1 > pct.max.cutoff & pct.2 < pct.min.cutoff, "UpDEG",ifelse(power > power.cutoff & avg_log2FC < -log2FC.cutoff & pct.1 < pct.min.cutoff & pct.2 > pct.max.cutoff, "DownDEG","NotDEG"))) %>% filter(Type !="NotDEG")%>% group_by(gene,Type) %>% summarise(nC=n()) %>% filter(nC >=nC.cutoff) %>% mutate(compair="em") %>% select(-nC) %>% ungroup()%>% select(-compair) %>% mutate(Type=gsub("DownDEG","LateEpi",Type)) %>% mutate(Type=gsub("UpDEG","EarlyEpi",Type)) %>% arrange(Type) %>% semi_join( hm.EarlyvsLate.Epi.out$D3post_Epi_Early %>% do.call("bind_rows",.) %>% mutate(Type=ifelse(power > power.cutoff & avg_log2FC > log2FC.cutoff & pct.1 > pct.max.cutoff & pct.2 < pct.min.cutoff, "UpDEG",ifelse(power > power.cutoff & avg_log2FC < -log2FC.cutoff & pct.1 < pct.min.cutoff & pct.2 > pct.max.cutoff, "DownDEG","NotDEG"))) %>% filter(Type !="NotDEG")%>% group_by(gene,Type) %>% summarise(nC=n()) %>% filter(nC >=nC.cutoff) %>% mutate(compair="em") %>% select(-nC) %>% ungroup()%>% select(-compair) %>% mutate(Type=gsub("DownDEG","LateEpi",Type)) %>% mutate(Type=gsub("UpDEG","EarlyEpi",Type)) %>% arrange(Type),by="gene") %>% left_join(Gene.anno,by="gene")

saveRDS(EarlyVSLate.Epi.mk,file=paste0("tmp_data/",TD,"/human.EarlyVSLate.Epi.sigMarker.rds"))

#' check the genes expression
#' ### VlnPlot for key marker genes
temp.sel.gene <- EarlyVSLate.Epi.mk$gene
temp.M <- as.list(names(cells.list)) %>% lapply(function(x) {data.frame(cell=cells.list[[x]],anno=x)}) %>% do.call("bind_rows",.) %>% tbl_df()

temp.sel.exp <- lognormExp.mBN[temp.sel.gene ,temp.M$cell]%>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,logExp,-Gene) %>% inner_join(temp.M,by="cell") %>% mutate(od=factor(anno,c("Blakeley_Epi","nBGuo_Epi" ,"SPH2016_Epi","D3post_Epi_Early","D3post_Epi_Late","CS7_Epi","IBD2_ELC.C0","IBD2_ELC.nC0","EBD2_ELC.C0","EBD2_ELC.nC0","nBGuo_ELC","nicolBla_ELC"),ordered = T)) %>% arrange(od)  

temp.plot <- list()
for (g in temp.sel.gene ) {
  temp.plot[[paste0("MK.Vln.",g)]] <- temp.sel.exp%>% filter(Gene==g) %>% ggplot(mapping=aes(x=od,y=logExp)) + geom_violin(mapping=aes(fill=anno),color="black",scale = "width")+geom_jitter(size=0.05)+ theme_classic() + theme(axis.text.x=element_text(angle = 40,hjust = 1)) +NoLegend()+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="bold"))+ylab("Expression level")+xlab("")+ylab("")
}
cowplot::plot_grid(plotlist = temp.plot)
