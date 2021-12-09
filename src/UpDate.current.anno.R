#' title: "UpDate published annotation"

# R4.0
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
suppressMessages(library(SeuratWrappers))
suppressMessages(library(readxl))
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

rename <- dplyr::rename
select<- dplyr::select
filter <- dplyr::filter

# loading meta data
ref.pj <-  c("D3post","CS7","SPH2016")
meta.filter<- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% mutate(EML=ifelse(EML %in% c("Ectoderm"),"Amnion",EML))

# loading whole integration information
data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/Pub.D2.umap.cord.rds"))
#data.ob <- readRDS(paste0("tmp_data/",TD,"/Pub.D2.mnn.rds"))
#data.ob@meta.data$seurat_clusters=(data.ob.umap %>% tibble::column_to_rownames("cell"))[colnames(data.ob),"seurat_clusters"]
#DimPlot(data.ob,group.by="seurat_clusters",label=T)+NoAxes()+NoLegend()


#' fix the Xiang et al "ICM" issue
meta.filter <- meta.filter%>% mutate(EML=ifelse(EML=="ICM" & pj=="D3post","TE",EML))

# New annotation from Meistermann_et_al for SPH2016
temp.anno <- read_xlsx(("data/Meistermann_et_al/1-s2.0-S1934590921001855-mmc2.xlsx"), sheet = 1) %>% tbl_df() %>% filter(Dataset=="Petropoulos2016") %>% rename(cell=Cell.Name,EML=UMAP.cluster) %>% select(cell,EML)
temp.anno <- temp.anno %>% mutate(EML=ifelse(EML %in% c("B1_B2","B1.EPI"),"EarlyBlastocyst",EML))%>% mutate(EML=ifelse(EML %in% c("EPI"),"EPI",EML)) %>% mutate(EML=ifelse(EML %in% c("EightCells","Morula"),"8C_Morula",EML)) %>% mutate(EML=ifelse(EML %in% c("EPI.PrE"),"EPI.PrE.INT",EML)) %>% mutate(EML=ifelse(EML %in% c("early_TE","late_TE","medium_TE",  "EPI.early_TE","EPI.PrE.TE","PrE.TE"),"TE",EML)) %>% mutate(EML=ifelse(EML %in% c("PrE"),"PE",EML))
meta.filter <- rows_update(meta.filter,temp.anno %>% filter(cell %in% meta.filter$cell),by="cell") %>% mutate(EML=ifelse(EML %in% c("EM_NA") & pj=="SPH2016","8C_Morula",EML))%>% mutate(EML=ifelse(EML %in% c("INT") & pj=="SPH2016","EPI",EML))
## New annotation from StirparotEtAl et al., for SPH2016  annotation
#temp.anno <- read_xlsx(("data/stirparo2018_tableS4.xlsx"), sheet = 1) %>% tbl_df() %>% filter(Study == "Petropoulos et al., 2016 (ERP012552)")%>% select("Cell","Revised lineage (this study)")
#colnames(temp.anno) <- c("cell","EML")
#temp.anno <- temp.anno %>% mutate(cell=gsub("_", ".", cell)) %>% mutate(EML=gsub("undefined","EM_NA",EML))%>% mutate(EML=gsub("epiblast","EPI",EML))%>% mutate(EML=gsub("Inner cell mass","ICM",EML))%>% mutate(EML=gsub("intermediate","INT",EML))%>% mutate(EML=gsub("primitive_endoderm","PE",EML))%>% mutate(EML=gsub("trophectoderm","TE",EML))
#meta.filter <- rows_update(meta.filter,temp.anno %>% filter(cell %in% meta.filter$cell),by="cell") %>% mutate(EML=gsub("EM_NA","EM_Early",EML)) %>% mutate(EML=gsub("INT","ICM_TE",EML))
#meta.filter <- rows_update(meta.filter,meta.filter %>% filter(pj %in% "SPH2016" & EML=="EM_Early" & devTime %in% c("E5","E6")) %>% mutate(EML="EarlyBlastocyst") %>% select(cell,EML),by="cell")

#' original sophie's dataset
#temp.anno <-  read.delim("~/My_project/ScRNA_analysis/data/Published/E-MTAB-3929.short.txt",stringsAsFactors = F,row.names = 1)%>% tibble::rownames_to_column("cell") %>% tbl_df()  %>% setNames(c("cell","Embryo","Days","Lineage1","Lineage2")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="primitive endoderm", "PE")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="epiblast", "Epiblast")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="trophectoderm", "TE")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="not applicable", "Prelineage")) %>% mutate(EML=Lineage1)%>% select(cell,EML)
#meta.filter %>% filter(pj=="SPH2016") %>% left_join(read.delim("~/My_project/ScRNA_analysis/data/Published/E-MTAB-3929.short.txt",stringsAsFactors = F,row.names = 1)%>% tibble::rownames_to_column("cell") %>% tbl_df()  %>% setNames(c("cell","Embryo","Days","Lineage1","Lineage2")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="primitive endoderm", "PE")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="epiblast", "EPI")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="trophectoderm", "TE")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="not applicable", "Prelineage")) %>% mutate(cp.anno=Lineage1)%>% select(cell,cp.anno),by="cell")  %>% left_join(read_xlsx(("data/Meistermann_et_al/1-s2.0-S1934590921001855-mmc2.xlsx"), sheet = 1) %>% tbl_df()  %>% rename(cell=Cell.Name) %>% select(cell,StirparotEtAl.lineage,Cell.state,UMAP.cluster),by="cell") %>% rename(JP.anno=EML) %>% group_by(cp.anno,JP.anno,StirparotEtAl.lineage,UMAP.cluster) %>% summarise(nCell=n_distinct(cell)) %>% arrange(cp.anno,UMAP.cluster,JP.anno) %>% write.table("tmp_data/SPH.anno.compair.tsv",col.names = T,row.names = F,quote=F,sep="\t")


#' fix the psa-epi cells issue from D3post 
temp.anno <- data.ob.umap %>% filter(pj=="D3post" & EML=='PSA-EPI') %>% mutate(EML="Epi")  %>% mutate(EML=ifelse(seurat_clusters %in% c("C3"), "TE",EML))%>% mutate(EML=ifelse(seurat_clusters %in% c("C6"), "Mes",EML))%>% mutate(EML=ifelse(seurat_clusters %in% c("C5"), "Endoderm",EML)) %>% select(cell,EML)
meta.filter <- rows_update(meta.filter,temp.anno %>% filter(cell %in% meta.filter$cell),by="cell") 
#' fix the mesdoerm cells from D3post 
temp.anno <- data.ob.umap %>% filter(pj=="D3post") %>%  filter(seurat_clusters %in% c("C4","C6")) %>% select(cell) %>% mutate(EML="Mes")
meta.filter <- rows_update(meta.filter,temp.anno,by="cell") 

#'  update the CS7 annotation 
CS7.meta <- readRDS(file=paste0("tmp_data/",TD,"/CS7.updata.anno.rds"))
meta.filter <- rows_update(meta.filter,CS7.meta  %>% select(cell,EML) ,by="cell") 


#' output reference meta data
ref.meta  <- meta.filter %>% filter(pj %in% ref.pj) %>% mutate(EML=ifelse(EML %in% c("AdvMes","AxMes","EmMes","NasMes"),"Mes",EML)) %>% mutate(EML=ifelse(EML %in% c("EPI","Epiblast"),"Epi",EML)) %>% mutate(EML=ifelse(EML %in% c("CTB","EarlyTE","EVT","LateTE","STB"),"TE",EML))%>% mutate(EML=ifelse(EML %in% c("Endoderm","PE","PrE"),"Endoderm",EML)) %>% rename(ref_cell=cell,ref_EML=EML,ref_pj=pj)

if (!file.exists(paste0("tmp_data/",TD,"/ref.meta.rds")) | rewrite) {
  print("save output")
  saveRDS(ref.meta,file=paste0("tmp_data/",TD,"/ref.meta.rds"))
}


