#' ---
#' title: "to plot"
#' output: 
#'  html_document:
#'    code_folding: hide
#'    
#' ---
#' 
#' 
#' R4.0
rm(list=ls())

condaENV <- "/home/chenzh/miniconda3/envs/R4.0"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(patchwork))
suppressMessages(library(cowplot))
suppressMessages(library(pheatmap))
suppressMessages(library(Seurat))
suppressMessages(library(ggalluvial))

# working directory
DIR <- "~/My_project/CheckBlastoids"
setwd(DIR)
knitr::opts_knit$set(root.dir=DIR)
#' Loading R functions
#source("~/PC/R_code/functions.R")
#source("~/PC/SnkM/SgCell.R")
source("src/local.quick.fun.R")
source("src/figures.setting.R")

#' loading font style
#extrafont::font_import(path="/home/chenzh/Downloads/fonts/",,prompt = F)
extrafont::loadfonts()
theme_set(theme_gray(base_family="Arial"))
#pdf.options(paper = "a4")
#theme_set(theme_gray(base_size = 6))
#+theme(text = element_text(size = 6))
#width = 210, height = 297, units = "mm"
#pdf("test.pdf",8,12)
#' loading options
TD="Nov14_2021"
MM="Pub"


#' loading data
load("tmp_data/gene.meta.Rdata",verbose=T)
meta.filter<- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% filter(pj %in% c("EBD2","IBD2","JPF2019","D3post","CS7","SPH2016","nBGuo","Blakeley") ) %>% mutate(EML=ifelse(EML %in% c("Ectoderm"),"Amnion",EML)) %>% mutate(EML=ifelse(EML=="EB_UI13","EB_ELC",EML))

#' loading the raw counts of CS7 and D3post
counts.filter <- (readRDS(paste0("tmp_data/",TD,"/counts.filter.rds")))[,meta.filter %>% filter(pj %in% c("D3post","CS7")) %>% pull(cell)]

#' loading normalized expression
lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/",MM,".lognormExp.mBN.rds"))

#' loading Amnion vs TE DEGs
hm.AMvsTE.out <- readRDS(paste0("tmp_data/",TD,"/human.AmnionVsTE.sigMarker.rds"))
#' loading all lineage markers
lineage.mk <- readRDS(paste0("tmp_data/",TD,"/",MM,".mk.allL.rds"))


#' loading whole human integration results
data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/Pub.D2.umap.cord.rds")) 


#' loading  human integration results except for JPF2019
data.ob.NoJPF2019.umap <- readRDS(paste0("tmp_data/",TD,"/Pub.D2.withoutAMLC.umap.cord.rds")) 

#' loading cross-species integration results
human.NHP.IT.umap.list <- readRDS(paste0("tmp_data/",TD,"/NHP.IT.umap.rds"))


#' loading NHP data
load(paste0("tmp_data/",TD,"/NHP.Rdata"),verbose=T)
NHP.lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/NHP.scran.norm.rds"))
NHP.lognormExp.mBN <- NHP.lognormExp.mBN[,NHP.meta %>% filter(EML %in% c("Amnion2","TE")) %>% pull(cell)]
#' loading psd-bulk expression data
load(paste0("tmp_data/",TD,"/Psd.norm.exp.Rdata"),verbose = T)

#'
ref.pj <-  c("D3post","CS7","SPH2016")

#' loading psd-bulk pca results
pca.out <- readRDS(paste0("tmp_data/",TD,"/psd.AmnVsTE.pca.out.rds"))


#' loading the one to check the CS7 re-annotation
load(paste0("tmp_data/",TD,"/check.CS7.sub.Rdata"),verbose=T)
#data.AEP.UMAP
#data.AEP.pure.DE

#' loading reference meta data
ref.meta <- readRDS(paste0("tmp_data/",TD,"/ref.meta.rds")) %>% select(ref_cell,ref_EML,ref_pj)

query_ref.raw <- readRDS(paste0("tmp_data/",TD,"/query_ref.raw.rds"))   %>% filter(query_pj %in% c("SPH2016","D3post","CS7","Blakeley","nBGuo","zhou2019"))%>% mutate(query_EML=FunForMatAnno(query_EML),ref_EML=FunForMatAnno(ref_EML)) %>% filter(! query_EML %in% c("ELC","EPI_Amnion","EPI.PrE.INT","EPI.PrE.INT","EPI_PriS","HLC","ICM","ICM-TE_trans","Undefined","TLC","ysTE","PriS_Amnion","Uncertain","8C_Morula","EarlyBlastocyst"))

#' loading prediction results
query_ref.anno.out <- readRDS(paste0("tmp_data/",TD,"/prediction.model.performance.rds"))

#query_ref.anno.out  %>% group_by(query_pj,ref_EML) %>% summarise(nC=n_distinct(query_cell)) %>% spread(ref_EML,nC) %>% replace(.,is.na(.),0)  %>% gather(ref_EML,nC,-query_pj) %>% mutate(Perp=nC/sum(nC))  %>% filter(ref_EML=="NoSigHits")  %>% filtertable(query_ref.anno.out  %>% filter(query_pj %in% c("Blakeley","CS7","D3post","EBD2","IBD2","JPF2019","nBGuo","nicolBla","zhou2019","SPH2016")) %>% filter(ref_EML=="NoSigHits") %>% mutate(cell=query_cell) %>% select(cell) %>% inner_join(meta.filter,by="cell") %>% pull(nGene) < 2500)

plot.results <- list()
plot.results.txt <- list()
plot.results.jpeg <- list()

unmatch.HLC.cells <- c("nBG_Gata3_Day3.1.pTEICM_8","nBG_Gata3_Day3.9.pTEICM_21","nBG_Gata3_Day3.1.pTEICM_16","nBG_Gata3_Day3.3.pTEICM_2","nBG_Gata3_Day3.7.pTEICM_22","nBG_Gata3_Day3.8.pTEICM_5","nBG_Gata3_Day4.6.mTE_2")
#' for the whole human dataset integration

# select the embryonic dataset
data.EM.umap.text.pos <- data.ob.umap  %>% filter(pj %in% c("D3post","CS7","SPH2016","nBGuo")) %>% filter(cellType=="EM") %>%group_by(rename_EML) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) %>% ungroup()  %>% filter(!rename_EML %in% c("3D_ICM","PSA-EPI","Unknown")) %>% mutate(UMAP_1=ifelse(rename_EML=="PriS",UMAP_1-1,UMAP_1)) #%>% mutate(rename_EML=gsub("ExE_Mes","Extraembryonic\nMesoderm",rename_EML))%>% mutate(rename_EML=gsub("PriS","Primitive\nStreak",rename_EML))

data.LC.umap.text.pos <- data.ob.umap  %>% filter(pj %in% c("EBD2","IBD2","JPF2019","nBGuo","nicolBla")) %>% filter(cellType!="EM") %>% filter(!cell %in% unmatch.HLC.cells) %>%group_by(rename_EML,pj) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) %>% ungroup()


plot.results$umap.em <- ggplot()+ geom_point(data.ob.umap %>% filter(pj %in% c("EBD2","IBD2","JPF2019","nBGuo","nicolBla"))%>% filter(cellType!="EM") ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.5)+geom_point( data.ob.umap %>% filter(pj %in% c("D3post","CS7","SPH2016","nBGuo")) %>% filter(cellType=="EM") %>% mutate(od=factor(rename_EML,c("ICM-TE_trans","ICM","Prelineage","TE","Mesoderm","Epiblast","Amnion","Endoderm","PriS","ExE_Mes","3D_ICM","PSA-EPI"),ordered = T)) %>% arrange(desc(od)) ,mapping=aes(x=UMAP_1,y=UMAP_2,fill=rename_EML,shape=pj),size=2,color="grey")+ggtitle("Embryonic dataset")+geom_text(data.EM.umap.text.pos ,mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")+theme_classic()+theme(plot.title = element_text(hjust=0.5))+scale_fill_manual(values=lineage.col.set)+scale_shape_manual(values = pj.shape.set,labels = c('SPH2016'="Petropoulos et al., 2016",'D3post'="Xiang et al., 2019",'CS7'="Tyser et al., 2020",'nBGuo'="Yanagida et al., 2021"))+NoAxes()+NoLegend()

plot.results$umap.em.legend <- ggplot()+ geom_point( data.ob.umap %>% filter(pj %in% c("D3post","CS7","SPH2016","nBGuo","nicolBla"))%>% filter(cellType=="EM") %>% filter(!rename_EML %in% c("PSA-EPI","3D_ICM","Unknown")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML,shape=pj),size=2)+theme_classic()+scale_color_manual("",values=lineage.col.set)+scale_shape_manual("",values = pj.shape.set.solid,labels = c('SPH2016'="Petropoulos et al., 2016",'D3post'="Xiang et al., 2020",'CS7'="Tyser et al., 2020",'nBGuo'="Yanagida et al., 2021")) + theme(legend.text = element_text( size = 6),legend.title = element_text( size = 1))
plot.results$umap.em.legend <- plot.results$umap.em.legend%>% ggpubr::get_legend() %>% ggpubr::as_ggplot()

umap.xlim=ggplot_build(plot.results$umap.em)$layout$panel_scales_x[[1]]$range$range
umap.ylim=ggplot_build(plot.results$umap.em)$layout$panel_scales_y[[1]]$range$range

plot.results.txt$umap.em <- ggplot()+ggtitle("Embryonic dataset")+geom_text(data.EM.umap.text.pos ,mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")+theme_classic()+theme(plot.title = element_text(hjust=0.5))+scale_fill_manual(values=lineage.col.set)+scale_shape_manual(values = pj.shape.set,labels = c('SPH2016'="Petropoulos et al., 2016",'D3post'="Xiang et al., 2019",'CS7'="Tyser et al., 2020",'nBGuo'="Yanagida et al., 2021"))+NoAxes()+NoLegend()+xlim(umap.xlim)+ylim(umap.ylim)
plot.results.jpeg$umap.em  <- ggplot()+ geom_point(data.ob.umap %>% filter(pj %in% c("EBD2","IBD2","JPF2019","nBGuo","nicolBla"))%>% filter(cellType!="EM") ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.5)+geom_point( data.ob.umap %>% filter(pj %in% c("D3post","CS7","SPH2016","nBGuo")) %>% filter(cellType=="EM") %>% mutate(od=factor(rename_EML,c("ICM-TE_trans","ICM","Prelineage","TE","Mesoderm","Epiblast","Amnion","Endoderm","PriS","ExE_Mes","3D_ICM","PSA-EPI"),ordered = T)) %>% arrange(desc(od)) ,mapping=aes(x=UMAP_1,y=UMAP_2,fill=rename_EML,shape=pj),size=2,color="grey")+ggtitle("")+theme_classic()+theme(plot.title = element_text(hjust=0.5))+scale_fill_manual(values=lineage.col.set)+scale_shape_manual(values = pj.shape.set,labels = c('SPH2016'="Petropoulos et al., 2016",'D3post'="Xiang et al., 2019",'CS7'="Tyser et al., 2020",'nBGuo'="Yanagida et al., 2021"))+NoAxes()+NoLegend()+xlim(umap.xlim)+ylim(umap.ylim)


plot.results$umap.JPF <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("JPF2019")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("JPF2019")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.5)+ggtitle("Zheng et al., 2019")+geom_text(data.EM.umap.text.pos %>% filter(rename_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3.5,color="Black")+  geom_text(data.LC.umap.text.pos %>% filter(pj=="JPF2019"),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")  +theme_classic()+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=lineage.col.set)+xlim(umap.xlim)+ylim(umap.ylim)
plot.results.txt$umap.JPF<-  ggplot()+ggtitle("Zheng et al., 2019")+geom_text(data.EM.umap.text.pos %>% filter(rename_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3.5,color="Black")+  geom_text(data.LC.umap.text.pos %>% filter(pj=="JPF2019"),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")  +theme_classic()+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=lineage.col.set)+xlim(umap.xlim)+ylim(umap.ylim)
plot.results.jpeg$umap.JPF <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("JPF2019")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("JPF2019")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.5)+ggtitle("") +theme_classic()+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=lineage.col.set)+xlim(umap.xlim)+ylim(umap.ylim)


plot.results$umap.EBD2 <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("EBD2")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("EBD2"))%>% mutate(od=factor(rename_EML,c("ELC","HLC","TLC","Undef"),ordered = T)) %>% arrange(desc(od)) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.5)+ggtitle("Stem-blastoids")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+geom_text(data.EM.umap.text.pos %>% filter(rename_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3.5,color="Black")+geom_text(data.LC.umap.text.pos %>% filter(pj=="EBD2") %>% filter(rename_EML %in% c("ELC","HLC","TLC")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")+xlim(umap.xlim)+ylim(umap.ylim)
plot.results.txt$umap.EBD2 <-  ggplot()+ggtitle("Stem-blastoids")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+geom_text(data.EM.umap.text.pos %>% filter(rename_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3.5,color="Black")+geom_text(data.LC.umap.text.pos %>% filter(pj=="EBD2") %>% filter(rename_EML %in% c("ELC","HLC","TLC")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")+xlim(umap.xlim)+ylim(umap.ylim)
plot.results.jpeg$umap.EBD2 <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("EBD2")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("EBD2"))%>% mutate(od=factor(rename_EML,c("ELC","HLC","TLC","Undef"),ordered = T)) %>% arrange(desc(od)) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.5)+ggtitle("")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+xlim(umap.xlim)+ylim(umap.ylim)

plot.results$umap.IBD2<-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("IBD2")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("IBD2"))  %>% mutate(od=factor(rename_EML,c("ELC","HLC","TLC","Undef"),ordered = T)) %>% arrange(desc(od)),mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.5)+ggtitle("Iblastoids")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+geom_text(data.EM.umap.text.pos %>% filter(rename_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3.5,color="Black")+geom_text(data.LC.umap.text.pos %>% filter(pj=="IBD2") %>% filter(rename_EML %in% c("ELC","HLC","TLC")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")+xlim(umap.xlim)+ylim(umap.ylim)
plot.results.txt$umap.IBD2 <-  ggplot()+ggtitle("Iblastoids")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+geom_text(data.EM.umap.text.pos %>% filter(rename_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3.5,color="Black")+geom_text(data.LC.umap.text.pos %>% filter(pj=="IBD2") %>% filter(rename_EML %in% c("ELC","HLC","TLC")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")+xlim(umap.xlim)+ylim(umap.ylim)
plot.results.jpeg$umap.IBD2  <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("IBD2")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("IBD2"))  %>% mutate(od=factor(rename_EML,c("ELC","HLC","TLC","Undef"),ordered = T)) %>% arrange(desc(od)),mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.5)+ggtitle("")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+xlim(umap.xlim)+ylim(umap.ylim)


plot.results$umap.nBGuoids <-  ggplot()+ geom_point(data.ob.umap %>% filter(! ((pj == "nBGuo") & (cellType!="EM"))) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter((pj%in% c("nBGuo")) & (cellType!="EM"))  %>% mutate(od=factor(rename_EML,c("ELC","HLC","TLC","Undef"),ordered = T)) %>% arrange(desc(od)),mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.5)+ggtitle("Guo-blastoids")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+geom_text(data.EM.umap.text.pos %>% filter(rename_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3.5,color="Black")+geom_text(data.LC.umap.text.pos %>% filter(pj=="nBGuo") %>% filter(rename_EML %in% c("ELC","HLC","TLC")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")+xlim(umap.xlim)+ylim(umap.ylim)
plot.results.txt$umap.nBGuoids  <-  ggplot()+ggtitle("Na-Blastoids")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+geom_text(data.EM.umap.text.pos %>% filter(rename_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3.5,color="Black")+geom_text(data.LC.umap.text.pos %>% filter(pj=="nBGuo") %>% filter(rename_EML %in% c("ELC","HLC","TLC")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")+xlim(umap.xlim)+ylim(umap.ylim)
plot.results.jpeg$umap.nBGuoids   <-  ggplot()+ geom_point(data.ob.umap %>% filter(! ((pj == "nBGuo") & (cellType!="EM"))) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter((pj%in% c("nBGuo")) & (cellType!="EM"))  %>% mutate(od=factor(rename_EML,c("ELC","HLC","TLC","Undef"),ordered = T)) %>% arrange(desc(od)),mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.5)+ggtitle("")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+xlim(umap.xlim)+ylim(umap.ylim)

plot.results$umap.nicolBla1 <-  ggplot()+ geom_point(data.ob.umap %>% filter( pj != "nicolBla" | subCT %in% c("okae_bts5","primed_H9","naive_H9") ) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("nicolBla") & rename_EML %in% c("ELC","HLC","TLC","Undef"))    %>% mutate(od=factor(rename_EML,c("ELC","HLC","TLC","Undef"),ordered = T)) %>% arrange(desc(od)),mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.5)+ggtitle("Nicolas-blastoids")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+geom_text(data.EM.umap.text.pos %>% filter(rename_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3.5,color="Black")+geom_text(data.LC.umap.text.pos %>% filter(pj=="nicolBla") %>% filter(rename_EML %in% c("ELC","HLC","TLC")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")+xlim(umap.xlim)+ylim(umap.ylim)
plot.results.txt$umap.nicolBla1 <-  ggplot()+ggtitle("Nicolas-blastoids")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+geom_text(data.EM.umap.text.pos %>% filter(rename_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3.5,color="Black")+geom_text(data.LC.umap.text.pos %>% filter(pj=="nicolBla") %>% filter(rename_EML %in% c("ELC","HLC","TLC")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")+xlim(umap.xlim)+ylim(umap.ylim)
plot.results.jpeg$umap.nicolBla1 <-  ggplot()+ geom_point(data.ob.umap %>% filter( pj != "nicolBla" | subCT %in% c("okae_bts5","primed_H9","naive_H9") ) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("nicolBla") & rename_EML %in% c("ELC","HLC","TLC","Undef"))    %>% mutate(od=factor(rename_EML,c("ELC","HLC","TLC","Undef"),ordered = T)) %>% arrange(desc(od)),mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.5)+ggtitle("Nicolas-blastoids")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+xlim(umap.xlim)+ylim(umap.ylim)


plot.results$umap.nicolBla2 <- ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("nicolBla") | devTime %in% c("24h","60h","90h")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("nicolBla") & devTime %in% c("naive_H9","okae_bts5","primed_H9")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.5)+ggtitle("Other cells from Nicolas dataset")+geom_text(data.EM.umap.text.pos %>% filter(rename_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3.5,color="Black")+  geom_text(data.LC.umap.text.pos %>% filter(pj=="nicolBla") %>% filter(!rename_EML %in% c("ELC","HLC","TLC","Undef")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")  +theme_classic()+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+xlim(umap.xlim)+ylim(umap.ylim)+scale_color_manual(values=lineage.col.set)

plot.results.txt$umap.nicolBla2 <- ggplot()+geom_text(data.EM.umap.text.pos %>% filter(rename_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3.5,color="Black")+  geom_text(data.LC.umap.text.pos %>% filter(pj=="nicolBla") %>% filter(!rename_EML %in% c("ELC","HLC","TLC","Undef")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")  +theme_classic()+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+xlim(umap.xlim)+ylim(umap.ylim)+scale_color_manual(values=lineage.col.set)

plot.results.jpeg$umap.nicolBla2 <- ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("nicolBla") | devTime %in% c("24h","60h","90h")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("nicolBla") & devTime %in% c("naive_H9","okae_bts5","primed_H9")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.5)+ggtitle("Other cells from Nicolas dataset") +theme_classic()+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+xlim(umap.xlim)+ylim(umap.ylim)+scale_color_manual(values=lineage.col.set)


#' integration results without PASE dataset
data.EM.umap.text.pos <- data.ob.NoJPF2019.umap  %>% filter(pj %in% c("D3post","CS7","SPH2016","nBGuo")) %>% filter(cellType=="EM") %>%group_by(rename_EML) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) %>% ungroup()  %>% filter(!rename_EML %in% c("3D_ICM","PSA-EPI","Unknown")) %>% mutate(UMAP_1=ifelse(rename_EML=="PriS",UMAP_1-1,UMAP_1)) 
data.LC.umap.text.pos <- data.ob.NoJPF2019.umap  %>% filter(pj %in% c("EBD2","IBD2","nBGuo","nicolBla")) %>% filter(cellType!="EM")  %>%group_by(rename_EML,pj) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) %>% ungroup()

# select the embryonic dataset
plot.results$umap.NOJPF.em <- ggplot()+ geom_point(data.ob.NoJPF2019.umap %>% filter(pj %in% c("EBD2","IBD2","JPF2019","nBGuo","nicolBla"))%>% filter(cellType!="EM") ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.5)+geom_point( data.ob.NoJPF2019.umap %>% filter(pj %in% c("D3post","CS7","SPH2016","nBGuo")) %>% filter(cellType=="EM") %>% mutate(od=factor(rename_EML,c("ICM-TE_trans","ICM","Prelineage","TE","Mesoderm","Epiblast","Amnion","Endoderm","PriS","ExE_Mes","3D_ICM","PSA-EPI"),ordered = T)) %>% arrange(desc(od)) ,mapping=aes(x=UMAP_1,y=UMAP_2,fill=rename_EML,shape=pj),size=2,color="grey")+ggtitle("Embryonic dataset")+geom_text(data.EM.umap.text.pos ,mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")+theme_classic()+theme(plot.title = element_text(hjust=0.5))+scale_fill_manual(values=lineage.col.set)+scale_shape_manual(values = pj.shape.set,labels = c('SPH2016'="Petropoulos et al., 2016",'D3post'="Xiang et al., 2019",'CS7'="Tyser et al., 2020",'nBGuo'="Yanagida et al., 2021"))+NoAxes()+NoLegend()
umap.xlim=ggplot_build(plot.results$umap.NOJPF.em)$layout$panel_scales_x[[1]]$range$range
umap.ylim=ggplot_build(plot.results$umap.NOJPF.em)$layout$panel_scales_y[[1]]$range$range
plot.results.txt$umap.NOJPF.em <- ggplot()+ggtitle("Embryonic dataset")+geom_text(data.EM.umap.text.pos ,mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")+theme_classic()+theme(plot.title = element_text(hjust=0.5))+scale_fill_manual(values=lineage.col.set)+scale_shape_manual(values = pj.shape.set,labels = c('SPH2016'="Petropoulos et al., 2016",'D3post'="Xiang et al., 2019",'CS7'="Tyser et al., 2020",'nBGuo'="Yanagida et al., 2021"))+NoAxes()+NoLegend()+xlim(umap.xlim)+ylim(umap.ylim)
plot.results.jpeg$umap.NOJPF.em  <- ggplot()+ geom_point(data.ob.NoJPF2019.umap %>% filter(pj %in% c("EBD2","IBD2","JPF2019","nBGuo","nicolBla"))%>% filter(cellType!="EM") ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.5)+geom_point( data.ob.NoJPF2019.umap %>% filter(pj %in% c("D3post","CS7","SPH2016","nBGuo")) %>% filter(cellType=="EM") %>% mutate(od=factor(rename_EML,c("ICM-TE_trans","ICM","Prelineage","TE","Mesoderm","Epiblast","Amnion","Endoderm","PriS","ExE_Mes","3D_ICM","PSA-EPI"),ordered = T)) %>% arrange(desc(od)) ,mapping=aes(x=UMAP_1,y=UMAP_2,fill=rename_EML,shape=pj),size=2,color="grey")+ggtitle("")+theme_classic()+theme(plot.title = element_text(hjust=0.5))+scale_fill_manual(values=lineage.col.set)+scale_shape_manual(values = pj.shape.set,labels = c('SPH2016'="Petropoulos et al., 2016",'D3post'="Xiang et al., 2019",'CS7'="Tyser et al., 2020",'nBGuo'="Yanagida et al., 2021"))+NoAxes()+NoLegend()+xlim(umap.xlim)+ylim(umap.ylim)


plot.results$umap.NOJPF.IBD2<-  ggplot()+ geom_point(data.ob.NoJPF2019.umap %>% filter(!pj %in% c("IBD2")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.NoJPF2019.umap %>% filter(pj%in% c("IBD2"))  %>% mutate(od=factor(rename_EML,c("ELC","HLC","TLC","Undef"),ordered = T)) %>% arrange(desc(od)),mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.5)+ggtitle("Iblastoids")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+geom_text(data.EM.umap.text.pos %>% filter(rename_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3.5,color="Black")+geom_text(data.LC.umap.text.pos %>% filter(pj=="IBD2") %>% filter(rename_EML %in% c("ELC","HLC","TLC")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")+xlim(umap.xlim)+ylim(umap.ylim)
plot.results.txt$umap.NOJPF.IBD2 <-  ggplot()+ggtitle("Iblastoids")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+geom_text(data.EM.umap.text.pos %>% filter(rename_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3.5,color="Black")+geom_text(data.LC.umap.text.pos %>% filter(pj=="IBD2") %>% filter(rename_EML %in% c("ELC","HLC","TLC")),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=4,color="Black",fontface="bold")+xlim(umap.xlim)+ylim(umap.ylim)
plot.results.jpeg$umap.NOJPF.IBD2  <-  ggplot()+ geom_point(data.ob.NoJPF2019.umap %>% filter(!pj %in% c("IBD2")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.NoJPF2019.umap %>% filter(pj%in% c("IBD2"))  %>% mutate(od=factor(rename_EML,c("ELC","HLC","TLC","Undef"),ordered = T)) %>% arrange(desc(od)),mapping=aes(x=UMAP_1,y=UMAP_2,col=rename_EML),size=0.5)+ggtitle("")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+xlim(umap.xlim)+ylim(umap.ylim)



#' Cross-species integration results
for ( n in names(human.NHP.IT.umap.list)) {
  data.temp.UMAP <- human.NHP.IT.umap.list[[n]]
  data.temp.UMAP <- data.temp.UMAP  %>% mutate(spec=ifelse(pj !="NHP","human","NHP")) %>% mutate(rename_EML=ifelse( rename_EML %in% c("Amnion"), paste0(rename_EML,"(",spec,")"),rename_EML))
  label1=unlist(strsplit(n,"\\."))[[1]]
  label2=unlist(strsplit(n,"\\."))[[2]]
  me=unlist(strsplit(n,"\\."))[[3]]
  temp.plot <- FunCSUMAP(data.temp.UMAP,label1,label2)
  plot.results[[n]] <- temp.plot[[1]]
  plot.results.txt[[n]] <- temp.plot[[3]]
  plot.results.jpeg[[n]] <- temp.plot[[2]]
  plot.results[[paste0(n,".spec")]] <- temp.plot[[4]]
}
plot.results$CS.legend <- data.frame(X=(1:(cross.lineage.col.set %>% length())),Y=1:(cross.lineage.col.set %>% length()),anno=names(cross.lineage.col.set)) %>% mutate(spec=ifelse(anno %in% c("Amnion(NHP)","E-Amnion"),"NHP","human")) %>% ggplot+geom_point(mapping=aes(x=X,y=Y,col=anno,shape=spec))+scale_color_manual(values=cross.lineage.col.set)+scale_shape_manual(values = species.shape.set)+theme_classic()
plot.results$CS.legend <- plot.results$CS.legend %>% ggpubr::get_legend() %>% ggpubr::as_ggplot()

#' AmVSTE marker gene expression in Amnion and all TE
temp.M <- data.ob.umap %>% filter(pj !="D3post") %>% select(cell,rename_EML,pj) %>% filter(rename_EML %in% c("Amnion","TE","TLC")) %>% mutate(anno=paste(pj,rename_EML,sep="_")) %>% select(cell,anno) %>% bind_rows(data.ob.umap %>% filter(pj =="D3post") %>% filter(rename_EML %in% "TE") %>% mutate(anno=ifelse(devTime %in% c("E6","E7"),"D3_PreTE","D3_PostTE")) %>% select(cell,anno)) %>% bind_rows(meta.filter %>% filter(EML=="TE" & pj=="Blakeley") %>% select(cell) %>% mutate(anno="Blakeley_TE"))

#ds
temp.M1 <- temp.M %>% filter(anno %in% c("CS7_Amnion","Blakeley_TE","nBGuo_TE","SPH2016_TE","D3_PreTE","D3_PostTE")) %>% split(.,.$anno) %>% lapply(function(x) FunMaSF(x,50)) %>% do.call("bind_rows",.)
temp.M2 <- temp.M %>% filter(anno %in% c("EBD2_TLC","IBD2_TLC","nBGuo_TLC","nicolBla_TLC")) %>% split(.,.$anno) %>% lapply(function(x) FunMaSF(x,50)) %>% do.call("bind_rows",.)

#' em heatmap
#+ fig.show='hide'
temp.anno <- temp.M1 %>% mutate(od=factor(anno,c("CS7_Amnion","Blakeley_TE","SPH2016_TE","nBGuo_TE","D3_PreTE","D3_PostTE"),ordered = T)) %>% arrange(od) %>% tibble::column_to_rownames("cell") 
temp.sel.cell <- rownames(temp.anno)
temp.sel.gene <- hm.AMvsTE.out %>% arrange(Type) %>% pull(gene)
row.gaps <-   hm.AMvsTE.out %>% arrange(Type)  %>% filter(Type=="Amnion") %>% nrow()
temp.exp <- lognormExp.mBN[temp.sel.gene,temp.sel.cell] 
temp.sel.exp <- t(apply(temp.exp,1,scale))
colnames(temp.sel.exp) <- colnames(temp.exp)
rownames(temp.sel.exp) <- rownames(temp.exp)
temp.sel.exp[temp.sel.exp>zs.limit] <- zs.limit
temp.sel.exp[temp.sel.exp<  (-1*zs.limit)] <- -1*zs.limit
gaps <- 10
plot.results$MK.heatmap.EM <- pheatmap(temp.sel.exp,cluster_rows=F,,cluster_cols=F,scale="none",annotation_col=temp.anno,show_colnames=F,show_rownames=F,color=heat.col, fontsize_row=7,border_color="NA",main="Embryo",gaps_row = rep(row.gaps,each=2))%>%ggplotify::as.ggplot()

#' TLC heatmap
#+ fig.show='hide'
temp.anno <- temp.M2 %>% mutate(od=factor(anno,c("IBD2_TLC","EBD2_TLC","nicolBla_TLC","nBGuo_TLC"),ordered = T)) %>% arrange(od) %>% tibble::column_to_rownames("cell") 
temp.sel.cell <- rownames(temp.anno)
temp.sel.gene <- hm.AMvsTE.out %>% arrange(Type) %>% pull(gene)
row.gaps <-   hm.AMvsTE.out %>% arrange(Type)  %>% filter(Type=="Amnion") %>% nrow()
temp.exp <- lognormExp.mBN[temp.sel.gene,temp.sel.cell] 
temp.sel.exp <- t(apply(temp.exp,1,scale))
colnames(temp.sel.exp) <- colnames(temp.exp)
rownames(temp.sel.exp) <- rownames(temp.exp)
temp.sel.exp[temp.sel.exp>zs.limit] <- zs.limit
temp.sel.exp[temp.sel.exp<  (-1*zs.limit)] <- -1*zs.limit
gaps <- 10
plot.results$MK.heatmap.blastoids <- pheatmap(temp.sel.exp,cluster_rows=F,,cluster_cols=F,scale="none",annotation_col=temp.anno,show_colnames=F,show_rownames=F,color=heat.col, fontsize_row=7,border_color="NA",main="Embryo",gaps_row = rep(row.gaps,each=2))%>%ggplotify::as.ggplot()

# combined the TLC cells together
temp.M3 <- temp.M %>% filter(anno %in% c("EBD2_TLC","IBD2_TLC","nBGuo_TLC","nicolBla_TLC"))
temp.sel.gene <- hm.AMvsTE.out %>% arrange(Type) %>% pull(gene)

row.gaps <-   hm.AMvsTE.out %>% arrange(Type)  %>% filter(Type=="Amnion") %>% nrow()
temp.exp <- temp.M3 %>% split(.,.$anno) %>% lapply(function(x) {lognormExp.mBN[temp.sel.gene,x$cell]%>% rowMeans() %>% as.data.frame()%>%  tibble::rownames_to_column("cell") %>% setNames(c("gene",unique(x$anno))) %>% tbl_df() %>% tibble::column_to_rownames("gene")%>% return()})  %>% do.call("bind_cols",.)
temp.exp <- temp.exp[,c("IBD2_TLC","EBD2_TLC","nicolBla_TLC","nBGuo_TLC")]
temp.sel.exp <- t(apply(temp.exp,1,scale))
colnames(temp.sel.exp) <- colnames(temp.exp)
rownames(temp.sel.exp) <- rownames(temp.exp)
temp.sel.exp[temp.sel.exp>zs.limit] <- zs.limit
temp.sel.exp[temp.sel.exp<  (-1*zs.limit)] <- -1*zs.limit

temp.sup.exp <- temp.sel.exp[1,,drop=F]
rownames(temp.sup.exp)="MakeColorKey"
temp.sup.exp <- c(-1*zs.limit,0,0,zs.limit)

gaps <- 10
plot.results$MK.heatmap.blastoids.meanlogExp <- pheatmap(temp.sup.exp%>% rbind(temp.sel.exp ),cluster_rows=F,,cluster_cols=F,scale="none",show_colnames=T,show_rownames=F,color=heat.col, fontsize_row=7,border_color="NA",main="Embryo",gaps_row = c(1,1,1,1,rep(row.gaps+1,each=2)))%>%ggplotify::as.ggplot()



#' AmVSTE marker gene expression in Amnion and all TE but in NHP
#+ fig.show='hide'
temp.anno <-NHP.meta %>% mutate(anno=paste(devTime,EML,sep="_")) %>% filter(anno %in% c("E14_Amnion2","E10_TE","E12_TE","E14_TE")) %>% select(cell,anno) 
#ds
temp.anno <- temp.anno  %>% split(.,.$anno) %>% lapply(function(x) FunMaSF(x,100)) %>% do.call("bind_rows",.)
temp.anno <- temp.anno %>% mutate(od=factor(anno,c("E14_Amnion2","E10_TE","E12_TE","E14_TE"),ordered=T)) %>% arrange(od)  %>% tibble::column_to_rownames("cell") 
temp.sel.cell <- rownames(temp.anno)
temp.sel.gene <- hm.AMvsTE.out %>% arrange(Type) %>% pull(gene) %>% intersect(GI.list$ov)
row.gaps <-   hm.AMvsTE.out %>% arrange(Type)  %>% filter(Type=="Amnion")%>% pull(gene) %>% intersect(GI.list$ov) %>% length()
temp.exp <- NHP.lognormExp.mBN[temp.sel.gene,temp.sel.cell] 
temp.sel.exp <- t(apply(temp.exp,1,scale))
colnames(temp.sel.exp) <- colnames(temp.exp)
rownames(temp.sel.exp) <- rownames(temp.exp)
temp.sel.exp[temp.sel.exp>zs.limit] <- zs.limit
temp.sel.exp[temp.sel.exp<  (-1*zs.limit)] <- -1*zs.limit

plot.results$NHP.MK.heatmap.EM <- pheatmap(temp.sel.exp,cluster_rows=F,,cluster_cols=F,scale="none",annotation_col=temp.anno,show_colnames=F,show_rownames=F,color=heat.col, fontsize_row=7,border_color="NA",main="Embryo(NHP)",gaps_row = rep(row.gaps,each=2))%>%ggplotify::as.ggplot()



#' VlnPlot for key marker genes
temp.sel.gene <- c("ISL1","GABRP","IGFBP5","GATA2","GCM1","BIN2")
temp.M <- data.ob.umap %>% filter(EML!="Tsw-AMLC")%>% filter(pj !="D3post") %>% select(cell,rename_EML,pj) %>% filter(rename_EML %in% c("Amnion","TE","TLC","AMLC")) %>% mutate(anno=paste(pj,rename_EML,sep="_")) %>% select(cell,anno) %>% bind_rows(data.ob.umap %>% filter(pj =="D3post") %>% filter(rename_EML %in% "TE") %>% mutate(anno=ifelse(devTime %in% c("E6","E7"),"D3_PreTE","D3_PostTE")) %>% select(cell,anno)) %>% bind_rows(meta.filter %>% filter(EML=="TE" & pj=="Blakeley") %>% select(cell) %>% mutate(anno="Blakeley_TE"))


temp.sel.exp <- lognormExp.mBN[temp.sel.gene ,temp.M$cell]%>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,logExp,-Gene) %>% inner_join(temp.M,by="cell") %>% mutate(od=factor(anno,c("CS7_Amnion","Blakeley_TE","SPH2016_TE","nBGuo_TE","D3_PreTE","D3_PostTE","JPF2019_AMLC","IBD2_TLC","EBD2_TLC","nicolBla_TLC","nBGuo_TLC"),ordered = T)) %>% arrange(od)  

#"CS7_Amnion","Blakeley_TE","SPH2016_TE","nBGuo_TE","D3_PreTE","D3_PostTE"
for (g in temp.sel.gene ) {
  plot.results[[paste0("MK.Vln.",g)]] <- temp.sel.exp%>% filter(Gene==g) %>% ggplot(mapping=aes(x=od,y=logExp)) + geom_violin(mapping=aes(fill=anno),color="black",scale = "width")+geom_jitter(size=0.05)+ theme_classic() + theme(axis.text.x=element_text(angle = 40,hjust = 1)) +NoLegend()+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="bold"))+ylab("Expression level")+xlab("")+ylab("")
  plot.results.jpeg[[paste0("MK.Vln.",g)]] <- temp.sel.exp%>% filter(Gene==g) %>% ggplot(mapping=aes(x=od,y=logExp)) + geom_violin(mapping=aes(fill=anno),color="black",scale = "width")+ theme_classic() + theme(axis.text.x=element_text(angle = 40,hjust = 1)) +NoLegend()+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="bold"))+ylab("Expression level")+xlab("")+ylab("")
}

#' VlnPlot for key marker genes in NHP
temp.sel.gene <- c("ISL1","GABRP","IGFBP5","GATA2","GCM1","BIN2")
temp.M <- NHP.meta %>% mutate(anno=paste(devTime,EML,sep="_")) %>% filter(anno %in% c("E14_Amnion2","E10_TE","E12_TE","E14_TE")) %>% select(cell,anno) 

temp.sel.exp <- NHP.lognormExp.mBN[temp.sel.gene ,temp.M$cell]%>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,logExp,-Gene) %>% inner_join(temp.M,by="cell") %>% mutate(od=factor(anno,c("E14_Amnion2","E10_TE","E12_TE","E14_TE"),ordered = T)) %>% arrange(od)  

for (g in temp.sel.gene ) {
  plot.results[[paste0("NHP.MK.Vln.",g)]] <- temp.sel.exp%>% filter(Gene==g) %>% ggplot(mapping=aes(x=od,y=logExp)) + geom_violin(mapping=aes(fill=anno),color="black",scale = "width")+geom_jitter(size=0.05)+ theme_classic() + theme(axis.text.x=element_text(angle = 40,hjust = 1)) +NoLegend()+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="bold"))+ylab("Expression level")+xlab("")+ylab("")
  plot.results.jpeg[[paste0("NHP.MK.Vln.",g)]] <- temp.sel.exp%>% filter(Gene==g) %>% ggplot(mapping=aes(x=od,y=logExp)) + geom_violin(mapping=aes(fill=anno),color="black",scale = "width")+ theme_classic() + theme(axis.text.x=element_text(angle = 40,hjust = 1)) +NoLegend()+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="bold"))+ylab("Expression level")+xlab("")+ylab("")
}

#' PCA results for pseudo-bulk sequencing dataset
pca.it <- pca.out$human
pc1.imp <- round(100*summary(pca.it)$importance[2,1],digits=2) 
pc2.imp <-round(100*summary(pca.it)$importance[2,2],digits=2) 
temp <- pca.it$x[,c("PC1","PC2")] %>% as.data.frame()%>% tibble::rownames_to_column("Sample") %>% mutate(Sample=ifelse(Sample=="D3post_TE_pre","D3_PreTE",Sample))%>% mutate(Sample=ifelse(Sample=="D3post_TE_post","D3_PostTE",Sample)) %>% separate(Sample,c("pj","anno"),sep="_",remove=F) %>% mutate(anno=ifelse(anno=="Tsw","AMLC",anno))%>% mutate(anno=ifelse(anno %in% c("PreTE","PostTE"),"TE",anno)) 
plot.results$psd.human.pca <- ggplot(temp,mapping=aes(x=PC1,y=PC2))+geom_point(mapping=aes(col=pj,shape=anno),size=3)+ggrepel::geom_text_repel(mapping=aes(label=Sample))+xlab(paste("PC1(",pc1.imp,"% Proportion of Variance)"))+ylab(paste("PC2(",pc2.imp,"% Proportion of Variance)"))+theme_classic()+ylim(-15,25)+xlim(-25,50)+scale_color_manual(values=pj.col)

pca.it <- pca.out$HuMonk
pc1.imp <- round(100*summary(pca.it)$importance[2,1],digits=2) 
pc2.imp <-round(100*summary(pca.it)$importance[2,2],digits=2) 
temp <- pca.it$x[,c("PC1","PC2")] %>% as.data.frame()%>% tibble::rownames_to_column("Sample") %>% mutate(Sample=ifelse(Sample=="D3post_TE_pre","D3_PreTE",Sample))%>% mutate(Sample=ifelse(Sample=="D3post_TE_post","D3_PostTE",Sample)) %>% mutate(Sample=gsub("NHP_D","NHPD",Sample))%>% separate(Sample,c("pj","anno"),sep="_",remove=F) %>% mutate(anno=ifelse(anno=="Tsw","AMLC",anno)) %>% mutate(anno=ifelse(anno %in% c("PreTE","PostTE"),"TE",anno))  %>% mutate(pj=ifelse(pj %in% c("NHPD10","NHPD12","NHPD14"),"NHP",pj))
plot.results$psd.NHP.human.pca <- ggplot(temp,mapping=aes(x=PC1,y=PC2))+geom_point(mapping=aes(col=pj,shape=anno))+ggrepel::geom_text_repel(mapping=aes(label=Sample))+xlab(paste("PC1(",pc1.imp,"% Proportion of Variance)"))+ylab(paste("PC2(",pc2.imp,"% Proportion of Variance)"))+theme_classic()+ylim(-25,50)+xlim(-25,50)+scale_color_manual(values=pj.col)



#temp.label.gene <- c("ISL1","GABRP","IGFBP5","RARRES2","LRRN1","GATA2","GCM1","BIN2","CYP11A1","CYP19A1")


plot.results$new.anno.EBD2 <- data.ob.umap %>% filter(pj %in% c("EBD2"))  %>% group_by(cluster_EML) %>% summarise(nCell=n_distinct(cell))  %>% mutate(Perc=nCell/sum(nCell)) %>% mutate(od=factor(cluster_EML,c("ELC","HLC","TLC","MeLC","AMLC","Undef"),ordered = T)) %>% arrange(od) %>% mutate(od=paste0(od,"\n(",nCell,")")) %>% mutate(od=factor(od,od,,ordered = T))  %>% ggplot()+geom_bar(mapping=aes(x=od,y=Perc*100,fill=cluster_EML),width=0.5,stat = "identity")+theme_classic()+ylab("% Proportion of cells")+xlab("")+ ggtitle("Stem-blastoids")+theme(plot.title = element_text(hjust=0.5))+ theme(legend.position="top",legend.title = element_blank())+scale_fill_manual(values=lineage.col.set)+ylim(0,67.5)+NoLegend()

plot.results$new.anno.IBD2 <- data.ob.umap %>% filter(pj %in% c("IBD2"))  %>% group_by(cluster_EML) %>% summarise(nCell=n_distinct(cell))  %>% mutate(Perc=nCell/sum(nCell)) %>% mutate(od=factor(cluster_EML,c("ELC","HLC","TLC","MeLC","AMLC","Undef"),ordered = T)) %>% arrange(od) %>% mutate(od=paste0(od,"\n(",nCell,")")) %>% mutate(od=factor(od,od,,ordered = T))  %>% ggplot()+geom_bar(mapping=aes(x=od,y=Perc*100,fill=cluster_EML),width=0.5,stat = "identity")+theme_classic()+ylab("% Proportion of cells")+xlab("")+ ggtitle("Iblastoids")+theme(plot.title = element_text(hjust=0.5))+ theme(legend.position="right",legend.title = element_blank())+scale_fill_manual(values=lineage.col.set)+ylim(0,67.5)

plot.results$new.anno.nBGuo <- data.ob.umap %>% filter(pj %in% c("nBGuo")) %>% filter(cellType!="EM")  %>% group_by(cluster_EML) %>% summarise(nCell=n_distinct(cell)) %>% full_join(data.frame(cluster_EML=c("ELC","HLC","TLC","MeLC","AMLC","Undef"),nBase=0) %>% tbl_df() %>% mutate_all(as.vector),by="cluster_EML") %>% mutate(nCell=ifelse(is.na(nCell),nBase,nCell)) %>% select(-nBase)  %>% mutate(Perc=nCell/sum(nCell))  %>% mutate(od=factor(cluster_EML,c("ELC","HLC","TLC","MeLC","AMLC","Undef"),ordered = T)) %>% arrange(od) %>% mutate(od=paste0(od,"\n(",nCell,")")) %>% mutate(od=factor(od,od,,ordered = T))  %>% ggplot()+geom_bar(mapping=aes(x=od,y=Perc*100,fill=cluster_EML),width=0.5,stat = "identity")+theme_classic()+ylab("% Proportion of cells")+xlab("")+ ggtitle("Nai-blastoids")+theme(plot.title = element_text(hjust=0.5))+ theme(legend.position="right",legend.title = element_blank())+scale_fill_manual(values=lineage.col.set)+ylim(0,67.5)


plot.results$new.anno.nicolBla <- data.ob.umap %>% filter(pj %in% c("nicolBla")) %>% filter(cellType=="blastoids")  %>% group_by(cluster_EML) %>% summarise(nCell=n_distinct(cell)) %>% full_join(data.frame(cluster_EML=c("ELC","HLC","TLC","MeLC","AMLC","Undef"),nBase=0) %>% tbl_df() %>% mutate_all(as.vector),by="cluster_EML") %>% mutate(nCell=ifelse(is.na(nCell),nBase,nCell)) %>% select(-nBase)  %>% mutate(Perc=nCell/sum(nCell))  %>% mutate(od=factor(cluster_EML,c("ELC","HLC","TLC","MeLC","AMLC","Undef"),ordered = T)) %>% arrange(od) %>% mutate(od=paste0(od,"\n(",nCell,")")) %>% mutate(od=factor(od,od,,ordered = T))  %>% ggplot()+geom_bar(mapping=aes(x=od,y=Perc*100,fill=cluster_EML),width=0.5,stat = "identity")+theme_classic()+ylab("% Proportion of cells")+xlab("")+ ggtitle("Nicolas-blastoids")+theme(plot.title = element_text(hjust=0.5))+ theme(legend.position="right",legend.title = element_blank())+scale_fill_manual(values=lineage.col.set)+ylim(0,67.5)

#' the cluster information
data.cluster.text.pos <- data.ob.umap %>%group_by(seurat_clusters) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) %>% ungroup()
plot.results$umap.cluster <- ggplot()+ geom_point(data.ob.umap ,mapping=aes(x=UMAP_1,y=UMAP_2,color=seurat_clusters),size=0.5)+geom_text(data.cluster.text.pos  ,mapping=aes(x=UMAP_1,y=UMAP_2,label=seurat_clusters),size=4,color="Black")+theme_classic()+NoAxes()+NoLegend()#+scale_color_manual(values=cluster.col.set)
plot.results.jpeg$umap.cluster <- ggplot()+ geom_point(data.ob.umap ,mapping=aes(x=UMAP_1,y=UMAP_2,color=seurat_clusters),size=0.5)+theme_classic()+NoAxes()+NoLegend()

#' dotplot for lineage markers based on new annotation

temp.M <-  data.ob.umap %>% filter(cluster_EML %in% c("Epiblast","ELC","TE","TLC","Endoderm","HLC","Mesoderm","MeLC","Amnion","AMLC"))
temp.M.EM <- temp.M %>% filter(seqType=="smt2") %>% filter(cluster_EML %in% c("Epiblast","TE","Endoderm","Mesoderm","Amnion"))
temp.mk.sel <- lineage.mk$sig %>% group_by(set)%>% filter(power >0.6) %>% top_n(10,power) %>%  dplyr::rename(Gene=gene) %>% mutate(od=factor(set,c("Epiblast","Endoderm","TE","Mesoderm","Amnion"),ordered = T)) %>% arrange(od,desc(power)) %>% select(-od)


temp <- lognormExp.mBN[temp.mk.sel $Gene ,temp.M.EM$cell]  %>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,logExp,-Gene) %>% inner_join(temp.M.EM %>% select(cell,cluster_EML) %>% dplyr::rename(cell.cluster=cluster_EML),by="cell") %>% inner_join(temp.mk.sel %>% select(Gene,set) %>% dplyr::rename(Gene.cluster=set) ,by="Gene") 


temp.input <- temp %>% group_by(cell.cluster,Gene.cluster,Gene) %>% summarise(meanlogExp=mean(logExp),nCell=n_distinct(cell)) %>% left_join(temp %>% group_by(cell.cluster,Gene.cluster,Gene) %>% filter(logExp >0)%>% summarise(nExpCell=n_distinct(cell)),by=c("cell.cluster","Gene.cluster","Gene"))  %>% mutate(nExpCell=ifelse(is.na(nExpCell),0,nExpCell))%>% mutate(nExp.ct=nExpCell/nCell) %>% group_by(Gene) %>% mutate(sv=(meanlogExp-mean(meanlogExp))/sd(meanlogExp)) %>% mutate(nExp.ct= as.vector(cut(nExp.ct*100,c(0,5,25,50,75,100),label=c(0,25,50,75,100)))) %>% replace(.,is.na(.),"0") %>% mutate(nExp.ct=as.numeric(nExp.ct))


plot.results$lm.em.exp.dot <- temp.input %>% mutate(Gene=factor(Gene,rev(temp.mk.sel $Gene),ordered = T)) %>% mutate(cell.cluster=factor(cell.cluster,c("Epiblast","Endoderm","TE","Mesoderm","Amnion"),ordered = T))%>% mutate(Gene.cluster=factor(Gene.cluster,c("Epiblast","Endoderm","TE","Mesoderm","Amnion"),ordered = T)) %>%  ggplot(aes(x=cell.cluster, y = Gene, color =sv, size = nExp.ct)) + geom_point() +facet_grid(Gene.cluster~.,,scale="free_y") + scale_color_viridis_c(name = '') + cowplot::theme_cowplot() + theme(axis.line  = element_blank()) +theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +ylab('') +theme(axis.ticks = element_blank()) +scale_size(name="",breaks = c(0,25,50,75,100),range = c(0.1, 1.5))+xlab("")+theme(axis.text.y = element_text(size = 5))+theme(axis.text.x = element_text(size = 7))+ theme(legend.position="top")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(strip.background=element_blank(),strip.text.y=element_blank())+NoLegend()+scale_y_discrete(position = "right")

plot.results$lm.em.exp.dot.legend <-temp.input  %>%  ggplot(aes(x=cell.cluster, y = Gene, color =sv, size = nExp.ct)) + geom_point()  + scale_color_viridis_c(name = '') + cowplot::theme_cowplot()  +scale_size(name="",breaks = c(0,25,50,75,100),range = c(0.1, 1.5))+ theme(legend.position="top")+theme(legend.text = element_text( size = 6),legend.box="vertical", legend.margin=margin())  #Scaled Expression  % cells expressing gene
plot.results$lm.em.exp.dot.legend<- plot.results$lm.em.exp.dot.legend%>% ggpubr::get_legend() %>% ggpubr::as_ggplot()


#' for blastoids
temp.M.blastoids <- temp.M  %>% filter(cluster_EML %in% c("ELC","TLC","HLC","MeLC","AMLC")) %>% filter(pj %in% c("EBD2","IBD2","nBGuo","nicolBla"))

temp.mk.sel <- lineage.mk$sig %>% group_by(set)%>% filter(power >0.6) %>% top_n(10,power) %>%  dplyr::rename(Gene=gene) %>% mutate(od=factor(set,c("Epiblast","Endoderm","TE","Mesoderm","Amnion"),ordered = T)) %>% arrange(od,desc(power)) %>% select(-od)


temp <- lognormExp.mBN[temp.mk.sel$Gene ,temp.M.blastoids$cell]  %>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,logExp,-Gene) %>% inner_join(temp.M %>% select(cell,cluster_EML,pj) %>% dplyr::rename(cell.cluster=cluster_EML),by="cell") %>% inner_join(temp.mk.sel%>% select(Gene,set) %>% dplyr::rename(Gene.cluster=set) ,by="Gene") 


temp.input <- temp %>% group_by(cell.cluster,Gene.cluster,pj,Gene) %>% summarise(meanlogExp=mean(logExp),nCell=n_distinct(cell)) %>% left_join(temp %>% group_by(cell.cluster,Gene.cluster,Gene,pj) %>% filter(logExp >0)%>% summarise(nExpCell=n_distinct(cell)),by=c("cell.cluster","Gene.cluster","Gene","pj"))  %>% mutate(nExpCell=ifelse(is.na(nExpCell),0,nExpCell))%>% mutate(nExp.ct=nExpCell/nCell) %>% group_by(Gene,pj) %>% mutate(sv=(meanlogExp-mean(meanlogExp))/sd(meanlogExp)) %>% mutate(nExp.ct= as.vector(cut(nExp.ct*100,c(0,5,25,50,75,100),label=c(0,25,50,75,100)))) %>% replace(.,is.na(.),"0") %>% mutate(nExp.ct=as.numeric(nExp.ct))


plot.results$lm.EBD2.exp.dot <- temp.input %>%filter(pj=="EBD2")  %>% mutate(Gene=factor(Gene,rev(temp.mk.sel$Gene),ordered = T)) %>% mutate(cell.cluster=factor(cell.cluster,c("ELC","HLC","TLC","MeLC","AMLC","Undef"),ordered = T))%>% mutate(Gene.cluster=factor(Gene.cluster,c("Epiblast","Endoderm","TE","Mesoderm","Amnion"),ordered = T)) %>%  ggplot(aes(x=cell.cluster, y = Gene, color =sv, size = nExp.ct)) + geom_point() +facet_grid(Gene.cluster~.,,scale="free_y") + scale_color_viridis_c(name = 'Scaled Expression ') + cowplot::theme_cowplot() + theme(axis.line  = element_blank()) +theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +ylab('') +theme(axis.ticks = element_blank()) +scale_size(name="% cells expressing",breaks = c(0,25,50,75,100),range = c(0.1,3))+xlab("")+theme(axis.text.y = element_text(size = 5))+theme(axis.text.x = element_text(size = 7))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(strip.background=element_blank(),strip.text.y=element_blank())+NoLegend()+ggtitle("Stem-Blastoids")+theme(plot.title = element_text(hjust=0.5,face="bold"))


plot.results$lm.IBD2.exp.dot <- temp.input %>%filter(pj=="IBD2")  %>% mutate(Gene=factor(Gene,rev(temp.mk.sel$Gene),ordered = T)) %>% mutate(cell.cluster=factor(cell.cluster,c("ELC","HLC","TLC","MeLC","AMLC","Undef"),ordered = T))%>% mutate(Gene.cluster=factor(Gene.cluster,c("Epiblast","Endoderm","TE","Mesoderm","Amnion"),ordered = T)) %>%  ggplot(aes(x=cell.cluster, y = Gene, color =sv, size = nExp.ct)) + geom_point() +facet_grid(Gene.cluster~.,,scale="free_y") + scale_color_viridis_c(name = 'Scaled Expression ') + cowplot::theme_cowplot() + theme(axis.line  = element_blank()) +theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +ylab('') +theme(axis.ticks = element_blank()) +scale_size(name="% cells expressing",breaks = c(0,25,50,75,100),range = c(0.1,3))+xlab("")+theme(axis.text.y = element_text(size = 5))+theme(axis.text.x = element_text(size = 7))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(strip.background=element_blank(),strip.text.y=element_blank())+NoLegend()+ggtitle("i-Blastoids")+theme(plot.title = element_text(hjust=0.5,face="bold"))


plot.results$lm.nBGuo.exp.dot <- temp.input %>%filter(pj=="nBGuo")  %>% mutate(Gene=factor(Gene,rev(temp.mk.sel$Gene),ordered = T)) %>% mutate(cell.cluster=factor(cell.cluster,c("ELC","HLC","TLC","MeLC","AMLC","Undef"),ordered = T))%>% mutate(Gene.cluster=factor(Gene.cluster,c("Epiblast","Endoderm","TE","Mesoderm","Amnion"),ordered = T)) %>%  ggplot(aes(x=cell.cluster, y = Gene, color =sv, size = nExp.ct)) + geom_point() +facet_grid(Gene.cluster~.,,scale="free_y") + scale_color_viridis_c(name = 'Scaled Expression ') + cowplot::theme_cowplot() + theme(axis.line  = element_blank()) +theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +ylab('') +theme(axis.ticks = element_blank()) +scale_size(name="% cells expressing",breaks = c(0,25,50,75,100),range = c(0.1,3))+xlab("")+theme(axis.text.y = element_text(size = 5))+theme(axis.text.x = element_text(size = 7))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(strip.background=element_blank(),strip.text.y=element_blank())+NoLegend()+ggtitle("Guo-Blastoids")+theme(plot.title = element_text(hjust=0.5,face="bold"))

plot.results$lm.nicolBla.exp.dot <- temp.input %>%filter(pj=="nicolBla")  %>% mutate(Gene=factor(Gene,rev(temp.mk.sel$Gene),ordered = T)) %>% mutate(cell.cluster=factor(cell.cluster,c("ELC","HLC","TLC","MeLC","AMLC","Undef"),ordered = T))%>% mutate(Gene.cluster=factor(Gene.cluster,c("Epiblast","Endoderm","TE","Mesoderm","Amnion"),ordered = T)) %>%  ggplot(aes(x=cell.cluster, y = Gene, color =sv, size = nExp.ct)) + geom_point() +facet_grid(Gene.cluster~.,,scale="free_y") + scale_color_viridis_c(name = 'Scaled Expression ') + cowplot::theme_cowplot() + theme(axis.line  = element_blank()) +theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +ylab('') +theme(axis.ticks = element_blank()) +scale_size(name="% cells expressing",breaks = c(0,25,50,75,100),range = c(0.1,3))+xlab("")+theme(axis.text.y = element_text(size = 5))+theme(axis.text.x = element_text(size = 7))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(strip.background=element_blank(),strip.text.y=element_blank())+NoLegend()+ggtitle("Nicolas-Blastoids")+theme(plot.title = element_text(hjust=0.5,face="bold"))

#' confirm the new annotation 
#' for "ICM" from D3post
temp.M <- meta.filter %>% filter(pj=="D3post") %>% inner_join(ref.meta %>% mutate(cell=ref_cell,new_EML=ref_EML) %>% select(cell,new_EML),by="cell") %>% mutate(anno=ifelse(EML=="ICM","TE(Prev-ICM)",new_EML))%>% mutate(anno=ifelse(EML=="EPI" & new_EML=="Mes","Mes(Prev-Epi)",anno)) %>% mutate(anno=ifelse(EML=="PSA-EPI" & new_EML=="Mes","Mes(Prev-PSA-EPI)",anno)) %>% mutate(anno=ifelse(EML=="PSA-EPI" & new_EML=="TE","TE(Prev-PSA-EPI)",anno))%>% mutate(anno=ifelse(EML=="PSA-EPI" & new_EML=="Epi","Epi(Prev-PSA-EPI)",anno))%>% mutate(anno=ifelse(EML=="PSA-EPI" & new_EML=="Endoderm","Endoderm(Prev-PSA-EPI)",anno))
data.temp <- CreateSeuratObject(counts.filter[,temp.M %>% pull(cell)], meta.data =(temp.M  %>% tibble::column_to_rownames("cell")))  %>% NormalizeData(verbose = FALSE)
temp.sel.gene <- c("NANOG","BMP2","GATA2","PMP22")#"POU5F1","NANOG","DPPA3","BMP2","GATA2","GATA3","PMP22","VIM"
temp.sel.exp <- data.temp@assays$RNA@data[temp.sel.gene ,temp.M$cell] %>% as.data.frame() %>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,logExp,-Gene) %>% inner_join(temp.M %>% select(cell,seqType,anno),by="cell")
temp.plot <- list()
for (g in temp.sel.gene ) {
  temp.plot[[g]] <- temp.sel.exp%>% filter(Gene==g) %>% ggplot(mapping=aes(x=anno,y=logExp)) + geom_violin(mapping=aes(fill=anno),color="black",scale = "width")+geom_jitter(size=0.05)+ theme_classic() + theme(axis.text.x=element_text(angle = 40,hjust = 1)) +NoLegend()+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="bold"))+ylab("Expression level")+xlab("")+ylab("") +scale_fill_manual(values=recheck.lineage.col.set)
  plot.results[[paste0("ReCheck.D3post.MK.Vln.",g)]] <- temp.plot[[g]] 
}


#' for Intermediate cells from CS7

plot.results[["CS7.recheck.old.anno"]] <- data.AEP.UMAP %>% ggplot()+geom_point(mapping=aes(x=UMAP_1,y=UMAP_2,col=EML))+theme_classic()+NoAxes()+NoLegend()+ggtitle("Previous annotation")+theme(plot.title = element_text(hjust=0.5))+geom_text(data=data.AEP.UMAP %>% group_by(EML) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) %>% ungroup(),mapping=aes(x=UMAP_1,y=UMAP_2,label=EML))+ggtitle("Previous annotation")+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=recheck.lineage.col.set)

plot.results[["CS7.recheck.new.anno"]] <- data.AEP.UMAP %>% ggplot()+geom_point(mapping=aes(x=UMAP_1,y=UMAP_2,col=new_EML))+theme_classic()+NoAxes()+NoLegend()+ggtitle("Previous annotation")+theme(plot.title = element_text(hjust=0.5))+geom_text(data=data.AEP.UMAP %>% group_by(new_EML) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) %>% ungroup(),mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML))+ggtitle("New annotation")+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=recheck.lineage.col.set)


temp.M <- data.AEP.UMAP 
data.temp <- CreateSeuratObject(counts.filter[,temp.M %>% pull(cell)], meta.data =(temp.M  %>% tibble::column_to_rownames("cell")))  %>% NormalizeData(verbose = FALSE)

zs.limit=2.5
temp.mk.list <-  data.AEP.pure.DE$sig %>% group_by(set) %>% filter(power >0.7) %>% top_n(15,power) 
temp.mk.list <- temp.mk.list $gene %>% split(temp.mk.list$set)
temp.exp <- data.temp@assays$RNA@data[unlist(temp.mk.list),]
temp.sel.exp <- t(apply(temp.exp,1,scale))
colnames(temp.sel.exp) <- colnames(temp.exp)
rownames(temp.sel.exp) <- rownames(temp.exp)
temp.sel.exp[temp.sel.exp>zs.limit] <- zs.limit
temp.sel.exp[temp.sel.exp<  (-1*zs.limit)] <- -1*zs.limit

temp.anno <- (temp.M %>% select(cell,EML,new_EML)%>% mutate(od=factor(new_EML,c("Amnion","EPI_Amnion","PriS_Amnion","Epiblast","EPI_PriS","PriS"),ordered=T))%>% arrange(od,EML)  %>% tibble::column_to_rownames("cell"))
#+ fig.show='hide'
plot.results$CS7.recheck.mk.INT <- pheatmap::pheatmap( temp.sel.exp[,rownames(temp.anno )],scale="row",show_rownames=T,show_colnames=F,cluster_cols=F,cluster_rows=F,main="",annotation_col = temp.anno[,c("EML","new_EML")],col=heat.col,annotation_colors = list(EML=recheck.lineage.col.set[unique(temp.anno$EML)],new_EML=recheck.lineage.col.set[unique(temp.anno$new_EML)]),fontsize_row = 5,gaps_row=head(temp.mk.list  %>% lapply(length) %>% cumsum(),2))%>%ggplotify::as.ggplot()

plot.results$CS7.recheck.mk.INT.nolegend <- pheatmap::pheatmap( temp.sel.exp[,rownames(temp.anno )],scale="row",show_rownames=T,show_colnames=F,cluster_cols=F,cluster_rows=F,main="",annotation_col = temp.anno[,c("EML","new_EML")],col=heat.col,annotation_colors = list(EML=recheck.lineage.col.set[unique(temp.anno$EML)],new_EML=recheck.lineage.col.set[unique(temp.anno$new_EML)]),fontsize_row = 5,gaps_row=head(temp.mk.list  %>% lapply(length) %>% cumsum(),2),legend=F,annotation_legend = FALSE)%>%ggplotify::as.ggplot()





#' check the general pvalue distribution
# temp <- query_ref.raw %>% group_by(query_cell,query_EML,query_pj,ref_EML) %>% summarise(pval=min(pval)) %>% ungroup() %>% filter(ref_EML %in% c("Amnion","TE","PriS","Mes","Epiblast","Endoderm","ExE_Mes"))%>% mutate(NLogPval=-log10(pval)) %>% mutate(NLogPval=ifelse(NLogPval>5,5,NLogPval)) %>% select(-pval)%>% spread(ref_EML,NLogPval)  %>% replace(.,is.na(.),0.1) 
# 
# temp %>% GGally::ggpairs(columns=4:10,aes(col=query_pj,alpha=0.5),progress=F,upper=NULL)
# temp %>% GGally::ggpairs(columns=4:10,aes(col=query_pj,alpha=0.5),progress=F,upper=NULL)


#' check the perfomance of model
temp.query_ref.anno.out <- query_ref.anno.out %>% filter(!query_EML %in% c("B1_B2","EarlyBlastocyst")) %>% filter(query_pj %in% c("SPH2016","D3post","CS7","Blakeley","nBGuo","zhou2019"))%>% mutate(query_EML=FunForMatAnno(query_EML),ref_EML=FunForMatAnno(ref_EML)) %>% filter(! query_EML %in% c("ELC","EPI_Amnion","EPI.PrE.INT","EPI.PrE.INT","EPI_PriS","HLC","ICM","ICM-TE_trans","Undefined","TLC","ysTE","PriS_Amnion","Uncertain","8C_Morula"))
temp.stat <- temp.query_ref.anno.out  %>% group_by(query_EML,query_pj) %>% summarise(TPFN=n_distinct(query_cell)) %>% ungroup() %>% left_join(temp.query_ref.anno.out %>% filter(query_EML==ref_EML)%>% group_by(query_EML,query_pj) %>% summarise(TP=n_distinct(query_cell)) %>% ungroup(),by=c("query_EML","query_pj")) %>% left_join(temp.query_ref.anno.out %>% filter(ref_EML!="Uncertain")%>% filter(query_EML!=ref_EML) %>% mutate(query_EML=ref_EML)%>% group_by(query_EML,query_pj) %>% summarise(FP=n_distinct(query_cell)) %>% ungroup(),by=c("query_EML","query_pj")) %>% replace(.,is.na(.),0) %>% mutate(sensitivity=TP/TPFN,specificity=TP/(TP+FP)) 
temp.stat <- temp.stat %>% select(query_EML,query_pj,sensitivity) %>% spread(query_pj,sensitivity) %>% replace(.,is.na(.),0) %>% gather(query_pj,sensitivity,-query_EML) %>% full_join(temp.stat %>% select(query_EML,query_pj,specificity) %>% spread(query_pj,specificity) %>% replace(.,is.na(.),0) %>% gather(query_pj,specificity,-query_EML),by=c("query_EML","query_pj")) %>% semi_join(query_ref.anno.out %>% mutate(query_EML=FunForMatAnno(query_EML)),by = c("query_EML", "query_pj"))

plot.results$model.perf <- temp.stat  %>% ggplot+geom_point(mapping=aes(x=specificity*100,y=sensitivity*100,shape=query_EML,color=query_pj),size=2.5)+theme_classic()+scale_color_manual(values=pj.col[c('SPH2016','D3post','CS7','nBGuo','Blakeley','zhou2019')],labels = c('SPH2016'="Petropoulos et al., 2016",'D3post'="Xiang et al., 2019",'CS7'="Tyser et al., 2020",'nBGuo'="Yanagida et al., 2021",'Blakeley'="Blakeley et al., 2015",'zhou2019'='Zhou et al., 2019'))+scale_shape_manual(values = EML.shape.set) + xlim(0,100)+ylim(0,100)+ggtitle("Performance of model")+theme(plot.title = element_text(hjust=0.5))+xlab("Specificity(%)")+ylab("Sensitivity(%)")

#' check the blastoids model
p="nBGuo"

temp.query <- query_ref.anno.out %>% filter(query_pj==p) %>% filter(!ref_EML %in% c("ambiguous","NoSigHits")) %>% mutate(query_EML=FunForMatAnno(query_EML),ref_EML=FunForMatAnno(ref_EML)) %>% group_by(ref_EML,query_EML) %>% filter(query_EML %in% c("ELC","HLC","TLC"))  %>% summarise(nCell=n_distinct(query_cell))  %>% to_lodes_form(axes=1:2,id="Cohort") %>% mutate(x=factor(x,c("query_EML","ref_EML"),ordered = T)) %>% mutate(stratum=as.vector(stratum)) %>% mutate(stratum=factor(stratum,c("ELC","HLC","TLC","EPI_Amnion","Amnion","Epiblast","Endoderm","TE"),ordered = T)) 

plot.results$nBGuo.model.perf <- ggplot(temp.query,aes(x = x, stratum = stratum, alluvium = Cohort, y = nCell,fill = stratum, label = stratum,col=stratum)) + scale_x_discrete(expand = c(.1, .1)) +geom_flow() +geom_stratum(alpha = .5,lwd=0.1) +geom_text(stat = "stratum", size = 3) +theme(legend.position = "none") +ggtitle("")+theme_void()+NoLegend()+ggtitle(p)+theme(plot.title = element_text(hjust=0.5))+scale_fill_manual(values=lineage.bind.col.set[temp.query$stratum %>% unique() %>% as.vector()])+scale_color_manual(values=lineage.bind.col.set[temp.query$stratum %>% unique() %>% as.vector()])+ggtitle("YAN-Blastoids")


p="IBD2"
temp.query <- query_ref.anno.out %>% filter(query_pj==p) %>% filter(!ref_EML %in% c("ambiguous","NoSigHits")) %>% mutate(query_EML=FunForMatAnno(query_EML),ref_EML=FunForMatAnno(ref_EML)) %>% group_by(ref_EML,query_EML) %>% filter(query_EML %in% c("ELC","HLC","TLC"))  %>% summarise(nCell=n_distinct(query_cell))  %>% to_lodes_form(axes=1:2,id="Cohort") %>% mutate(x=factor(x,c("query_EML","ref_EML"),ordered = T)) %>% mutate(stratum=as.vector(stratum)) %>% mutate(stratum=factor(stratum,c("ELC","HLC","TLC","Epiblast","EPI_PriS","PriS","PriS_Amnion","EPI_Amnion","Amnion","Mes","Endoderm","TE"),ordered = T)) 
plot.results$IBD2.model.perf <- ggplot(temp.query,aes(x = x, stratum = stratum, alluvium = Cohort, y = nCell,fill = stratum, label = stratum,col=stratum)) + scale_x_discrete(expand = c(.1, .1)) +geom_flow() +geom_stratum(alpha = .5,lwd=0.1) +geom_text(stat = "stratum", size = 3) +theme(legend.position = "none") +ggtitle("")+theme_void()+NoLegend()+ggtitle(p)+theme(plot.title = element_text(hjust=0.5))+scale_fill_manual(values=lineage.bind.col.set[temp.query$stratum %>% unique() %>% as.vector()])+scale_color_manual(values=lineage.bind.col.set[temp.query$stratum %>% unique() %>% as.vector()])+ggtitle("IBlastoids")


p="EBD2"
temp.query <- query_ref.anno.out %>% filter(query_pj==p) %>% filter(!ref_EML %in% c("ambiguous","NoSigHits")) %>% mutate(query_EML=FunForMatAnno(query_EML),ref_EML=FunForMatAnno(ref_EML)) %>% group_by(ref_EML,query_EML) %>% filter(query_EML %in% c("ELC","HLC","TLC"))  %>% summarise(nCell=n_distinct(query_cell))  %>% to_lodes_form(axes=1:2,id="Cohort") %>% mutate(x=factor(x,c("query_EML","ref_EML"),ordered = T)) %>% mutate(stratum=as.vector(stratum)) %>% mutate(stratum=factor(stratum,c("ELC","HLC","TLC","Epiblast","EPI_PriS","PriS","PriS_Amnion","Amnion","Endoderm","TE"),ordered = T)) 
plot.results$EBD2.model.perf <- ggplot(temp.query,aes(x = x, stratum = stratum, alluvium = Cohort, y = nCell,fill = stratum, label = stratum,col=stratum)) + scale_x_discrete(expand = c(.1, .1)) +geom_flow() +geom_stratum(alpha = .5,lwd=0.1) +geom_text(stat = "stratum", size = 3) +theme(legend.position = "none") +ggtitle("")+theme_void()+NoLegend()+ggtitle(p)+theme(plot.title = element_text(hjust=0.5))+scale_fill_manual(values=lineage.bind.col.set[temp.query$stratum %>% unique() %>% as.vector()])+scale_color_manual(values=lineage.bind.col.set[temp.query$stratum %>% unique() %>% as.vector()])+ggtitle("Yu-Blastoids")

p="nicolBla"
temp.query <- query_ref.anno.out %>% filter(query_pj==p) %>% filter(!ref_EML %in% c("ambiguous","NoSigHits")) %>% mutate(query_EML=FunForMatAnno(query_EML),ref_EML=FunForMatAnno(ref_EML)) %>% group_by(ref_EML,query_EML) %>% filter(query_EML %in% c("ELC","HLC","TLC"))  %>% summarise(nCell=n_distinct(query_cell))  %>% to_lodes_form(axes=1:2,id="Cohort") %>% mutate(x=factor(x,c("query_EML","ref_EML"),ordered = T)) %>% mutate(stratum=as.vector(stratum)) %>% mutate(stratum=factor(stratum,c("ELC","HLC","TLC","Epiblast","PriS","PriS_Amnion","EPI_Amnion","Amnion","Mes","Endoderm","TE"),ordered = T)) 
plot.results$nicolBla.model.perf <- ggplot(temp.query,aes(x = x, stratum = stratum, alluvium = Cohort, y = nCell,fill = stratum, label = stratum,col=stratum)) + scale_x_discrete(expand = c(.1, .1)) +geom_flow() +geom_stratum(alpha = .5,lwd=0.1) +geom_text(stat = "stratum", size = 3) +theme(legend.position = "none") +ggtitle("")+theme_void()+NoLegend()+ggtitle(p)+theme(plot.title = element_text(hjust=0.5))+scale_fill_manual(values=lineage.bind.col.set[temp.query$stratum %>% unique() %>% as.vector()])+scale_color_manual(values=lineage.bind.col.set[temp.query$stratum %>% unique() %>% as.vector()])+ggtitle("Nicolas-Blastoids")


plot.results$qr.bind.legend <- data.frame(X=(1:(lineage.bind.col.set %>% length())),Y=1:(lineage.bind.col.set %>% length()),anno=names(lineage.bind.col.set)) %>% ggplot+geom_point(mapping=aes(x=X,y=Y,col=anno),shape=15)+scale_color_manual(values=lineage.bind.col.set)+theme_classic()
plot.results$qr.bind.legend <- plot.results$qr.bind.legend %>% ggpubr::get_legend() %>% ggpubr::as_ggplot()



# The detailed cell number with prediction
# print(
#   query_ref.anno.out %>% filter(query_pj %in% c("EBD2","IBD2","nBGuo","nicolBla")) %>% filter(!ref_EML %in% c("ambiguous","NoSigHits")) %>% mutate(query_EML=FunForMatAnno(query_EML),ref_EML=FunForMatAnno(ref_EML))  %>% group_by(query_pj,ref_EML) %>% summarise(nCell=n_distinct(query_cell)) %>% as.data.frame()
# )

# print(
#   query_ref.anno.out %>% filter(query_pj %in% c("EBD2","IBD2","nBGuo","nicolBla")) %>% filter(!ref_EML %in% c("ambiguous","NoSigHits")) %>% mutate(query_EML=FunForMatAnno(query_EML),ref_EML=FunForMatAnno(ref_EML))  %>% group_by(query_EML,query_pj) %>% summarise(nCell=n_distinct(query_cell)) %>% as.data.frame() %>% arrange(query_pj,query_EML)
# )

#' check the Meserderm Marker gene expression in MeLC cells in different blastoids
a="<br>"
#' VlnPlot for Mes marker genes
temp.sel.gene <- c("LIX1","TMEM88","RGS5","PMP22","VIM")#"ANXA1","COL1A1",
for (n in c("IBD2","EBD2","nicolBla")) {
  temp.M <- data.ob.umap %>% filter(pj==n) %>% filter(cluster_EML %in% c("ELC","TLC","HLC","MeLC","AMLC")) %>% select(cell,cluster_EML)
  
  temp.sel.exp <- lognormExp.mBN[temp.sel.gene ,temp.M$cell]%>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,logExp,-Gene) %>% inner_join(temp.M,by="cell") %>% mutate(od=factor(cluster_EML,c("MeLC","ELC","HLC","TLC","AMLC"),ordered = T)) %>% arrange(od)  
  
  for (g in temp.sel.gene ) {
    plot.results[[paste0("Mes.MK.Vln.",n,".",g)]] <- temp.sel.exp%>% filter(Gene==g) %>% ggplot(mapping=aes(x=od,y=logExp)) + geom_violin(mapping=aes(fill=cluster_EML),color="black",scale = "width")+geom_jitter(size=0.05)+ theme_classic() + theme(axis.text.x=element_text(angle = 40,hjust = 1)) +NoLegend()+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="bold"))+ylab("Expression level")+xlab("")+ylab("")+scale_fill_manual(values=lineage.col.set[c("ELC","TLC","HLC","MeLC","AMLC")])
    plot.results.jpeg[[paste0("Mes.MK.Vln.",n,".",g)]] <- temp.sel.exp%>% filter(Gene==g) %>% ggplot(mapping=aes(x=od,y=logExp)) + geom_violin(mapping=aes(fill=anno),color="black",scale = "width")+ theme_classic() + theme(axis.text.x=element_text(angle = 40,hjust = 1)) +NoLegend()+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="bold"))+ylab("Expression level")+xlab("")+ylab("")+scale_fill_manual(values=lineage.col.set[c("ELC","TLC","HLC","MeLC","AMLC")])
  }
}

save(plot.results.jpeg,plot.results.txt,plot.results,file=paste0("tmp_data/",TD,"/figplot.Rdata"))


for (n in names(plot.results)) {
  print(plot.results[[n]])
}


