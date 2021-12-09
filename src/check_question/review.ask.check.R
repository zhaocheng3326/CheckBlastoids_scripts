#' R4.0
rm(list=ls())

condaENV <- "/home/chenzh/miniconda3/envs/R4.0"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(pheatmap))
suppressMessages(library(Seurat))


# working directory
DIR <- "~/My_project/CheckBlastoids"
setwd(DIR)

#' Loading R functions
source("~/PC/R_code/functions.R")
source("~/PC/SnkM/SgCell.R")
source("loca.quick.fun.R")
source("figures.setting.R")


#' loading options
TD="Nov14_2021"
MM="Pub"


#' loading data
load("tmp_data/gene.meta.Rdata",verbose=T)
meta.filter<- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% filter(pj %in% c("EBD2","IBD2","JPF2019","D3post","CS7","SPH2016","nBGuo","Blakeley") ) %>% mutate(EML=ifelse(EML %in% c("Ectoderm"),"Amnion",EML)) %>% mutate(EML=ifelse(EML=="EB_UI13","EB_ELC",EML))

#' counts expression
counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))

#' loading normalized expression
lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/",MM,".lognormExp.mBN.rds"))

#' loading Amnion vs TE DEGs
hm.AMvsTE.out <- readRDS(paste0("tmp_data/",TD,"/human.AmnionVsTE.sigMarker.rds"))
lineage.mk <- readRDS(paste0("tmp_data/",TD,"/",MM,".mk.allL.rds"))
#' loading all lineage markers

#' loading whole human integration results
data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/Pub.D2.umap.cord.rds")) %>%  mutate(rename_EML=ifelse(EML=="EB_UI13","ELC",rename_EML)) %>% mutate(EML=ifelse(EML=="EB_UI13","EB_ELC",EML))
data.ob <- readRDS(paste0("tmp_data/",TD,"/Pub.D2.mnn.rds"))

#' loading cross-species integration results
data.umap.list <- readRDS(paste0("tmp_data/",TD,"/NHP.IT.umap.rds"))


#' loading prediction results
query_ref.anno.out <- readRDS(paste0("tmp_data/",TD,"/prediction.model.performance.rds")) 


#' type of related cells
cells.list <- list()

cells.list$Blakeley_TE <- meta.filter %>% filter(pj=="Blakeley" & EML=="TE")  %>% pull(cell)
cells.list$SPH2016_TE <- meta.filter %>% filter(pj=="SPH2016" & EML=="TE")  %>% pull(cell)
cells.list$D3post_TE <- meta.filter %>% filter(pj=="D3post" & EML %in% c("CTB","STB","EVT"))  %>% pull(cell)
cells.list$D3post_TE_pre <- meta.filter %>% filter(pj=="D3post" & EML %in% c("CTB","STB","EVT") & devTime %in% c("E6","E7"))  %>% pull(cell)
cells.list$D3post_TE_post <- meta.filter %>% filter(pj=="D3post" & EML %in% c("CTB","STB","EVT") & devTime %in% c("E8","E9","E10","E12","E14"))  %>% pull(cell)
cells.list$nBGuo_TE <- meta.filter %>% filter(pj=="nBGuo" & EML %in% c("EarlyTE","TE"))  %>% pull(cell)
cells.list$CS7_Amnion <- meta.filter %>% filter(pj=="CS7" & EML %in% c("Amnion")) %>% pull(cell)

cells.list$IBD2_TLC <-  meta.filter %>% filter(pj=="IBD2" & EML %in% c("IB_TE")) %>% pull(cell)
cells.list$EBD2_TLC <-  meta.filter %>% filter(pj=="EBD2" & EML %in% c("EB_TLC")) %>% pull(cell)
cells.list$EBD2_TLC.mod <- data.ob.umap %>% filter(pj=="EBD2" & cluster_EML %in% c("TLC")) %>% pull(cell) # purified EB_TLC
cells.list$nBGuo_TLC <-  meta.filter %>% filter(pj=="nBGuo" & EML %in% c("TLC")) %>% pull(cell)
cells.list$JPF2019_AMLC <-  meta.filter %>% filter(pj=="JPF2019" & EML %in% c("AMLC")) %>% pull(cell)
cells.list$JPF2019_Tsw_AMLC <-  meta.filter %>% filter(pj=="JPF2019" & EML %in% c("Tsw-AMLC")) %>% pull(cell)

#" check ISL1 and GABRP
for (g in c("ISL1","GABRP")) {
  for (n in names(cells.list)) {
    
    a1=length(cells.list[[n]])
    b1=sum(counts.filter[g,cells.list[[n]]] >0 )
    r1=round(b1/a1,3)
    print(paste(n,g,b1,a1,r1))
  }
  
}

for (g in c("ISL1","GABRP")) {
  print(
    fisher.test(
      matrix(c(
        sum(counts.filter[g,cells.list[["IBD2_TLC"]]] >0),
        length(cells.list[["IBD2_TLC"]])-sum(counts.filter[g,cells.list[["IBD2_TLC"]]] >0),
        sum(counts.filter[g,cells.list[["SPH2016_TE"]]] >0),
        length(cells.list[["SPH2016_TE"]])-sum(counts.filter[g,cells.list[["SPH2016_TE"]]] >0)
      ),nrow = 2)
    )
  )
}


#' check the "Mes" of Iblastoids
IBD2.MeLC <- data.ob.umap %>% filter(cluster_EML=="MeLC" & pj=="IBD2") %>% pull(cell)
IBD2.IM1 <-  data.ob.umap %>% filter(EML=="IB_IM1" & pj=="IBD2") %>% pull(cell)
IBD2.IM2 <-  data.ob.umap %>% filter(EML=="IB_IM2" & pj=="IBD2") %>% pull(cell)

data.ob.umap %>% mutate(selected=ifelse(cell %in% IBD2.MeLC,"sel","no")) %>% mutate(selected=factor(selected,c("no","sel",ordered=T))) %>% arrange(selected) %>% ggplot()+geom_point(mapping=aes(x=UMAP_1,y=UMAP_2,col=selected),size=0.25)+theme_void()+scale_color_manual(values=c("no"="grey","sel"="red"))+ theme(legend.position = "none")+ggtitle("MeLC from Iblastoids")+theme(plot.title = element_text(hjust=0.5,size=15,face="bold"))

data.umap.list$NHP.IBD2.mnn %>% mutate(selected=ifelse(cell %in% IBD2.MeLC,"sel","no")) %>% mutate(selected=factor(selected,c("no","sel",ordered=T))) %>% arrange(selected) %>% ggplot()+geom_point(mapping=aes(x=UMAP_1,y=UMAP_2,col=selected),size=0.25)+theme_void()+scale_color_manual(values=c("no"="grey","sel"="red"))+ theme(legend.position = "none")+ggtitle("MeLC from Iblastoids")+theme(plot.title = element_text(hjust=0.5,size=15,face="bold"))+ geom_text(data.umap.list$NHP.IBD2.mnn %>% group_by(EML,pj) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)),mapping=aes(x=UMAP_1,y=UMAP_2,label=EML))

data.umap.list$NHP.IBD2.mnn %>% mutate(selected=ifelse(cell %in% IBD2.IM1,"sel","no")) %>% mutate(selected=factor(selected,c("no","sel",ordered=T))) %>% arrange(selected) %>% ggplot()+geom_point(mapping=aes(x=UMAP_1,y=UMAP_2,col=selected),size=0.25)+theme_void()+scale_color_manual(values=c("no"="grey","sel"="red"))+ theme(legend.position = "none")+ggtitle("MeLC from Iblastoids")+theme(plot.title = element_text(hjust=0.5,size=15,face="bold"))+ geom_text(data.umap.list$NHP.IBD2.mnn %>% group_by(EML,pj) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)),mapping=aes(x=UMAP_1,y=UMAP_2,label=EML))

data.umap.list$NHP.CS7.mnn %>% mutate(selected=ifelse(EML %in% (c("Mes1","Mes2")),"sel","no")) %>% mutate(selected=factor(selected,c("no","sel",ordered=T))) %>% arrange(selected) %>% ggplot()+geom_point(mapping=aes(x=UMAP_1,y=UMAP_2,col=selected),size=0.25)+theme_void()+scale_color_manual(values=c("no"="grey","sel"="red"))+ theme(legend.position = "none")+theme(plot.title = element_text(hjust=0.5,size=15,face="bold"))+ geom_text(data.umap.list$NHP.CS7.mnn %>% group_by(EML,pj) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)),mapping=aes(x=UMAP_1,y=UMAP_2,label=EML))+ggtitle("Mes1 and Mes2 from NHP dataset")


a=read.delim("/home/chenzh/My_project/CheckBlastoids/data/GSE136447/EPI_vs_AME.tsv",stringsAsFactors = F,head=F) %>% tbl_df()


ct <- 0
table(counts.filter[c("ISL1","GABRP","GATA2","GATA3"),meta.filter  %>% filter(pj=="IBD2" & EML=="IB_TE") %>% pull(cell)] %>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,counts,-Gene) %>% filter(Gene=="GATA2") %>% pull(counts) > ct)
table(counts.filter[c("ISL1","GABRP","GATA2","GATA3"),meta.filter  %>% filter(pj=="IBD2" & EML=="IB_TE") %>% pull(cell)] %>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,counts,-Gene) %>% filter(Gene=="GATA3") %>% pull(counts) > ct)

table(counts.filter[c("ISL1","GABRP","GATA2","GATA3"),meta.filter  %>% filter(pj=="IBD2" & EML=="IB_TE") %>% pull(cell)] %>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,counts,-Gene) %>% filter(Gene=="ISL1") %>% pull(counts) > ct)
table(counts.filter[c("ISL1","GABRP","GATA2","GATA3"),meta.filter  %>% filter(pj=="IBD2" & EML=="IB_TE") %>% pull(cell)] %>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,counts,-Gene) %>% filter(Gene=="GABRP") %>% pull(counts) > ct)
#1234



ct <- 0
table(counts.filter[c("ISL1","GABRP","GATA2","GATA3"),data.ob.umap %>% filter(new_EML=="TE" & pj %in% c("SPH2016","D3post","CS7")) %>% pull(cell)] %>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,counts,-Gene) %>% filter(Gene=="GATA2") %>% pull(counts) > ct)

table(counts.filter[c("ISL1","GABRP","GATA2","GATA3"),data.ob.umap %>% filter(new_EML=="TE" & pj %in% c("SPH2016","D3post","CS7")) %>% pull(cell)] %>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,counts,-Gene) %>% filter(Gene=="GATA3") %>% pull(counts) > ct)

table(counts.filter[c("ISL1","GABRP","GATA2","GATA3"),data.ob.umap %>% filter(new_EML=="TE" & pj %in% c("SPH2016","D3post","CS7")) %>% pull(cell)] %>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,counts,-Gene) %>% filter(Gene=="ISL1") %>% pull(counts) > ct)


table(counts.filter[c("ISL1","GABRP","GATA2","GATA3"),data.ob.umap %>% filter(new_EML=="TE" & pj %in% c("SPH2016","D3post","CS7")) %>% pull(cell)] %>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,counts,-Gene) %>% filter(Gene=="GABRP") %>% pull(counts) > ct)



#temp.sel.gene <- c("ISL1","GABRP","IGFBP2","GATA2","GTSF1","TGFBR3") # from the dot plot 

temp.sel.gene <- c("TNMD","GABRP","ITGB6","HOXD9","VTCN1","KRT24","MUC16","TNC", "VIM", "VTCN1","HLA-A","HLA-B","GATA2","GATA3")
temp.M <- data.ob.umap %>% filter(EML %in% c("TE","IB_TE","EB_TLC","STB","EVT","CTB","Amnion","Tsw-AMLC","AMLC"))
temp.M.smt2 <- temp.M %>% filter(seqType=="smt2")
temp.M.10X <- temp.M %>% filter(seqType=="10X")


temp.sel.exp <- lognormExp.mBN[temp.sel.gene ,c(temp.M.10X$cell,temp.M.smt2$cell)]%>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,logExp,-Gene) %>% inner_join(data.ob.umap %>% select(cell,seqType,EML,cluster_EML) %>% mutate( bigEML= ifelse(EML %in% c("TE","STB","EVT","CTB"),'TE',EML)),by="cell") %>% mutate(od=factor(bigEML,c("Amnion","TE","AMLC","Tsw-AMLC","EB_TLC","IB_TE"),label=c("Amnion","TE","AMLC","Tsw-AMLC","TLC(SB)","TLC(IB)"),ordered = T)) %>% mutate(EML=ifelse(EML=="EB_TLC","TLC(SB)",EML))%>% mutate(EML=ifelse(EML=="IB_TE","TLC(IB)",EML))%>% mutate(seqType=factor(seqType, levels=c('smt2','10X'))) %>% arrange(od)  %>% mutate(EML=ifelse(EML%in% c("TE","STB","EVT","CTB"),'TE',EML))


# In a recent report by Io et al., Cell Stem Cell 2021, they report some genes that seem to be enriched in the amnion: TNC, VIM, VTCN1 (Fig. 6). Though they have not yet been validated at a protein level on human in vivo amnion so need to be taken with caution, are these genes also expressed in the populations claimed to be amnion? Would it be helpful to include/compare with the amnion dataset used in his paper though they are based on the data from the cynomolgus monkey?
#In the same report by Io et al., HLA Class I expression is absent in TE cells (Fig. S5D), similarly to the first trimester trophoblast (Lee et al., Stem cell Reports 2016). Could the expression of HLA-A,B be examined in the amnion and iBlastoid datasets?
# In Fig. 1c, the expression levels of ISL1 and GABRP are too sparse to be considered as being expressed in TLC(IB). What is the percentage of TLC(IB) that expressed these two genes at a reasonable expression threshold? Is their expression in TLC(IB) higher than expected by chance? # 

pdf("temp.pdf")
for (g in temp.sel.gene ) {
  print(
    temp.sel.exp%>% filter(Gene==g) %>% ggplot(mapping=aes(x=od,y=logExp)) + geom_violin(mapping=aes(fill=EML),color="black",scale = "width")+geom_jitter(size=0.01)+ theme_classic() + theme(axis.text.x=element_text(angle = 40,hjust = 1)) +NoLegend()+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="bold"))+ylab("Expression level")+xlab("")+facet_grid(~seqType,scale="free_x")+theme(strip.background = element_blank(),strip.text.x = element_blank())+scale_fill_manual(values=vln.col.set)+theme(axis.text.x = element_text(size = 6))+ylab("")
  )
}
dev.off()

#' some extra check 
#' ### VlnPlot for key marker genes
#temp.sel.gene <- c("TNMD","GABRP","ITGB6","HOXD9","VTCN1","KRT24","MUC16")
temp.sel.gene <- c("DNMT3L","DPPA3","GTSF1","CA3","RARRES2")
temp.M <- data.ob.umap %>% filter(EML!="Tsw-AMLC")%>% filter(pj !="D3post") %>% select(cell,rename_EML,pj) %>% filter(rename_EML %in% c("Amnion","TE","TLC","AMLC")) %>% mutate(anno=paste(pj,rename_EML,sep="_")) %>% select(cell,anno) %>% bind_rows(data.ob.umap %>% filter(pj =="D3post") %>% filter(rename_EML %in% "TE") %>% mutate(anno=ifelse(devTime %in% c("E6","E7"),"D3_PreTE","D3_PostTE")) %>% select(cell,anno)) %>% bind_rows(meta.filter %>% filter(EML=="TE" & pj=="Blakeley") %>% select(cell) %>% mutate(anno="Blakeley_TE"))

temp.sel.exp <- lognormExp.mBN[temp.sel.gene ,temp.M$cell]%>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,logExp,-Gene) %>% inner_join(temp.M,by="cell") %>% mutate(od=factor(anno,c("CS7_Amnion","Blakeley_TE","SPH2016_TE","nBGuo_TE","D3_PreTE","D3_PostTE","JPF2019_AMLC","IBD2_TLC","EBD2_TLC","nicolBla_TLC","nBGuo_TLC"),ordered = T)) %>% arrange(od)  

#"CS7_Amnion","Blakeley_TE","SPH2016_TE","nBGuo_TE","D3_PreTE","D3_PostTE"
temp.list <- list()
for (g in temp.sel.gene ) {
  temp.list[[g]] <- temp.sel.exp%>% filter(Gene==g) %>% ggplot(mapping=aes(x=od,y=logExp)) + geom_violin(mapping=aes(fill=anno),color="black",scale = "width")+geom_jitter(size=0.05)+ theme_classic() + theme(axis.text.x=element_text(angle = 40,hjust = 1)) +NoLegend()+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="bold"))+ylab("Expression level")+xlab("")+ylab("")

}
cowplot::plot_grid(plotlist=temp.list)
