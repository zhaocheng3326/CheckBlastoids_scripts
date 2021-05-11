#' ---
#' title: "Generating figures"
#' output: 
#'  html_document:
#'    code_folding: hide
#'    
#' ---


#  R3.6
rm(list=ls())

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(patchwork))
suppressMessages(library(cowplot))
suppressMessages(library(pheatmap))
suppressMessages(library(Seurat))

#suppressMessages(library(scran))
#suppressMessages(library(batchelor))
#suppressMessages(library(readxl))

# working directory
DIR <- "~/My_project/JP_project"
setwd(DIR)

#' Loading R functions
source("~/PC/R_code/functions.R")
source("~/PC/SnkM/SgCell.R")
source("loca.quick.fun.R")

suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)

#' #### setting 
options(digits = 4)
options(future.globals.maxSize= 3001289600)
TD="D2_pub"
zs.limit <- 2.5
heat.col <- colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100)
rename <- dplyr::rename
select<- dplyr::select
filter <- dplyr::filter

lineage.col.set <-c("Amnion"="#F8766D" ,"AMLC"="#F8766D" , "ELC"="#00BADE" , "Epiblast"="#00BADE","Endoderm"="#B385FF","HLC"="#B385FF","TE"= "#64B200", "TLC"= "#64B200","MeLC"= "#EF67EB","Mesoderm"= "#EF67EB", "Prelineage"= "#7B554E", "PriS"="#DB8E00",  "Undef"= "grey50", "ExE_Mes"="#C4423E","hES"="#51C0CC","hPGCLC"="#CD8BBC","Stem-blastoids"="#F8766D","Iblastoids"="#00BA38","Embryonic"="royalblue3","ICM"="grey50","Stem-blastoids"="grey33","Iblastoids"="grey33")

heat.col <- colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100)

vln.col.set <- c("TE"="#DA8F00","CTB"= "#86AC00","EVT"="#00C094","STB"= "#00B6EB" ,"Amnion"="#F8766D","TLC(SB)"="#DA8F00","TLC(IB)"="#FF6B96","AMLC"= "#F8766D","Tsw-AMLC"= "#F8766D","ELC"="#00BADE" ,"TLC"= "#64B200","HLC"="#B385FF","MeLC"= "#EF67EB")

cluster.col.set <- c("C0"="#00C094","C1"="#E38900","C2"="#99A800","C3"="#00BFC4","C4"="#00BC56","C5"="#F8766D","C6"="#53B400","C7"="#C49A00","C8"="#06A4FF","C9"="#A58AFF","C10"="#DF70F8","C11"="#FB61D7","C12"="#FF66A8","Others"="grey66")

pj.shape.set <- c(CS7=21,D3post=22,SPH2016=24)
pj.shape.set.solid <- c(CS7=16,D3post=15,SPH2016=17)

label.pj <- list()
label.pj$SPH2016="Petropoulos et al., 2016"
label.pj$D3post="Xiang et al., 2020"
label.pj$CS7="Tyser et al., 2021"
label.pj$JPF2019="Zheng et al., 2019"


#' #### loading data
load("tmp_data/gene.meta.Rdata",verbose=T)
# expression
lognormExp.mBN <- readRDS(paste0("tmp_data/",TD,"/lognormExp.mBN.rds"))
# meta data
meta.filter<- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% filter(pj %in% c("EBD2","IBD2","JPF2019","D3post","CS7","SPH2016") ) %>% mutate(EML=ifelse(EML %in% c("Ectoderm"),"Amnion",EML))
meta.filter.ds<- readRDS(paste0("tmp_data/",TD,"/meta.filter.ds.rds"))%>% filter(pj %in% c("EBD2","IBD2","JPF2019","D3post","CS7","SPH2016") ) 
# #### loading top DE among AMN, TE, Epiblast and End
load(paste0("tmp_data/",TD,"/DE.allL.Rdata"))
# #### loading umap result
data.ob.umap <- readRDS(paste0("tmp_data/",TD,"/Pub.D2.umap.cord.rds"))



#' text position embryonic dataset
data.EM.umap.text.pos <- data.ob.umap  %>% filter(pj %in% c("D3post","CS7","SPH2016"))  %>%group_by(new_EML) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) %>% ungroup()  %>% filter(!new_EML %in% c("3D_ICM","PSA-EPI")) %>% mutate(UMAP_1=ifelse(new_EML=="PriS",UMAP_1-1,UMAP_1)) %>% mutate(new_EML=gsub("ExE_Mes","Extraembryonic\nMesoderm",new_EML))%>% mutate(new_EML=gsub("PriS","Primitive\nStreak",new_EML))

#' text position of other datasets
data.LC.umap.text.pos <- data.ob.umap  %>% filter(pj %in% c("EBD2","IBD2","JPF2019"))  %>%group_by(new_EML,pj) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) %>% ungroup() %>% mutate(UMAP_2=ifelse(new_EML=="HLC",UMAP_2+0.5,UMAP_2)) # UMAP_2 of HLC adding 0.5 to avoid overlapping



Plot.umap <- list() 
Plot.umap$em <- ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("D3post","CS7","SPH2016")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.5)+geom_point( data.ob.umap %>% filter(pj %in% c("D3post","CS7","SPH2016"))  %>% mutate(od=factor(new_EML,c("Prelineage","TE","Mesoderm","Epiblast","Amnion","Endoderm","PriS","ExE_Mes","3D_ICM","PSA-EPI"),ordered = T)) %>% arrange(desc(od)) ,mapping=aes(x=UMAP_1,y=UMAP_2,fill=new_EML,shape=pj),size=2,color="grey")+ggtitle("Embryonic dataset")+geom_text(data.EM.umap.text.pos ,mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=4,color="Black",fontface="bold")+theme_classic()+theme(plot.title = element_text(hjust=0.5))+scale_fill_manual(values=lineage.col.set)+scale_shape_manual(values = pj.shape.set,labels = c('SPH2016'="Petropoulos et al., 2016",'D3post'="Xiang et al., 2019",'CS7'="Tyser et al., 2020"))+NoAxes()+NoLegend()

Plot.umap$em_legend <- ggplot()+ geom_point( data.ob.umap %>% filter(pj %in% c("D3post","CS7","SPH2016")) %>% filter(!new_EML %in% c("PSA-EPI","3D_ICM")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=new_EML,shape=pj),size=2)+theme_classic()+scale_color_manual("",values=lineage.col.set)+scale_shape_manual("",values = pj.shape.set.solid,labels = c('SPH2016'="Petropoulos et al., 2016",'D3post'="Xiang et al., 2020",'CS7'="Tyser et al., 2020")) + theme(legend.text = element_text( size = 6),legend.title = element_text( size = 1))
Plot.umap$em_legend <- Plot.umap$em_legend%>% ggpubr::get_legend() %>% ggpubr::as_ggplot()

umap.xlim=ggplot_build(Plot.umap$em)$layout$panel_scales_x[[1]]$range$range
umap.ylim=ggplot_build(Plot.umap$em)$layout$panel_scales_y[[1]]$range$range
Plot.umap$em.tiff1 <- ggplot()+ ggtitle("Embryonic dataset")+geom_text(data.EM.umap.text.pos ,mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=4,color="Black",fontface="bold")+theme_classic()+theme(plot.title = element_text(hjust=0.5))+scale_fill_manual(values=lineage.col.set)+scale_shape_manual(values = pj.shape.set,labels = c('SPH2016'="Petropoulos et al., 2016",'D3post'="Xiang et al., 2019",'CS7'="Tyser et al., 2020"))+NoAxes()+NoLegend()+xlim(umap.xlim)+ylim(umap.ylim)
Plot.umap$em.tiff2 <- ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("D3post","CS7","SPH2016")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.5)+geom_point( data.ob.umap %>% filter(pj %in% c("D3post","CS7","SPH2016"))  %>% mutate(od=factor(new_EML,c("Prelineage","TE","Mesoderm","Epiblast","Amnion","Endoderm","PriS","ExE_Mes","3D_ICM","PSA-EPI"),ordered = T)) %>% arrange(desc(od)) ,mapping=aes(x=UMAP_1,y=UMAP_2,fill=new_EML,shape=pj),size=2,color="grey")+ggtitle("")+theme_classic()+theme(plot.title = element_text(hjust=0.5))+scale_fill_manual(values=lineage.col.set)+scale_shape_manual(values = pj.shape.set,labels = c('SPH2016'="Petropoulos et al., 2016",'D3post'="Xiang et al., 2019",'CS7'="Tyser et al., 2020"))+NoAxes()+NoLegend()+xlim(umap.xlim)+ylim(umap.ylim)


Plot.umap$JPF <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("JPF2019")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("JPF2019")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=new_EML),size=0.5)+ggtitle("Zheng et al., 2019")+geom_text(data.EM.umap.text.pos %>% filter(new_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=3.5,color="Black")+  geom_text(data.LC.umap.text.pos %>% filter(pj=="JPF2019"),mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=4,color="Black",fontface="bold")  +theme_classic()+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=lineage.col.set)+xlim(umap.xlim)+ylim(umap.ylim)
Plot.umap$JPF.tiff1 <-  ggplot()+ggtitle("Zheng et al., 2019")+geom_text(data.EM.umap.text.pos %>% filter(new_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=3.5,color="Black")+  geom_text(data.LC.umap.text.pos %>% filter(pj=="JPF2019"),mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=4,color="Black",fontface="bold")  +theme_classic()+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=lineage.col.set)+xlim(umap.xlim)+ylim(umap.ylim)
Plot.umap$JPF.tiff2 <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("JPF2019")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("JPF2019")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=new_EML),size=0.5)+ggtitle("") +theme_classic()+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=lineage.col.set)+xlim(umap.xlim)+ylim(umap.ylim)


Plot.umap$EBD2 <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("EBD2")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("EBD2"))%>% mutate(od=factor(new_EML,c("ELC","HLC","TLC","Undef"),ordered = T)) %>% arrange(desc(od)) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=new_EML),size=0.5)+ggtitle("Stem-blastoids")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+geom_text(data.EM.umap.text.pos %>% filter(new_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=3.5,color="Black")+geom_text(data.LC.umap.text.pos %>% filter(pj=="EBD2") %>% filter(new_EML %in% c("ELC","HLC","TLC")),mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=4,color="Black",fontface="bold")+xlim(umap.xlim)+ylim(umap.ylim)
Plot.umap$EBD2.tiff1 <-  ggplot()+ggtitle("Stem-blastoids")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+geom_text(data.EM.umap.text.pos %>% filter(new_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=3.5,color="Black")+geom_text(data.LC.umap.text.pos %>% filter(pj=="EBD2") %>% filter(new_EML %in% c("ELC","HLC","TLC")),mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=4,color="Black",fontface="bold")+xlim(umap.xlim)+ylim(umap.ylim)
Plot.umap$EBD2.tiff2 <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("EBD2")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("EBD2"))%>% mutate(od=factor(new_EML,c("ELC","HLC","TLC","Undef"),ordered = T)) %>% arrange(desc(od)) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=new_EML),size=0.5)+ggtitle("")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+xlim(umap.xlim)+ylim(umap.ylim)

Plot.umap$IBD2 <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("IBD2")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("IBD2"))  %>% mutate(od=factor(new_EML,c("ELC","HLC","TLC","Undef"),ordered = T)) %>% arrange(desc(od)),mapping=aes(x=UMAP_1,y=UMAP_2,col=new_EML),size=0.5)+ggtitle("Iblastoids")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+geom_text(data.EM.umap.text.pos %>% filter(new_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=3.5,color="Black")+geom_text(data.LC.umap.text.pos %>% filter(pj=="IBD2") %>% filter(new_EML %in% c("ELC","HLC","TLC")),mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=4,color="Black",fontface="bold")+xlim(umap.xlim)+ylim(umap.ylim)
Plot.umap$IBD2.tiff1 <-  ggplot()+ggtitle("Iblastoids")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+geom_text(data.EM.umap.text.pos %>% filter(new_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=3.5,color="Black")+geom_text(data.LC.umap.text.pos %>% filter(pj=="IBD2") %>% filter(new_EML %in% c("ELC","HLC","TLC")),mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=4,color="Black",fontface="bold")+xlim(umap.xlim)+ylim(umap.ylim)
Plot.umap$IBD2.tiff2 <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("IBD2")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("IBD2"))  %>% mutate(od=factor(new_EML,c("ELC","HLC","TLC","Undef"),ordered = T)) %>% arrange(desc(od)),mapping=aes(x=UMAP_1,y=UMAP_2,col=new_EML),size=0.5)+ggtitle("")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+xlim(umap.xlim)+ylim(umap.ylim)


pdf("D2.UMAP.plot.pdf",7,6)
(
  plot_spacer()+Plot.umap$em+Plot.umap$em_legend+plot_spacer()+plot_layout(widths = c(0.5,2,1,0.5))
)/(
  (Plot.umap$JPF|Plot.umap$EBD2|Plot.umap$IBD2)
  
)+plot_layout(heights = c(2,1.5))
dev.off()
pdf("D2.UMAP.plot.text.pdf",7,6)
(
  plot_spacer()+Plot.umap$em.tiff1+Plot.umap$em_legend+plot_spacer()+plot_layout(widths = c(0.5,2,1,0.5))
)/(
  (Plot.umap$JPF.tiff1|Plot.umap$EBD2.tiff1|Plot.umap$IBD2.tiff1)
  
)+plot_layout(heights = c(2,1.5))
dev.off()
pdf("D2.UMAP.plot.point.pdf",7,6)
(
  plot_spacer()+Plot.umap$em.tiff2+plot_spacer()+plot_spacer()+plot_layout(widths = c(0.5,2,1,0.5))
)/(
  (Plot.umap$JPF.tiff2|Plot.umap$EBD2.tiff2|Plot.umap$IBD2.tiff2)
  
)+plot_layout(heights = c(2,1.5))
dev.off()

#' #### Greate the heatmap and dotplot 
temp.M <-  data.ob.umap %>% filter(new_EML %in% c("Epiblast","ELC","TE","TLC","Endoderm","HLC","Mesoderm","MeLC","Amnion","AMLC"))
temp.M.smt2 <- temp.M %>% filter(seqType=="smt2")
temp.M.10X <- temp.M %>% filter(seqType=="10X")
de.gene.sel <- temp.DM$sig %>% group_by(set)%>% filter(power >0.6) %>% top_n(10,power) %>%  dplyr::rename(Gene=gene) %>% mutate(od=factor(set,c("Epiblast","Endoderm","TE","Mesoderm","Amnion"),ordered = T)) %>% arrange(od,desc(power)) %>% select(-od)


temp <- lognormExp.mBN[de.gene.sel$Gene ,temp.M.smt2$cell]  %>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,logExp,-Gene) %>% inner_join(temp.M.smt2 %>% select(cell,new_EML) %>% dplyr::rename(cell.cluster=new_EML),by="cell") %>% inner_join(de.gene.sel%>% select(Gene,set) %>% dplyr::rename(Gene.cluster=set) ,by="Gene") 


temp.input <- temp %>% group_by(cell.cluster,Gene.cluster,Gene) %>% summarise(meanlogExp=mean(logExp),nCell=n_distinct(cell)) %>% left_join(temp %>% group_by(cell.cluster,Gene.cluster,Gene) %>% filter(logExp >0)%>% summarise(nExpCell=n_distinct(cell)),by=c("cell.cluster","Gene.cluster","Gene"))  %>% mutate(nExpCell=ifelse(is.na(nExpCell),0,nExpCell))%>% mutate(nExp.ct=nExpCell/nCell) %>% group_by(Gene) %>% mutate(sv=(meanlogExp-mean(meanlogExp))/sd(meanlogExp)) %>% mutate(nExp.ct= as.vector(cut(nExp.ct*100,c(-1,5,25,50,75,100),label=c(0,25,50,75,100)))) %>% mutate(nExp.ct=as.numeric(nExp.ct))

Plot.Exp <- list()
Plot.Exp$dot <- temp.input %>% mutate(Gene=factor(Gene,rev(de.gene.sel$Gene),ordered = T)) %>% mutate(cell.cluster=factor(cell.cluster,c("Epiblast","Endoderm","TE","Mesoderm","Amnion"),ordered = T))%>% mutate(Gene.cluster=factor(Gene.cluster,c("Epiblast","Endoderm","TE","Mesoderm","Amnion"),ordered = T)) %>%  ggplot(aes(x=cell.cluster, y = Gene, color =sv, size = nExp.ct)) + geom_point() +facet_grid(Gene.cluster~.,,scale="free_y") + scale_color_viridis_c(name = '') + cowplot::theme_cowplot() + theme(axis.line  = element_blank()) +theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +ylab('') +theme(axis.ticks = element_blank()) +scale_size(name="",breaks = c(0,25,50,75,100),range = c(0.1, 1.5))+xlab("")+theme(axis.text.y = element_text(size = 5))+theme(axis.text.x = element_text(size = 7))+ theme(legend.position="top")+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(strip.background=element_blank(),strip.text.y=element_blank())+NoLegend()+scale_y_discrete(position = "right")

Plot.Exp$dot_legend <-temp.input  %>%  ggplot(aes(x=cell.cluster, y = Gene, color =sv, size = nExp.ct)) + geom_point()  + scale_color_viridis_c(name = '') + cowplot::theme_cowplot()  +scale_size(name="",breaks = c(0,25,50,75,100),range = c(0.1, 1.5))+ theme(legend.position="top")+theme(legend.text = element_text( size = 6),legend.box="vertical", legend.margin=margin())  #Scaled Expression  % cells expressing gene
Plot.Exp$dot_legend <- Plot.Exp$dot_legend%>% ggpubr::get_legend() %>% ggpubr::as_ggplot()

# heatmap for embryonic cells
temp.sel.cell <- temp.M.smt2%>% pull(cell)
temp.exp <- lognormExp.mBN[de.gene.sel$Gene ,temp.sel.cell] 
temp.sel.exp <- t(apply(temp.exp,1,scale))
colnames(temp.sel.exp) <- colnames(temp.exp)
rownames(temp.sel.exp) <- rownames(temp.exp)
temp.sel.exp[temp.sel.exp>zs.limit] <- zs.limit
temp.sel.exp[temp.sel.exp<  (-1*zs.limit)] <- -1*zs.limit
gaps <- 10

temp.anno <-temp.M %>% filter(cell %in% temp.sel.cell) %>% select(cell,new_EML,pj)  %>% dplyr::rename(cluster=new_EML)%>% arrange(cluster) %>% mutate(od=factor(cluster,c("Epiblast","Endoderm","TE","Mesoderm","Amnion"),ordered = T)) %>% arrange(od) %>% select(-od)  %>% tibble::column_to_rownames("cell") 
pheatmap(temp.sel.exp[,rownames(temp.anno)],cluster_rows=F,,cluster_cols=F,scale="none",annotation_col=temp.anno[,c("cluster","pj"),drop=F],show_colnames=F,show_rownames=T,color=heat.col, fontsize_row=7,border_color="NA",main="\nEmbryonic cells",gaps_row = rep(c(gaps,2*gaps,3*gaps,4*gaps),each=2),annotation_colors=list(cluster=lineage.col.set[c("Epiblast","Endoderm","TE","Mesoderm","Amnion")]))


#' pheatmap of Blastoids (all cells)
temp.sel.cell <- temp.M.10X %>% filter(pj %in% c("IBD2","EBD2") )%>% pull(cell)
temp.exp <- lognormExp.mBN[de.gene.sel$Gene ,temp.sel.cell] 
temp.sel.exp <- t(apply(temp.exp,1,scale))
colnames(temp.sel.exp) <- colnames(temp.exp)
rownames(temp.sel.exp) <- rownames(temp.exp)
temp.sel.exp[temp.sel.exp>zs.limit] <- zs.limit
temp.sel.exp[temp.sel.exp<  (-1*zs.limit)] <- -1*zs.limit
gaps <- 10

temp.anno <-temp.M %>% filter(cell %in% temp.sel.cell) %>% select(cell,new_EML,pj)  %>% dplyr::rename(cluster=new_EML)%>% arrange(cluster) %>% mutate(od=factor(cluster,c("ELC","HLC","TLC"),ordered = T)) %>% arrange(od) %>% select(-od)  %>% tibble::column_to_rownames("cell") 
Plot.Exp$phIB <- pheatmap(temp.sel.exp[,rownames(temp.anno)[temp.anno$pj=="IBD2"]],cluster_rows=F,,cluster_cols=F,scale="none",annotation_col=temp.anno[,c("cluster"),drop=F],show_colnames=F,show_rownames=F,color=heat.col, fontsize_row=7,border_color="NA",main="\nIblastoids",gaps_row = rep(c(gaps,2*gaps,3*gaps,4*gaps),each=2),annotation_colors=list(cluster=lineage.col.set[c("ELC","HLC","TLC")]),legend=F,annotation_legend =F,annotation_names_col =F)%>%ggplotify::as.ggplot()

Plot.Exp$phEB <- pheatmap(temp.sel.exp[,rownames(temp.anno)[temp.anno$pj=="EBD2"]],cluster_rows=F,,cluster_cols=F,scale="none",annotation_col=temp.anno[,c("cluster"),drop=F],show_colnames=F,show_rownames=F,color=heat.col, fontsize_row=7,border_color="NA",main="\nStem-blastoids",gaps_row = rep(c(gaps,2*gaps,3*gaps,4*gaps),each=2),annotation_colors=list(cluster=lineage.col.set[c("ELC","HLC","TLC")]),legend=F,annotation_legend =F,annotation_names_col =F)%>%ggplotify::as.ggplot()

pheatmap(temp.sel.exp[,rownames(temp.anno)[temp.anno$pj=="EBD2"]],cluster_rows=F,,cluster_cols=F,scale="none",annotation_col=temp.anno[,c("cluster"),drop=F],show_colnames=F,show_rownames=F,color=heat.col, fontsize_row=7,border_color="NA",main="\nStem-blastoids",gaps_row = rep(c(gaps,2*gaps,3*gaps,4*gaps),each=2),annotation_colors=list(cluster=lineage.col.set[c("ELC","HLC","TLC")]),filename = "temp.pdf")

#' pheatmap using downsampled dataset
temp.sel.cell <- temp.M.10X %>% filter(pj %in% c("IBD2","EBD2") ) %>% filter(cell %in% meta.filter.ds$cell)%>% pull(cell)
temp.exp <- lognormExp.mBN[de.gene.sel$Gene ,temp.sel.cell] 
temp.sel.exp <- t(apply(temp.exp,1,scale))
colnames(temp.sel.exp) <- colnames(temp.exp)
rownames(temp.sel.exp) <- rownames(temp.exp)
temp.sel.exp[temp.sel.exp>zs.limit] <- zs.limit
temp.sel.exp[temp.sel.exp<  (-1*zs.limit)] <- -1*zs.limit
gaps <- 10

temp.anno <-temp.M %>% filter(cell %in% temp.sel.cell) %>% select(cell,new_EML,pj)  %>% dplyr::rename(cluster=new_EML)%>% arrange(cluster) %>% mutate(od=factor(cluster,c("ELC","HLC","TLC"),ordered = T)) %>% arrange(od) %>% select(-od)  %>% tibble::column_to_rownames("cell") 


Plot.Exp$phIB.ds <- pheatmap(temp.sel.exp[,rownames(temp.anno)[temp.anno$pj=="IBD2"]],cluster_rows=F,,cluster_cols=F,scale="none",annotation_col=temp.anno[,c("cluster"),drop=F],show_colnames=F,show_rownames=F,color=heat.col, fontsize_row=7,border_color="NA",main="\nIblastoids",gaps_row = rep(c(gaps,2*gaps,3*gaps,4*gaps),each=2),annotation_colors=list(cluster=lineage.col.set[c("ELC","HLC","TLC")]),legend=F,annotation_legend =F,annotation_names_col =F)%>%ggplotify::as.ggplot()

Plot.Exp$phEB.ds <- pheatmap(temp.sel.exp[,rownames(temp.anno)[temp.anno$pj=="EBD2"]],cluster_rows=F,,cluster_cols=F,scale="none",annotation_col=temp.anno[,c("cluster"),drop=F],show_colnames=F,show_rownames=F,color=heat.col, fontsize_row=7,border_color="NA",main="\nStem-blastoids",gaps_row = rep(c(gaps,2*gaps,3*gaps,4*gaps),each=2),annotation_colors=list(cluster=lineage.col.set[c("ELC","HLC","TLC")]),legend=F,annotation_legend =F,annotation_names_col =F)%>%ggplotify::as.ggplot()

pdf("D2.dot.plot.exp.pdf",7,5)
(
  (Plot.Exp$dot_legend/Plot.Exp$dot+plot_layout(heights = c(1,7)))|Plot.Exp$phEB.ds|Plot.Exp$phIB.ds
)
dev.off()

#' Gene module score analysis

temp.sel.cell <- temp.M.10X %>% filter(pj %in% c("IBD2","EBD2") )%>% pull(cell)
Lineage.markers <- temp.DM$sig %>% group_by(set)%>% filter(power >0.6) %>% top_n(10,power) %>% split(.$set)

temp.exp <- lognormExp.mBN[unlist(lapply(Lineage.markers,function(x){return(x$gene)})),temp.sel.cell] %>% as.data.frame() %>% tibble::rownames_to_column("gene") %>% tbl_df() %>% gather(cell,logExp,-gene) %>%inner_join(temp.M %>% select(cell,pj,new_EML),by="cell") %>% inner_join(temp.DM$sig  %>% select(gene,set),by="gene") 

temp.exp.10X <-temp.exp %>% group_by(gene) %>% summarise(logExp.mean=mean(logExp),logExp.sd=sd(logExp))

temp.exp.10X.mk.score <- temp.exp %>% inner_join(temp.exp.10X,by="gene") %>% mutate(sv=(logExp-logExp.mean)/logExp.sd) %>% mutate(sv=ifelse(sv> zs.limit,zs.limit,sv)) %>% mutate(sv=ifelse(sv< (-1*zs.limit),-1*zs.limit,sv))  %>% group_by(cell,pj,new_EML,set) %>% summarise(sv=mean(sv,na.rm = T)) 

temp.plot <- list()
for (p in c("IBD2","EBD2")) {
  for (n in unique(temp.exp.10X.mk.score$new_EML)) {
    temp.plot[[paste(n,p)]] <- temp.exp.10X.mk.score %>% filter(new_EML==n & pj==p) %>% mutate(x=set) %>% ggplot(mapping=aes(x=x,y=sv))+geom_violin(mapping=aes(fill=set),color="black",scale = "width")+geom_boxplot(color="grey50",fill="black",width=.2, outlier.shape = NA,coef = 0)+ theme_classic()  +NoLegend()+ggtitle(paste("10X",p))+theme(axis.text.x=element_text(angle = 40,hjust = 1),plot.title = element_text(hjust=0.5,face="bold"),panel.background = element_rect(fill = NA,colour="black"),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),panel.grid.major.x = element_line(colour = NA),strip.text = element_text(face = "bold"),strip.background = element_blank())+facet_wrap(. ~ new_EML,ncol = 3,scales="free")+ylab("Gene module score")+xlab("") +ylim(-1.5,2.25)+scale_fill_manual(values=lineage.col.set)
  }
}
print(cowplot::plot_grid(plotlist=temp.plot))





#c("RARRES2","ISL1","LRRN1","SERPING1","OCIAD2","GABRP")
#c("DPPA3","GTSF1","TGFBR3","AC022784.1","AP000547.3","REEP1")

#' violin plot for the topDE (TE vs AMN)
#temp.sel.gene <- c("ISL1","RARRES2","GABRP","PTGES","ADAM15","REEP1") # from the dot plot 
temp.sel.gene <- c("ISL1","GABRP","IGFBP2","GATA2","GTSF1","TGFBR3") # from the dot plot 

temp.sel.gene.big <-  c(temp.DM$sig %>% group_by(set)%>% filter(power >0.6)  %>% filter(set=="TE") %>% pull(gene)%>% head(100))

temp.M <- data.ob.umap %>% filter(EML %in% c("TE","IB_TE","EB_TLC","STB","EVT","CTB","Amnion","Tsw-AMLC","AMLC"))
temp.M.smt2 <- temp.M %>% filter(seqType=="smt2")
temp.M.10X <- temp.M %>% filter(seqType=="10X")


temp.sel.exp <- lognormExp.mBN[temp.sel.gene ,c(temp.M.10X$cell,temp.M.smt2$cell)]%>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,logExp,-Gene) %>% inner_join(data.ob.umap %>% select(cell,seqType,EML,cluster_EML) %>% mutate( bigEML= ifelse(EML %in% c("TE","STB","EVT","CTB"),'TE',EML)),by="cell") %>% mutate(od=factor(bigEML,c("Amnion","TE","AMLC","Tsw-AMLC","EB_TLC","IB_TE"),label=c("Amnion","TE","AMLC","Tsw-AMLC","TLC(SB)","TLC(IB)"),ordered = T)) %>% mutate(EML=ifelse(EML=="EB_TLC","TLC(SB)",EML))%>% mutate(EML=ifelse(EML=="IB_TE","TLC(IB)",EML))%>% mutate(seqType=factor(seqType, levels=c('smt2','10X'))) %>% arrange(od)  %>% mutate(EML=ifelse(EML%in% c("TE","STB","EVT","CTB"),'TE',EML))


Plot.exp.viol <- list()
Plot.exp.viol.tiff <- list()
for (g in temp.sel.gene ) {
  Plot.exp.viol[[g]] <- temp.sel.exp%>% filter(Gene==g) %>% ggplot(mapping=aes(x=od,y=logExp)) + geom_violin(mapping=aes(fill=EML),color="black",scale = "width")+geom_jitter(size=0.01)+ theme_classic() + theme(axis.text.x=element_text(angle = 40,hjust = 1)) +NoLegend()+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="bold"))+ylab("Expression level")+xlab("")+facet_grid(~seqType,scale="free_x")+theme(strip.background = element_blank(),strip.text.x = element_blank())+scale_fill_manual(values=vln.col.set)+theme(axis.text.x = element_text(size = 6))+ylab("")
  Plot.exp.viol.tiff[[g]] <- temp.sel.exp%>% filter(Gene==g) %>% ggplot(mapping=aes(x=od,y=logExp)) + geom_violin(mapping=aes(fill=EML),color="black",scale = "width")+ theme_classic() + theme(axis.text.x=element_text(angle = 40,hjust = 1)) +NoLegend()+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="bold"))+ylab("Expression level")+xlab("")+facet_grid(~seqType,scale="free_x")+theme(strip.background = element_blank(),strip.text.x = element_blank())+scale_fill_manual(values=vln.col.set)+theme(axis.text.x = element_text(size = 6))+ylab("")
}
pdf("D2.AmnVSTE.Vin.pdf",13.5,4.5)
print(
  plot_grid(plotlist=Plot.exp.viol,nrow=1,ncol=6)
)
print(
  plot_grid(plotlist=Plot.exp.viol.tiff,nrow=1,ncol=6)
)
dev.off()


#' check the TE marker gene expresion in 3D_ICM cells
temp.sel.M <- meta.filter %>% filter(pj=="D3post")
temp.sel.gene <- c("POU5F1","NANOG","DPPA3","BMP2","GATA2","GATA3")
temp.sel.exp <- lognormExp.mBN[temp.sel.gene ,temp.sel.M$cell]%>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,logExp,-Gene) %>% inner_join(data.ob.umap %>% select(cell,seqType,EML,EML),by="cell")
temp.plot <- list()
for (g in temp.sel.gene ) {
  temp.plot[[g]] <- temp.sel.exp%>% filter(Gene==g) %>% ggplot(mapping=aes(x=EML,y=logExp)) + geom_violin(mapping=aes(fill=EML),color="black",scale = "width")+geom_jitter(size=0.05)+ theme_classic() + theme(axis.text.x=element_text(angle = 40,hjust = 1)) +NoLegend()+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="bold"))+ylab("Expression level")+xlab("")+ylab("")
}
plot_grid(plotlist=temp.plot)



#' create the extra dotplot for the new annotation

data.cluster.text.pos <- data.ob.umap %>%group_by(seurat_clusters) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) %>% ungroup()
nCluster <- c(0:12)
nCluster.fill=data.frame(seurat_clusters=paste0("C",nCluster)) %>%tbl_df()%>% mutate(seurat_clusters=factor(seurat_clusters,paste0("C",nCluster),ordered = T))


Plot.umap$cluster <- ggplot()+ geom_point(data.ob.umap ,mapping=aes(x=UMAP_1,y=UMAP_2,color=seurat_clusters),size=0.5)+geom_text(data.cluster.text.pos  ,mapping=aes(x=UMAP_1,y=UMAP_2,label=seurat_clusters),size=4,color="Black")+theme_classic()+NoAxes()+NoLegend()+scale_color_manual(values=cluster.col.set)

Plot.umap$cluster.tiff <- ggplot()+ geom_point(data.ob.umap ,mapping=aes(x=UMAP_1,y=UMAP_2,color=seurat_clusters),size=0.5)+theme_classic()+NoAxes()+NoLegend()+scale_color_manual(values=cluster.col.set)

temp <- (data.ob.umap %>% filter(pj %in% c("CS7","D3post","SPH2016")) %>% group_by(seurat_clusters) %>% summarise(nCell=n_distinct(cell))%>% mutate(Type="Embryonic"))%>% bind_rows(data.ob.umap %>% filter(pj %in% c("CS7","D3post","SPH2016") & new_EML=="Epiblast")  %>% group_by(seurat_clusters) %>% summarise(nCell=n_distinct(cell))%>% mutate(Type="Epiblast")) %>% bind_rows(data.ob.umap %>% filter(pj %in% c("CS7","D3post","SPH2016") & new_EML=="Endoderm")  %>% group_by(seurat_clusters) %>% summarise(nCell=n_distinct(cell))%>% mutate(Type="Endoderm")) %>% bind_rows(data.ob.umap %>% filter(pj %in% c("CS7","D3post","SPH2016") & new_EML=="TE")  %>% group_by(seurat_clusters) %>% summarise(nCell=n_distinct(cell))%>% mutate(Type="TE")) %>% bind_rows(data.ob.umap %>% filter(pj %in% c("CS7","D3post","SPH2016") & new_EML=="Mesoderm")  %>% group_by(seurat_clusters) %>% summarise(nCell=n_distinct(cell))%>% mutate(Type="Mesoderm"))%>% bind_rows(data.ob.umap %>% filter(pj %in% c("CS7","D3post","SPH2016") & new_EML=="Amnion")  %>% group_by(seurat_clusters) %>% summarise(nCell=n_distinct(cell))%>% mutate(Type="Amnion")) %>% bind_rows(data.ob.umap %>% filter(pj %in% c("JPF2019") & new_EML=="AMLC")  %>% group_by(seurat_clusters) %>% summarise(nCell=n_distinct(cell))%>% mutate(Type="AMLC"))%>% bind_rows(data.ob.umap %>% filter(pj %in% c("JPF2019") & new_EML=="MeLC") %>% group_by(seurat_clusters) %>% summarise(nCell=n_distinct(cell))%>% mutate(Type="MeLC"))%>% bind_rows(data.ob.umap %>% filter(pj %in% c("EBD2") )  %>% group_by(seurat_clusters) %>% summarise(nCell=n_distinct(cell))%>% mutate(Type="Stem-blastoids"))%>% bind_rows(data.ob.umap %>% filter(pj %in% c("IBD2") )  %>% group_by(seurat_clusters) %>% summarise(nCell=n_distinct(cell))%>% mutate(Type="Iblastoids")) %>% spread(seurat_clusters,nCell) %>% replace(.,is.na(.),0) %>% gather(seurat_clusters,nCell,-Type) %>%  group_by(Type)   %>% mutate(Perc=nCell/sum(nCell))%>% mutate(seurat_clusters=factor(seurat_clusters,paste0("C",nCluster),ordered = T)) %>% filter(!Type %in% c("AMLC","MeLC")) %>% mutate(wrap_od=factor(Type,c("Embryonic","Epiblast","Endoderm","TE","Mesoderm","Amnion","Stem-blastoids","Iblastoids"),ordered = T))


Plot.cluster.prop <-  list()
Plot.cluster.prop$em <-  temp  %>% ggplot()+geom_bar(mapping=aes(x=seurat_clusters,y=Perc*100,fill=Type),width=0.5,stat = "identity")+theme_classic()+ylab("% Proportion of cells")+xlab("")+ ggtitle("")+theme(plot.title = element_text(hjust=0.5))+coord_flip() +facet_wrap(~wrap_od,nrow=1)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(strip.background=element_blank(),strip.text.y=element_blank())+scale_fill_manual(values=lineage.col.set)+NoLegend()+theme(strip.placement = NULL)#+ theme(legend.position="top",legend.title = "") +scale_y_discrete(position = "right")


temp.M <-  data.ob.umap %>% filter(pj %in% c("IBD2","EBD2"))
de.gene.sel <- temp.DM$sig %>% group_by(set)%>% filter(power >0.6) %>% top_n(10,power) %>%  dplyr::rename(Gene=gene) %>% mutate(od=factor(set,c("Epiblast","Endoderm","TE","Mesoderm","Amnion"),ordered = T)) %>% arrange(od,desc(power)) %>% select(-od)


temp <- lognormExp.mBN[de.gene.sel$Gene ,temp.M$cell]  %>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,logExp,-Gene) %>% inner_join(temp.M %>% select(cell,cluster_EML,pj) %>% dplyr::rename(cell.cluster=cluster_EML),by="cell") %>% inner_join(de.gene.sel%>% select(Gene,set) %>% dplyr::rename(Gene.cluster=set) ,by="Gene") 


temp.input <- temp %>% group_by(cell.cluster,Gene.cluster,pj,Gene) %>% summarise(meanlogExp=mean(logExp),nCell=n_distinct(cell)) %>% left_join(temp %>% group_by(cell.cluster,Gene.cluster,Gene,pj) %>% filter(logExp >0)%>% summarise(nExpCell=n_distinct(cell)),by=c("cell.cluster","Gene.cluster","Gene","pj"))  %>% mutate(nExpCell=ifelse(is.na(nExpCell),0,nExpCell))%>% mutate(nExp.ct=nExpCell/nCell) %>% group_by(Gene,pj) %>% mutate(sv=(meanlogExp-mean(meanlogExp))/sd(meanlogExp)) %>% mutate(nExp.ct= as.vector(cut(nExp.ct*100,c(-1,5,25,50,75,100),label=c(0,25,50,75,100)))) %>% mutate(nExp.ct=as.numeric(nExp.ct))


Plot.Exp$dot.EBD2 <- temp.input %>%filter(pj=="EBD2")  %>% mutate(Gene=factor(Gene,rev(de.gene.sel$Gene),ordered = T)) %>% mutate(cell.cluster=factor(cell.cluster,c("ELC","HLC","TLC","MeLC","AMLC","Undef"),ordered = T))%>% mutate(Gene.cluster=factor(Gene.cluster,c("Epiblast","Endoderm","TE","Mesoderm","Amnion"),ordered = T)) %>%  ggplot(aes(x=cell.cluster, y = Gene, color =sv, size = nExp.ct)) + geom_point() +facet_grid(Gene.cluster~.,,scale="free_y") + scale_color_viridis_c(name = 'Scaled Expression ') + cowplot::theme_cowplot() + theme(axis.line  = element_blank()) +theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +ylab('') +theme(axis.ticks = element_blank()) +scale_size(name="% cells expressing",breaks = c(0,25,50,75,100),range = c(0.1,3))+xlab("")+theme(axis.text.y = element_text(size = 6))+theme(axis.text.x = element_text(size = 7))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(strip.background=element_blank(),strip.text.y=element_blank())+NoLegend()


Plot.Exp$dot.IBD2 <- temp.input %>%filter(pj=="IBD2")  %>% mutate(Gene=factor(Gene,rev(de.gene.sel$Gene),ordered = T)) %>% mutate(cell.cluster=factor(cell.cluster,c("ELC","HLC","TLC","MeLC","AMLC","Undef"),ordered = T))%>% mutate(Gene.cluster=factor(Gene.cluster,c("Epiblast","Endoderm","TE","Mesoderm","Amnion"),ordered = T)) %>%  ggplot(aes(x=cell.cluster, y = Gene, color =sv, size = nExp.ct)) + geom_point() +facet_grid(Gene.cluster~.,,scale="free_y") + scale_color_viridis_c(name = 'Scaled Expression ') + cowplot::theme_cowplot() + theme(axis.line  = element_blank()) +theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +ylab('') +theme(axis.ticks = element_blank()) +scale_size(name="% cells expressing",breaks = c(0,25,50,75,100),range = c(0.1,3))+xlab("")+theme(axis.text.y = element_text(size = 6))+theme(axis.text.x = element_text(size = 7))+ theme(legend.position="right",legend.text = element_text(size=6),legend.title = element_text(size=6))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(strip.background=element_blank(),strip.text.y=element_blank())



Plot.cluster.num <- list()
Plot.cluster.num$EBD2 <- data.ob.umap %>% filter(pj %in% c("EBD2"))  %>% group_by(cluster_EML) %>% summarise(nCell=n_distinct(cell))  %>% mutate(Perc=nCell/sum(nCell)) %>% mutate(od=factor(cluster_EML,c("ELC","HLC","TLC","MeLC","AMLC","Undef"),ordered = T)) %>% arrange(od) %>% mutate(od=paste0(od,"\n(",nCell,")")) %>% mutate(od=factor(od,od,,ordered = T))  %>% ggplot()+geom_bar(mapping=aes(x=od,y=Perc*100,fill=cluster_EML),width=0.5,stat = "identity")+theme_classic()+ylab("% Proportion of cells")+xlab("")+ ggtitle("Stem-blastoids")+theme(plot.title = element_text(hjust=0.5))+ theme(legend.position="top",legend.title = element_blank())+scale_fill_manual(values=lineage.col.set)+ylim(0,50)+NoLegend()

Plot.cluster.num$IBD2 <- data.ob.umap %>% filter(pj %in% c("IBD2"))  %>% group_by(cluster_EML) %>% summarise(nCell=n_distinct(cell))  %>% mutate(Perc=nCell/sum(nCell)) %>% mutate(od=factor(cluster_EML,c("ELC","HLC","TLC","MeLC","AMLC","Undef"),ordered = T)) %>% arrange(od) %>% mutate(od=paste0(od,"\n(",nCell,")")) %>% mutate(od=factor(od,od,,ordered = T))  %>% ggplot()+geom_bar(mapping=aes(x=od,y=Perc*100,fill=cluster_EML),width=0.5,stat = "identity")+theme_classic()+ylab("% Proportion of cells")+xlab("")+ ggtitle("Iblastoids")+theme(plot.title = element_text(hjust=0.5))+ theme(legend.position="right",legend.title = element_blank())+scale_fill_manual(values=lineage.col.set)+ylim(0,50)


pdf("D2.figure.raw.pdf",9,14)
(
  plot_spacer()+Plot.umap$em+Plot.umap$em_legend+plot_spacer()+plot_layout(widths = c(0.5,2,1,0.5))
)/(
  (Plot.umap$JPF|Plot.umap$EBD2|Plot.umap$IBD2)
  
)/(
  (Plot.Exp$dot_legend/Plot.Exp$dot+plot_layout(heights = c(1,7)))|Plot.Exp$phIB.ds|Plot.Exp$phEB.ds
)/ (
  Plot.exp.viol[[1]]|Plot.exp.viol[[2]]|Plot.exp.viol[[3]]|
    Plot.exp.viol[[4]]|Plot.exp.viol[[5]]|Plot.exp.viol[[6]]
)+plot_layout(heights = c(2,1.5,4,0.8))
dev.off()


pdf("D2.figure.sup.raw.pdf",9,14)
wrap_plots(
  ( 
    plot_spacer()/Plot.umap$cluster/ plot_spacer()+plot_layout(heights = c(0.2,2,0.2))
  ) ,Plot.cluster.prop$em ,design="ABB"
)/(
  Plot.cluster.num$EBD2+Plot.cluster.num$IBD2+ Plot.Exp$dot.EBD2+ Plot.Exp$dot.IBD2 +plot_layout(design="AB\nCD\nCD")
)+plot_layout(heights = c(1,3))
wrap_plots(
  ( 
    plot_spacer()/ plot_spacer()/ plot_spacer()+plot_layout(heights = c(0.2,2,0.2))
  ) ,Plot.cluster.prop$em ,design="ABB"
)/(
  Plot.cluster.num$EBD2+Plot.cluster.num$IBD2+ Plot.Exp$dot.EBD2+ Plot.Exp$dot.IBD2 +plot_layout(design="AB\nCD\nCD")
)+plot_layout(heights = c(1,3))
dev.off()


#' cell contribution to blastoids

#' cells fell into C3, C8,  C5 as the blastoids cells
pie.plot <- list()

temp <- data.ob.umap %>% filter(devTime %in% c("E6","E7")) %>% filter(EML!="3D_ICM") %>% mutate(seurat_clusters=ifelse(seurat_clusters %in% c("C3","C5","C8","C10"),seurat_clusters,"Others")) %>% group_by(seurat_clusters,cluster_EML)%>% summarise(nCell=n_distinct(cell)) %>% mutate(seurat_clusters=factor(seurat_clusters,c("C3","C5","C8","C10","Others"),ordered=T)) %>% mutate(cluster_EML=factor(cluster_EML,c("Epiblast","Endoderm","TE"))) %>% arrange(cluster_EML,seurat_clusters)  %>% ungroup() %>% mutate(nCell.prop=nCell/sum(nCell)) %>% mutate(ymax=cumsum(nCell.prop)) %>% mutate(xmax=ifelse(cluster_EML=="Epiblast",2,NA)) %>% mutate(xmax=ifelse(cluster_EML=="Endoderm",3,xmax))%>% mutate(xmax=ifelse(cluster_EML=="TE",4,xmax)) %>% mutate(xmin=ifelse(cluster_EML=="Epiblast",1,NA)) %>% mutate(xmin=ifelse(cluster_EML=="Endoderm",2,xmin))%>% mutate(xmin=ifelse(cluster_EML=="TE",3,xmin)) 
temp$ymin=c(0, head(temp$ymax, n=-1)) 

temp.EML <- data.ob.umap %>% filter(devTime %in% c("E6","E7")) %>% filter(EML!="3D_ICM")  %>% group_by(cluster_EML)%>% summarise(nCell=n_distinct(cell))  %>% mutate(cluster_EML=factor(cluster_EML,c("Epiblast","Endoderm","TE"))) %>% arrange(cluster_EML)  %>% ungroup() %>% mutate(nCell.prop=nCell/sum(nCell))
temp.EML$ymin=c(0, head(temp.EML$ymax, n=-1))

pdf("E6_E7.prop.pdf",4.5,3)
print(
  ggplot() +geom_rect(temp,mapping= aes(ymax=ymax, ymin=ymin, xmax=xmax, xmin=xmin, fill=seurat_clusters),color="black") +scale_fill_manual(values=cluster.col.set)+ coord_polar(theta="y") + xlim(c(1, 4)) +theme_void()+ggtitle("E6/E7 cells")+FunTitle()
)
dev.off()

#data.ob.umap %>% filter(devTime %in% c("E6","E7")) %>% mutate(Type=ifelse(seurat_clusters %in% c("C3","C5","C8","C10"),"Blastoids clusters","Others")) %>% group_by(Type)%>% summarise(nCell=n_distinct(cell))  %>% mutate(nCell.perc=nCell/sum(nCell))%>%  ggplot( mapping=aes(x="", y=nCell, fill=Type))+ geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+ scale_fill_manual(values=c("Blastoids clusters"="royalblue3", "Others"="grey66"))+blank_theme +theme(axis.text.x=element_blank()) +geom_text(aes(y = nCell/2+c(0, cumsum(nCell)[-length(nCell)]),label = scales::percent(nCell.perc,accuracy=0.01)), size=5)+theme(legend.title = element_blank())+ggtitle("Blastoids(E6-E7)")+theme(plot.title = element_text(hjust=0.5))


pie.plot$EBD2 <- data.ob.umap %>% filter(pj %in% "EBD2") %>% mutate(Type=ifelse(seurat_clusters %in% c("C3","C5","C8","C10"),"Blastoids clusters","Others")) %>% group_by(Type)%>% summarise(nCell=n_distinct(cell))  %>% mutate(nCell.perc=nCell/sum(nCell)) %>% arrange((nCell.perc))%>%  ggplot( mapping=aes(x="", y=nCell, fill=Type))+ geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+ scale_fill_manual(values=c("Blastoids clusters"="royalblue3", "Others"="grey66"))+blank_theme +theme(axis.text.x=element_blank()) +geom_text(aes(y = nCell/2+c(0, cumsum(nCell)[-length(nCell)]),label = paste0(scales::percent(nCell.perc,accuracy=0.01),"\n(",nCell,")")), size=5)+theme(legend.title = element_blank())+ggtitle("Stem-blastoids")+theme(plot.title = element_text(hjust=0.5))

pie.plot$EBD2.cluster <- data.ob.umap %>% filter(pj %in% "EBD2") %>% mutate(Type=ifelse(seurat_clusters %in% c("C3","C5","C8","C10"),seurat_clusters,"Others")) %>% group_by(Type)%>% summarise(nCell=n_distinct(cell))  %>% mutate(nCell.perc=nCell/sum(nCell)) %>% arrange(desc(Type))%>%  ggplot( mapping=aes(x="", y=nCell, fill=Type))+ geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+ scale_fill_manual(values=cluster.col.set)+blank_theme +theme(axis.text.x=element_blank()) +geom_text(aes(y = nCell/2+c(0, cumsum(nCell)[-length(nCell)]),label = paste0(scales::percent(nCell.perc,accuracy=0.01),"\n(",nCell,")")), size=5)+theme(legend.title = element_blank())+ggtitle("Stem-blastoids")+theme(plot.title = element_text(hjust=0.5))



pie.plot$IBD2 <- data.ob.umap %>% filter(pj %in% "IBD2") %>% mutate(Type=ifelse(seurat_clusters %in% c("C3","C5","C8","C10"),"Blastoids clusters","Others")) %>% group_by(Type)%>% summarise(nCell=n_distinct(cell))  %>% mutate(nCell.perc=nCell/sum(nCell)) %>% arrange(desc(nCell.perc))%>%  ggplot( mapping=aes(x="", y=nCell, fill=Type))+ geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+ scale_fill_manual(values=c("Blastoids clusters"="royalblue3", "Others"="grey66"))+blank_theme +theme(axis.text.x=element_blank()) +geom_text(aes(y = nCell/2+c(0, cumsum(nCell)[-length(nCell)]),label = paste0(scales::percent(nCell.perc,accuracy=0.01),"\n(",nCell,")")), size=5)+theme(legend.title = element_blank())+ggtitle("Iblastoids")+theme(plot.title = element_text(hjust=0.5))

pie.plot$IBD2.cluster <- data.ob.umap %>% filter(pj %in% "IBD2") %>% mutate(Type=ifelse(seurat_clusters %in% c("C3","C5","C8","C10"),seurat_clusters,"Others")) %>% group_by(Type)%>% summarise(nCell=n_distinct(cell))  %>% mutate(nCell.perc=nCell/sum(nCell)) %>% arrange(desc(Type))%>%  ggplot( mapping=aes(x="", y=nCell, fill=Type))+ geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+ scale_fill_manual(values=cluster.col.set)+blank_theme +theme(axis.text.x=element_blank()) +geom_text(aes(y = nCell/2+c(0, cumsum(nCell)[-length(nCell)]),label = paste0(scales::percent(nCell.perc,accuracy=0.01),"\n(",nCell,")")), size=5)+theme(legend.title = element_blank())+ggtitle("Iblastoids")+theme(plot.title = element_text(hjust=0.5))+NoLegend()


pdf("D2.pie.plot.pdf",9,3)
pie.plot[[1]]+pie.plot[[2]]+pie.plot[[3]]+plot_layout(guides = 'collect')
dev.off()


pie.plot$EBD2.cluster+pie.plot$IBD2.cluster+plot_layout(guides = 'collect')



