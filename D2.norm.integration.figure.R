#' ---
#' title: "normalization, integration, defining markers and figure generating"
#' output: 
#'  html_document:
#'    code_folding: show
#' ---
 
rm(list=ls())
suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(scran))
suppressMessages(library(batchelor))
suppressMessages(library(Seurat))
suppressMessages(library(ggplot2))
suppressMessages(library(patchwork))
suppressMessages(library(pheatmap))

# working directory
DIR <- "~/My_project/JP_project"
setwd(DIR)


options(digits = 4)
options(future.globals.maxSize= 3001289600)

suppressMessages(library(foreach))
suppressMessages(library(doParallel))
numCores <- 10
registerDoParallel(numCores)

#' ## options 
opts <- list()
TD="D2_pub"

zs.limit <- 2.5




#' ### define I/O
opts$io <- list()
opts$io$counts.filter <- paste0("tmp_data/",TD,"/counts.filter.rds")
opts$io$meta.filter <- paste0("tmp_data/",TD,"/meta.filter.rds")
opts$io$meta.filter.ds <- paste0("tmp_data/",TD,"/meta.filter.ds.rds")
opts$io$lognormExp.mBN <- paste0("tmp_data/",TD,"/lognormExp.mBN.rds")
opts$io$data.ob <-  paste0("tmp_data/",TD,"/Pub.D2.mnn.rds")
opts$io$data.ob.umap <- paste0("tmp_data/",TD,"/Pub.D2.umap.cord.rds")
opts$io$markers <- paste0("tmp_data/",TD,"/DE.allL.Rdata")
opts$it <- list()
opts$it$nGene <- 2000
opts$it$opts$pc <- 25

#' #### loading counts and meta information
counts.filter <- readRDS(opts$io$counts.filter)
meta.filter<- readRDS(opts$io$meta.filter)%>% mutate(EML=ifelse(EML %in% c("Ectoderm"),"Amnion",EML))
meta.filter.ds <- readRDS(opts$io$meta.filter.ds)

#' #### loading figure setting
source("figures.setting.R")

#' #### loading local functions 
source("loca.quick.fun.R")
rename <- dplyr::rename
select<- dplyr::select
filter <- dplyr::filter

#' #### normalization 
if (file.exists(opts$io$lognormExp.mBN)) {
  lognormExp.mBN <- readRDS(opts$io$lognormExp.mBN)
}else{
  #' get the expressed genes ( expressed in at least 2 datatsets)
  expG.set <- list()
  for (b in unique(meta.filter$pj  %>% unique() %>% as.vector())) { 
    temp.cell <- meta.filter %>% filter(pj==b) %>% pull(cell)
    expG.set[[b]] <- rownames(counts.filter )[rowSums(counts.filter[,temp.cell] >=1) >=5]
  }
  sel.expG <-unlist(expG.set) %>% unique() %>% as.vector()
  
  # run scran normalization followed by multiBatchNorm 
  sce.ob <- list()
  
  for (b in unique(meta.filter$pj  %>% unique() %>% as.vector())) { 
    print(b)
    temp.M <- meta.filter %>% filter(pj==b) 
    temp.sce <-  SingleCellExperiment(list(counts=as.matrix(counts.filter[sel.expG,temp.M$cell])),colData=(temp.M %>% tibble::column_to_rownames("cell"))) %>% computeSumFactors()
    sce.ob[[b]] <- temp.sce
  }
  mBN.sce.ob <- multiBatchNorm(sce.ob$SPH2016,sce.ob$D3post,sce.ob$CS7,sce.ob$JPF2019,sce.ob$IBD2,sce.ob$EBD2)
  names(mBN.sce.ob) <- c("SPH2016","D3post","CS7","JPF2019","IBD2","EBD2")   
  lognormExp.mBN<- mBN.sce.ob %>% lapply(function(x) {logcounts(x) %>% as.data.frame()  %>% return()}) %>% do.call("bind_cols",.)
  saveRDS(lognormExp.mBN,opts$io$lognormExp.mBN)
  #' release memory
  rm(mBN.sce.ob)
}


#'  #### integrating all datasets
if (file.exists(opts$io$data.ob)) {
  data.ob <- readRDS(opts$io$data.ob)
  data.ob.umap <- readRDS(opts$io$data.ob.umap)
}else{
  
  data.merge <- CreateSeuratObject(counts.filter[rownames(lognormExp.mBN),meta.filter$cell], meta.data = (meta.filter %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) 
  data.merge@assays$RNA@data <- as.matrix(lognormExp.mBN[rownames(lognormExp.mBN),colnames(data.merge)])
  data.spt <- SplitObject(data.merge, split.by = "pj")%>% lapply(function(x){x=FindVariableFeatures(x,verbose=F,nfeatures=2000)})
  
  #' release memory
  rm(data.merge)
  
  nGene <- opts$it$nGene
  pc <- opts$it$opts$pc
  data.spt <- data.spt[c("JPF2019","IBD2","EBD2","D3post","CS7","SPH2016")]
  set.seed(123)
  data.ob <- SeuratWrappers::RunFastMNN(data.spt,verbose=F,features=nGene) %>% RunUMAP( reduction = "mnn", dims = 1:pc,verbose=F) %>% FindNeighbors( reduction = "mnn", dims = 1:pc) 
  data.temp <- data.ob %>% FindClusters(reso=0.4,verbose=F)
  rm(data.spt)
  
  
  data.ob.umap <- data.temp@meta.data %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(c(cell,SID:mt.perc)) %>% mutate(seurat_clusters=paste0("C",as.vector(Idents(data.temp))))  %>% inner_join(data.temp@reductions$umap@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df(),by="cell") %>% inner_join(data.temp@reductions$mnn@cell.embeddings %>% as.data.frame() %>% tibble::rownames_to_column("cell") %>% tbl_df() %>% select(cell:mnn_10),by="cell")

  #' update sph2016 annotation
  cell_lineage_data <- read.delim("~/My_project/ScRNA_analysis/data/Published/E-MTAB-3929.short.txt",stringsAsFactors = F,row.names = 1)%>% tibble::rownames_to_column("cell") %>% tbl_df()  %>% setNames(c("cell","Embryo","Days","Lineage1","Lineage2")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="primitive endoderm", "PrE")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="epiblast", "Epiblast")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="trophectoderm", "TE")) %>% mutate(Lineage1=replace(Lineage1, Lineage1=="not applicable", "Prelineage")) %>% mutate(new_lineage=Lineage1)%>% select(cell,new_lineage) ## using previous published annotation
  
  data.ob.umap <- data.ob.umap %>% left_join(cell_lineage_data,by="cell") %>% mutate(EML=ifelse(is.na(new_lineage),EML,new_lineage)) %>% select(-new_lineage)
  
  #' rename EML annotation
  data.ob.umap <- data.ob.umap %>% mutate(EML=ifelse(EML=="EM_NA","Prelineage",EML)) %>% mutate(EML=ifelse(EML=="EPI","Epiblast",EML))%>% mutate(EML=ifelse(EML=="PE","PrE",EML))  %>% mutate(EML=ifelse(pj=="D3post" & EML=="ICM","3D_ICM",EML))
  
  #' create big annotation
  data.ob.umap <- data.ob.umap %>% mutate(new_EML=EML) %>% mutate(new_EML=ifelse(new_EML%in% c("AdvMes","AxMes","EmMes","NasMes"),"Mesoderm",new_EML)) %>% mutate(new_EML=ifelse(new_EML%in% c("CTB","EVT","STB","TE"),"TE",new_EML)) %>% mutate(new_EML=ifelse(new_EML%in% c("MeLC2","MeLC1"),"MeLC",new_EML)) %>% mutate(new_EML=ifelse(new_EML%in% c("Tsw-AMLC","AMLC"),"AMLC",new_EML)) %>% mutate(new_EML=ifelse(new_EML%in% c("PrE"),"Endoderm",new_EML))  %>% mutate(new_EML=ifelse(new_EML=="EB_ELC","ELC",new_EML)) %>% mutate(new_EML=ifelse(new_EML=="EB_HLC","HLC",new_EML))%>% mutate(new_EML=ifelse(new_EML=="EB_TLC","TLC",new_EML)) %>% mutate(new_EML=ifelse(new_EML=="IB_Epi","ELC",new_EML))%>% mutate(new_EML=ifelse(new_EML=="IB_PE","HLC",new_EML))%>% mutate(new_EML=ifelse(new_EML=="IB_TE","TLC",new_EML))  %>% mutate(new_EML=ifelse(new_EML%in% c("EB_U10","EB_U5","EB_U6","EB_UE7","EB_UI13","IB_IM1","IB_IM2","IB_IM3","IB_uc","IB_NR"),"Undef",new_EML))
  
  #' create cluster information
  data.ob.umap <- data.ob.umap %>% filter(pj %in% c("EBD2","IBD2")) %>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C4","C8"),"ELC","Undef")) %>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C5"),"HLC",cluster_EML))%>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C3"),"TLC",cluster_EML))%>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C0","C2","C11"),"MeLC",cluster_EML)) %>% mutate(cluster_EML=ifelse(seurat_clusters %in% c("C1","C6"),"AMLC",cluster_EML)) %>% bind_rows(data.ob.umap %>% filter(!pj %in% c("EBD2","IBD2")) %>% mutate(cluster_EML=new_EML))
  
  saveRDS(data.ob,file=opts$io$data.ob)
  saveRDS(data.ob.umap,file=opts$io$data.ob.umap)
}

#' #### extract marker genes from embryonic cells based on published annotation
if (file.exists(opts$io$markers)) {
  load(opts$io$markers,verbose=T)
}else{
  temp.M <-  data.ob.umap %>% filter(new_EML %in% c("Epiblast","ELC","TE","TLC","Endoderm","HLC","Mesoderm","MeLC","Amnion","AMLC"))
  temp.M.smt2 <- temp.M %>% filter(seqType=="smt2")
  temp.M.10X <- temp.M %>% filter(seqType=="10X")
  
  expG.set <- list()
  expG.set$smt2 <- rownames(counts.filter )[rowSums(counts.filter[,temp.M.smt2$cell] >=1) >=5] # the most frequent expressed gene
  expG.set$X10 <- rownames(counts.filter )[rowSums(counts.filter[,temp.M.10X$cell] >=1) >=5]# the most frequent expressed gene
  
  temp.expG <- intersect(expG.set$smt2,expG.set$X10) %>% intersect(rownames(lognormExp.mBN))
  
  #' for smart-seq dataset only
  data.temp <- CreateSeuratObject(counts.filter[temp.expG ,temp.M.smt2$cell], meta.data = (temp.M.smt2 %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) 
  data.temp@assays$RNA@data <- as.matrix(lognormExp.mBN[rownames(data.temp),colnames(data.temp)])
  
  Idents(data.temp) <- as.factor(data.temp@meta.data$new_EML)
  
  #' using the roc test
  de.out <- FindAllMarkers(data.temp,test.use="roc",verbose=F)  %>% tbl_df() %>%rename(Gene=gene)
  temp.DM <- FunRF_FindAllMarkers_para(data.temp)
  save(de.out,temp.DM, file=opts$io$markers)
}

#' ## generating figures
# text position embryonic dataset
data.EM.umap.text.pos <- data.ob.umap  %>% filter(pj %in% c("D3post","CS7","SPH2016"))  %>%group_by(new_EML) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) %>% ungroup()  %>% filter(!new_EML %in% c("3D_ICM","PSA-EPI")) %>% mutate(UMAP_1=ifelse(new_EML=="PriS",UMAP_1-1,UMAP_1)) %>% mutate(new_EML=gsub("ExE_Mes","Extraembryonic\nMesoderm",new_EML))%>% mutate(new_EML=gsub("PriS","Primitive\nStreak",new_EML))

# text position of other datasets
data.LC.umap.text.pos <- data.ob.umap  %>% filter(pj %in% c("EBD2","IBD2","JPF2019"))  %>%group_by(new_EML,pj) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) %>% ungroup() %>% mutate(UMAP_2=ifelse(new_EML=="HLC",UMAP_2+0.5,UMAP_2)) # UMAP_2 of HLC adding 0.5 to avoid overlapping

#' #### umaps
Plot.umap <- list() 
Plot.umap$em <- ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("D3post","CS7","SPH2016")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.5)+geom_point( data.ob.umap %>% filter(pj %in% c("D3post","CS7","SPH2016"))  %>% mutate(od=factor(new_EML,c("Prelineage","TE","Mesoderm","Epiblast","Amnion","Endoderm","PriS","ExE_Mes","3D_ICM","PSA-EPI"),ordered = T)) %>% arrange(desc(od)) ,mapping=aes(x=UMAP_1,y=UMAP_2,fill=new_EML,shape=pj),size=2,color="grey")+ggtitle("Embryonic dataset")+geom_text(data.EM.umap.text.pos ,mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=4,color="Black",fontface="bold")+theme_classic()+theme(plot.title = element_text(hjust=0.5))+scale_fill_manual(values=lineage.col.set)+scale_shape_manual(values = pj.shape.set,labels = c('SPH2016'="Petropoulos et al., 2016",'D3post'="Xiang et al., 2019",'CS7'="Tyser et al., 2020"))+NoAxes()+NoLegend()

umap.xlim=ggplot_build(Plot.umap$em)$layout$panel_scales_x[[1]]$range$range
umap.ylim=ggplot_build(Plot.umap$em)$layout$panel_scales_y[[1]]$range$range

Plot.umap$em_legend <- ggplot()+ geom_point( data.ob.umap %>% filter(pj %in% c("D3post","CS7","SPH2016")) %>% filter(!new_EML %in% c("PSA-EPI","3D_ICM")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=new_EML,shape=pj),size=2)+theme_classic()+scale_color_manual("",values=lineage.col.set)+scale_shape_manual("",values = pj.shape.set.solid,labels = c('SPH2016'="Petropoulos et al., 2016",'D3post'="Xiang et al., 2020",'CS7'="Tyser et al., 2020")) + theme(legend.text = element_text( size = 6),legend.title = element_text( size = 1))
Plot.umap$em_legend <- Plot.umap$em_legend%>% ggpubr::get_legend() %>% ggpubr::as_ggplot()

Plot.umap$JPF <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("JPF2019")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("JPF2019")) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=new_EML),size=0.5)+ggtitle("Zheng et al., 2019")+geom_text(data.EM.umap.text.pos %>% filter(new_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=3.5,color="Black")+  geom_text(data.LC.umap.text.pos %>% filter(pj=="JPF2019"),mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=4,color="Black",fontface="bold")  +theme_classic()+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+scale_color_manual(values=lineage.col.set)+xlim(umap.xlim)+ylim(umap.ylim)


Plot.umap$EBD2 <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("EBD2")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("EBD2"))%>% mutate(od=factor(new_EML,c("ELC","HLC","TLC","Undef"),ordered = T)) %>% arrange(desc(od)) ,mapping=aes(x=UMAP_1,y=UMAP_2,col=new_EML),size=0.5)+ggtitle("Stem-blastoids")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+geom_text(data.EM.umap.text.pos %>% filter(new_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=3.5,color="Black")+geom_text(data.LC.umap.text.pos %>% filter(pj=="EBD2") %>% filter(new_EML %in% c("ELC","HLC","TLC")),mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=4,color="Black",fontface="bold")+xlim(umap.xlim)+ylim(umap.ylim)

Plot.umap$IBD2 <-  ggplot()+ geom_point(data.ob.umap %>% filter(!pj %in% c("IBD2")) ,mapping=aes(x=UMAP_1,y=UMAP_2),color="grey",shape=16,size=0.25)+geom_point( data.ob.umap %>% filter(pj%in% c("IBD2"))  %>% mutate(od=factor(new_EML,c("ELC","HLC","TLC","Undef"),ordered = T)) %>% arrange(desc(od)),mapping=aes(x=UMAP_1,y=UMAP_2,col=new_EML),size=0.5)+ggtitle("Iblastoids")+theme_classic()+scale_color_manual(values=lineage.col.set)+NoLegend()+NoAxes()+theme(plot.title = element_text(hjust=0.5))+geom_text(data.EM.umap.text.pos %>% filter(new_EML %in% c("Amnion","Endoderm","Epiblast","Mesoderm","TE")),mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=3.5,color="Black")+geom_text(data.LC.umap.text.pos %>% filter(pj=="IBD2") %>% filter(new_EML %in% c("ELC","HLC","TLC")),mapping=aes(x=UMAP_1,y=UMAP_2,label=new_EML),size=4,color="Black",fontface="bold")+xlim(umap.xlim)+ylim(umap.ylim)

#' #### dot plot and heatmap
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
# temp.sel.cell <- temp.M.smt2%>% pull(cell)
# temp.exp <- lognormExp.mBN[de.gene.sel$Gene ,temp.sel.cell] 
# temp.sel.exp <- t(apply(temp.exp,1,scale))
# colnames(temp.sel.exp) <- colnames(temp.exp)
# rownames(temp.sel.exp) <- rownames(temp.exp)
# temp.sel.exp[temp.sel.exp>zs.limit] <- zs.limit
# temp.sel.exp[temp.sel.exp<  (-1*zs.limit)] <- -1*zs.limit
# gaps <- 10
# 
# temp.anno <-temp.M %>% filter(cell %in% temp.sel.cell) %>% select(cell,new_EML,pj)  %>% dplyr::rename(cluster=new_EML)%>% arrange(cluster) %>% mutate(od=factor(cluster,c("Epiblast","Endoderm","TE","Mesoderm","Amnion"),ordered = T)) %>% arrange(od) %>% select(-od)  %>% tibble::column_to_rownames("cell") 
# pheatmap(temp.sel.exp[,rownames(temp.anno)],cluster_rows=F,,cluster_cols=F,scale="none",annotation_col=temp.anno[,c("cluster","pj"),drop=F],show_colnames=F,show_rownames=T,color=heat.col, fontsize_row=7,border_color="NA",main="\nEmbryonic cells",gaps_row = rep(c(gaps,2*gaps,3*gaps,4*gaps),each=2),annotation_colors=list(cluster=lineage.col.set[c("Epiblast","Endoderm","TE","Mesoderm","Amnion")]))


temp.sel.cell <-data.ob.umap %>% filter(new_EML %in% c("Epiblast","ELC","TE","TLC","Endoderm","HLC","Mesoderm","MeLC","Amnion","AMLC")) %>% filter(seqType=="10X") %>% filter(pj %in% c("IBD2","EBD2") ) %>% filter(cell %in% meta.filter.ds$cell)%>% pull(cell)
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


#' #### vlnplot
temp.sel.gene <- c("ISL1","GABRP","IGFBP2","GATA2","GTSF1","TGFBR3") # from the dot plot 

temp.M <- data.ob.umap %>% filter(EML %in% c("TE","IB_TE","EB_TLC","STB","EVT","CTB","Amnion","Tsw-AMLC","AMLC"))
temp.M.smt2 <- temp.M %>% filter(seqType=="smt2")
temp.M.10X <- temp.M %>% filter(seqType=="10X")


temp.sel.exp <- lognormExp.mBN[temp.sel.gene ,c(temp.M.10X$cell,temp.M.smt2$cell)]%>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,logExp,-Gene) %>% inner_join(data.ob.umap %>% select(cell,seqType,EML,cluster_EML) %>% mutate( bigEML= ifelse(EML %in% c("TE","STB","EVT","CTB"),'TE',EML)),by="cell") %>% mutate(od=factor(bigEML,c("Amnion","TE","AMLC","Tsw-AMLC","EB_TLC","IB_TE"),label=c("Amnion","TE","AMLC","Tsw-AMLC","TLC(SB)","TLC(IB)"),ordered = T)) %>% mutate(EML=ifelse(EML=="EB_TLC","TLC(SB)",EML))%>% mutate(EML=ifelse(EML=="IB_TE","TLC(IB)",EML))%>% mutate(seqType=factor(seqType, levels=c('smt2','10X'))) %>% arrange(od)  %>% mutate(EML=ifelse(EML%in% c("TE","STB","EVT","CTB"),'TE',EML))


Plot.exp.viol <- list()
for (g in temp.sel.gene ) {
  Plot.exp.viol[[g]] <- temp.sel.exp%>% filter(Gene==g) %>% ggplot(mapping=aes(x=od,y=logExp)) + geom_violin(mapping=aes(fill=EML),color="black",scale = "width")+geom_jitter(size=0.01)+ theme_classic() + theme(axis.text.x=element_text(angle = 40,hjust = 1)) +NoLegend()+ggtitle(g)+theme(plot.title = element_text(hjust=0.5,face="bold"))+ylab("Expression level")+xlab("")+facet_grid(~seqType,scale="free_x")+theme(strip.background = element_blank(),strip.text.x = element_blank())+scale_fill_manual(values=vln.col.set)+theme(axis.text.x = element_text(size = 6))+ylab("")
}

#' #### summary reannotation

data.cluster.text.pos <- data.ob.umap %>%group_by(seurat_clusters) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)) %>% ungroup()
nCluster <- c(0:12)
nCluster.fill=data.frame(seurat_clusters=paste0("C",nCluster)) %>%tbl_df()%>% mutate(seurat_clusters=factor(seurat_clusters,paste0("C",nCluster),ordered = T))


Plot.umap$cluster <- ggplot()+ geom_point(data.ob.umap ,mapping=aes(x=UMAP_1,y=UMAP_2,color=seurat_clusters),size=0.5)+geom_text(data.cluster.text.pos  ,mapping=aes(x=UMAP_1,y=UMAP_2,label=seurat_clusters),size=4,color="Black")+theme_classic()+NoAxes()+NoLegend()+scale_color_manual(values=cluster.col.set)

temp <- (data.ob.umap %>% filter(pj %in% c("CS7","D3post","SPH2016")) %>% group_by(seurat_clusters) %>% summarise(nCell=n_distinct(cell))%>% mutate(Type="Embryonic"))%>% bind_rows(data.ob.umap %>% filter(pj %in% c("CS7","D3post","SPH2016") & new_EML=="Epiblast")  %>% group_by(seurat_clusters) %>% summarise(nCell=n_distinct(cell))%>% mutate(Type="Epiblast")) %>% bind_rows(data.ob.umap %>% filter(pj %in% c("CS7","D3post","SPH2016") & new_EML=="Endoderm")  %>% group_by(seurat_clusters) %>% summarise(nCell=n_distinct(cell))%>% mutate(Type="Endoderm")) %>% bind_rows(data.ob.umap %>% filter(pj %in% c("CS7","D3post","SPH2016") & new_EML=="TE")  %>% group_by(seurat_clusters) %>% summarise(nCell=n_distinct(cell))%>% mutate(Type="TE")) %>% bind_rows(data.ob.umap %>% filter(pj %in% c("CS7","D3post","SPH2016") & new_EML=="Mesoderm")  %>% group_by(seurat_clusters) %>% summarise(nCell=n_distinct(cell))%>% mutate(Type="Mesoderm"))%>% bind_rows(data.ob.umap %>% filter(pj %in% c("CS7","D3post","SPH2016") & new_EML=="Amnion")  %>% group_by(seurat_clusters) %>% summarise(nCell=n_distinct(cell))%>% mutate(Type="Amnion")) %>% bind_rows(data.ob.umap %>% filter(pj %in% c("JPF2019") & new_EML=="AMLC")  %>% group_by(seurat_clusters) %>% summarise(nCell=n_distinct(cell))%>% mutate(Type="AMLC"))%>% bind_rows(data.ob.umap %>% filter(pj %in% c("JPF2019") & new_EML=="MeLC") %>% group_by(seurat_clusters) %>% summarise(nCell=n_distinct(cell))%>% mutate(Type="MeLC"))%>% bind_rows(data.ob.umap %>% filter(pj %in% c("EBD2") )  %>% group_by(seurat_clusters) %>% summarise(nCell=n_distinct(cell))%>% mutate(Type="Stem-blastoids"))%>% bind_rows(data.ob.umap %>% filter(pj %in% c("IBD2") )  %>% group_by(seurat_clusters) %>% summarise(nCell=n_distinct(cell))%>% mutate(Type="Iblastoids")) %>% spread(seurat_clusters,nCell) %>% replace(.,is.na(.),0) %>% gather(seurat_clusters,nCell,-Type) %>%  group_by(Type)   %>% mutate(Perc=nCell/sum(nCell))%>% mutate(seurat_clusters=factor(seurat_clusters,paste0("C",nCluster),ordered = T)) %>% filter(!Type %in% c("AMLC","MeLC")) %>% mutate(wrap_od=factor(Type,c("Embryonic","Epiblast","Endoderm","TE","Mesoderm","Amnion","Stem-blastoids","Iblastoids"),ordered = T))


Plot.cluster.prop <-  list()
Plot.cluster.prop$em <-  temp  %>% ggplot()+geom_bar(mapping=aes(x=seurat_clusters,y=Perc*100,fill=Type),width=0.5,stat = "identity")+theme_classic()+ylab("% Proportion of cells")+xlab("")+ ggtitle("")+theme(plot.title = element_text(hjust=0.5))+coord_flip() +facet_wrap(~wrap_od,nrow=1)+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(strip.background=element_blank(),strip.text.y=element_blank())+scale_fill_manual(values=lineage.col.set)+NoLegend()+theme(strip.placement = NULL)#+ theme(legend.position="top",legend.title = "") +scale_y_discrete(position = "right")


Plot.cluster.num <- list()
Plot.cluster.num$EBD2 <- data.ob.umap %>% filter(pj %in% c("EBD2"))  %>% group_by(cluster_EML) %>% summarise(nCell=n_distinct(cell))  %>% mutate(Perc=nCell/sum(nCell)) %>% mutate(od=factor(cluster_EML,c("ELC","HLC","TLC","MeLC","AMLC","Undef"),ordered = T)) %>% arrange(od) %>% mutate(od=paste0(od,"\n(",nCell,")")) %>% mutate(od=factor(od,od,,ordered = T))  %>% ggplot()+geom_bar(mapping=aes(x=od,y=Perc*100,fill=cluster_EML),width=0.5,stat = "identity")+theme_classic()+ylab("% Proportion of cells")+xlab("")+ ggtitle("Stem-blastoids")+theme(plot.title = element_text(hjust=0.5))+ theme(legend.position="top",legend.title = element_blank())+scale_fill_manual(values=lineage.col.set)+ylim(0,50)+NoLegend()

Plot.cluster.num$IBD2 <- data.ob.umap %>% filter(pj %in% c("IBD2"))  %>% group_by(cluster_EML) %>% summarise(nCell=n_distinct(cell))  %>% mutate(Perc=nCell/sum(nCell)) %>% mutate(od=factor(cluster_EML,c("ELC","HLC","TLC","MeLC","AMLC","Undef"),ordered = T)) %>% arrange(od) %>% mutate(od=paste0(od,"\n(",nCell,")")) %>% mutate(od=factor(od,od,,ordered = T))  %>% ggplot()+geom_bar(mapping=aes(x=od,y=Perc*100,fill=cluster_EML),width=0.5,stat = "identity")+theme_classic()+ylab("% Proportion of cells")+xlab("")+ ggtitle("Iblastoids")+theme(plot.title = element_text(hjust=0.5))+ theme(legend.position="right",legend.title = element_blank())+scale_fill_manual(values=lineage.col.set)+ylim(0,50)


#' #### dot plot for reannotation of two balstoids datasets
temp.M <-  data.ob.umap %>% filter(pj %in% c("IBD2","EBD2"))
de.gene.sel <- temp.DM$sig %>% group_by(set)%>% filter(power >0.6) %>% top_n(10,power) %>%  dplyr::rename(Gene=gene) %>% mutate(od=factor(set,c("Epiblast","Endoderm","TE","Mesoderm","Amnion"),ordered = T)) %>% arrange(od,desc(power)) %>% select(-od)
temp <- lognormExp.mBN[de.gene.sel$Gene ,temp.M$cell]  %>% tibble::rownames_to_column("Gene") %>% tbl_df() %>% gather(cell,logExp,-Gene) %>% inner_join(temp.M %>% select(cell,cluster_EML,pj) %>% dplyr::rename(cell.cluster=cluster_EML),by="cell") %>% inner_join(de.gene.sel%>% select(Gene,set) %>% dplyr::rename(Gene.cluster=set) ,by="Gene") 
temp.input <- temp %>% group_by(cell.cluster,Gene.cluster,pj,Gene) %>% summarise(meanlogExp=mean(logExp),nCell=n_distinct(cell)) %>% left_join(temp %>% group_by(cell.cluster,Gene.cluster,Gene,pj) %>% filter(logExp >0)%>% summarise(nExpCell=n_distinct(cell)),by=c("cell.cluster","Gene.cluster","Gene","pj"))  %>% mutate(nExpCell=ifelse(is.na(nExpCell),0,nExpCell))%>% mutate(nExp.ct=nExpCell/nCell) %>% group_by(Gene,pj) %>% mutate(sv=(meanlogExp-mean(meanlogExp))/sd(meanlogExp)) %>% mutate(nExp.ct= as.vector(cut(nExp.ct*100,c(-1,5,25,50,75,100),label=c(0,25,50,75,100)))) %>% mutate(nExp.ct=as.numeric(nExp.ct))


Plot.Exp$dot.EBD2 <- temp.input %>%filter(pj=="EBD2")  %>% mutate(Gene=factor(Gene,rev(de.gene.sel$Gene),ordered = T)) %>% mutate(cell.cluster=factor(cell.cluster,c("ELC","HLC","TLC","MeLC","AMLC","Undef"),ordered = T))%>% mutate(Gene.cluster=factor(Gene.cluster,c("Epiblast","Endoderm","TE","Mesoderm","Amnion"),ordered = T)) %>%  ggplot(aes(x=cell.cluster, y = Gene, color =sv, size = nExp.ct)) + geom_point() +facet_grid(Gene.cluster~.,,scale="free_y") + scale_color_viridis_c(name = 'Scaled Expression ') + cowplot::theme_cowplot() + theme(axis.line  = element_blank()) +theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +ylab('') +theme(axis.ticks = element_blank()) +scale_size(name="% cells expressing",breaks = c(0,25,50,75,100),range = c(0.1,3))+xlab("")+theme(axis.text.y = element_text(size = 6))+theme(axis.text.x = element_text(size = 7))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(strip.background=element_blank(),strip.text.y=element_blank())+NoLegend()


Plot.Exp$dot.IBD2 <- temp.input %>%filter(pj=="IBD2")  %>% mutate(Gene=factor(Gene,rev(de.gene.sel$Gene),ordered = T)) %>% mutate(cell.cluster=factor(cell.cluster,c("ELC","HLC","TLC","MeLC","AMLC","Undef"),ordered = T))%>% mutate(Gene.cluster=factor(Gene.cluster,c("Epiblast","Endoderm","TE","Mesoderm","Amnion"),ordered = T)) %>%  ggplot(aes(x=cell.cluster, y = Gene, color =sv, size = nExp.ct)) + geom_point() +facet_grid(Gene.cluster~.,,scale="free_y") + scale_color_viridis_c(name = 'Scaled Expression ') + cowplot::theme_cowplot() + theme(axis.line  = element_blank()) +theme(axis.text.x = element_text(angle = 60, vjust = 1, hjust=1)) +ylab('') +theme(axis.ticks = element_blank()) +scale_size(name="% cells expressing",breaks = c(0,25,50,75,100),range = c(0.1,3))+xlab("")+theme(axis.text.y = element_text(size = 6))+theme(axis.text.x = element_text(size = 7))+ theme(legend.position="right",legend.text = element_text(size=6),legend.title = element_text(size=6))+theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())+theme(strip.background=element_blank(),strip.text.y=element_blank())

# ## fig1
#+ fig.width=9, fig.height=14
(
  plot_spacer()+Plot.umap$em+Plot.umap$em_legend+plot_spacer()+plot_layout(widths = c(0.5,2,1,0.5))
)/(
  (Plot.umap$JPF|Plot.umap$EBD2|Plot.umap$IBD2)
  
)/(
  (Plot.Exp$dot_legend/Plot.Exp$dot+plot_layout(heights = c(1,7)))|Plot.Exp$phEB.ds|Plot.Exp$phIB.ds
)/ (
  Plot.exp.viol[[1]]|Plot.exp.viol[[2]]|Plot.exp.viol[[3]]|
    Plot.exp.viol[[4]]|Plot.exp.viol[[5]]|Plot.exp.viol[[6]]
)+plot_layout(heights = c(2,1.5,4,0.8))


#' ## sup figures
#+ fig.width=9, fig.height=14
wrap_plots(
  ( 
    plot_spacer()/Plot.umap$cluster/ plot_spacer()+plot_layout(heights = c(0.2,2,0.2))
  ) ,Plot.cluster.prop$em ,design="ABB"
)/(
  Plot.cluster.num$EBD2+Plot.cluster.num$IBD2+ Plot.Exp$dot.EBD2+ Plot.Exp$dot.IBD2 +plot_layout(design="AB\nCD\nCD")
)+plot_layout(heights = c(1,3))



#' ## proportion of clusters of blastcysts

# cells fell into C3, C8,  C5 as the blastoids cells
pie.plot <- list()

temp <- data.ob.umap %>% filter(devTime %in% c("E6","E7")) %>% filter(EML!="3D_ICM") %>% mutate(seurat_clusters=ifelse(seurat_clusters %in% c("C3","C5","C8","C10"),seurat_clusters,"Others")) %>% group_by(seurat_clusters,cluster_EML)%>% summarise(nCell=n_distinct(cell)) %>% mutate(seurat_clusters=factor(seurat_clusters,c("C3","C5","C8","C10","Others"),ordered=T)) %>% mutate(cluster_EML=factor(cluster_EML,c("Epiblast","Endoderm","TE"))) %>% arrange(cluster_EML,seurat_clusters)  %>% ungroup() %>% mutate(nCell.prop=nCell/sum(nCell)) %>% mutate(ymax=cumsum(nCell.prop)) %>% mutate(xmax=ifelse(cluster_EML=="Epiblast",2,NA)) %>% mutate(xmax=ifelse(cluster_EML=="Endoderm",3,xmax))%>% mutate(xmax=ifelse(cluster_EML=="TE",4,xmax)) %>% mutate(xmin=ifelse(cluster_EML=="Epiblast",1,NA)) %>% mutate(xmin=ifelse(cluster_EML=="Endoderm",2,xmin))%>% mutate(xmin=ifelse(cluster_EML=="TE",3,xmin)) 
temp$ymin=c(0, head(temp$ymax, n=-1)) 

temp.EML <- data.ob.umap %>% filter(devTime %in% c("E6","E7")) %>% filter(EML!="3D_ICM")  %>% group_by(cluster_EML)%>% summarise(nCell=n_distinct(cell))  %>% mutate(cluster_EML=factor(cluster_EML,c("Epiblast","Endoderm","TE"))) %>% arrange(cluster_EML)  %>% ungroup() %>% mutate(nCell.prop=nCell/sum(nCell))
temp.EML$ymin=c(0, head(temp.EML$ymax, n=-1))

#+ fig.width=4.5, fig.height=3
print(
  ggplot() +geom_rect(temp,mapping= aes(ymax=ymax, ymin=ymin, xmax=xmax, xmin=xmin, fill=seurat_clusters),color="black") +scale_fill_manual(values=cluster.col.set)+ coord_polar(theta="y") + xlim(c(1, 4)) +theme_void()+ggtitle("E6/E7 cells")+ theme(plot.title = element_text(hjust=0.5,size=15,face="bold"))
)

sessionInfo()

