FunTitle <- function() {
  p <- theme(plot.title = element_text(hjust=0.5,size=15,face="bold"))
  return(p)
}
FunHandleINT <- function(x) {
  x.int <- x %>% filter(ref_EML %in% c("EPI_Amnion","EPI.PrE.INT","PriS_Amnion","EPI_PriS"))
  
  x.unique <-  x %>% filter(!ref_EML %in% c("EPI_Amnion","EPI.PrE.INT","PriS_Amnion","EPI_PriS")) %>% filter(ref_EML %in% c("Amnion","PriS","Epi") & ref_pj=="CS7") %>% bind_rows(x %>% filter(!ref_EML %in% c("EPI_Amnion","EPI.PrE.INT","PriS_Amnion","EPI_PriS")) %>% filter(ref_EML %in% c("Epi","Endoderm") & ref_pj=="SPH2016"))
                                                                                                            
  
  x.basic <- x   %>% anti_join(x.int, by = c("query_cell","mean_pval","ref_pj","nc_hits")) %>% anti_join(x.unique,by = c("query_cell","mean_pval","ref_pj","nc_hits"))
  
  if (nrow(x.int) >0 & nrow(x.unique) >0) {
    x.int.type <- x.int %>% mutate(INT.ref_EML=ref_EML) %>%  mutate(ref_EML=ifelse(ref_EML=="EPI_Amnion","Epi",ref_EML)) %>% mutate(ref_EML=ifelse(ref_EML=="PriS_Amnion","PriS",ref_EML)) %>% mutate(ref_EML=ifelse(ref_EML=="EPI_PriS","Epi",ref_EML))  %>% mutate(ref_EML=ifelse(ref_EML=="EPI.PrE.INT","Epi",ref_EML))  %>% bind_rows(x.int %>% mutate(INT.ref_EML=ref_EML)  %>% mutate(ref_EML=ifelse(ref_EML=="EPI_Amnion","Amnion",ref_EML)) %>% mutate(ref_EML=ifelse(ref_EML=="PriS_Amnion","Amnion",ref_EML))%>% mutate(ref_EML=ifelse(ref_EML=="EPI_PriS","PriS",ref_EML))   %>% mutate(ref_EML=ifelse(ref_EML=="EPI.PrE.INT","Endoderm",ref_EML))) %>% select(query_cell,ref_EML,ref_pj,INT.ref_EML) %>% unique()
    
    x.int.dub <- x.int  %>% mutate(ref_EML=ifelse(ref_EML=="EPI_Amnion","Epi",ref_EML)) %>% mutate(ref_EML=ifelse(ref_EML=="PriS_Amnion","PriS",ref_EML)) %>% mutate(ref_EML=ifelse(ref_EML=="EPI_PriS","Epi",ref_EML))  %>% mutate(ref_EML=ifelse(ref_EML=="EPI.PrE.INT","Epi",ref_EML))  %>% bind_rows(x.int %>% mutate(ref_EML=ifelse(ref_EML=="EPI_Amnion","Amnion",ref_EML)) %>% mutate(ref_EML=ifelse(ref_EML=="PriS_Amnion","Amnion",ref_EML))%>% mutate(ref_EML=ifelse(ref_EML=="EPI_PriS","PriS",ref_EML))   %>% mutate(ref_EML=ifelse(ref_EML=="EPI.PrE.INT","Endoderm",ref_EML))) %>% mutate(nc_pva=mean_pval*nc_hits)%>% group_by(query_cell,query_EML,query_pj,ref_pj,ref_EML) %>% summarise(nc_hits=sum(nc_hits),nc_pval=sum(nc_pva),min_pval=min(min_pval),n_hits=mean(n_hits)) %>% ungroup() %>% mutate(majority_frac=nc_hits/n_hits,mean_pval=nc_pval/nc_hits) %>% select(-nc_pval)
    x.unique <- x.unique %>% group_by(query_cell,ref_pj) %>% top_n(1,-1*mean_pval) %>% top_n(1,nc_hits) %>% ungroup()
    
    x.unique <- x.unique %>% inner_join(x.int.dub,by = c("query_cell", "query_EML", "query_pj", "ref_pj", "ref_EML")) %>% mutate(nc_hits=nc_hits.x+nc_hits.y,n_hits=n_hits.x,mean_pval=(mean_pval.x*nc_hits.x+mean_pval.y*nc_hits.y)/(nc_hits.x+nc_hits.y),min_pval=min(min_pval.x,min_pval.y),majority_frac=nc_hits/n_hits) %>% select(-c(nc_hits.x,mean_pval.x,n_hits.x,majority_frac.x,nc_hits.y,mean_pval.y,n_hits.y,majority_frac.y,min_pval.x,min_pval.y)) %>% bind_rows( x.unique %>% anti_join(x.int.dub,by = c("query_cell", "query_EML", "query_pj", "ref_pj", "ref_EML")) )
    x.out <- x.int %>% anti_join(x.int.type %>% semi_join(x.unique,by = c("query_cell", "ref_EML", "ref_pj")) %>% mutate(ref_EML=INT.ref_EML), by = c("query_cell", "ref_pj", "ref_EML")) %>% bind_rows(x.unique) 
    
  }else if (nrow(x.int) ==0) {
    x.out <- x.unique
  }else if (nrow(x.uinque) ==0) {
    x.out <- x.int
  }
  x.mod <- x.out %>% bind_rows(x.basic)
  
 return(x.mod) 
}


FunHandleINT_min <- function(x) {
  x.int <- x %>% filter(ref_EML %in% c("EPI_Amnion","EPI.PrE.INT","PriS_Amnion"))
  x.unique <-  x %>% filter(!ref_EML %in% c("EPI_Amnion","EPI.PrE.INT","PriS_Amnion")) %>% filter(ref_EML %in% c("Amnion","PriS","Epi","Endoderm"))
  x.basic <- x  %>% filter(!ref_EML %in% c("EPI_Amnion","EPI.PrE.INT","PriS_Amnion")) %>% filter(!ref_EML %in% c("Amnion","PriS","Epi","Endoderm"))
  
  x.int1 <- x.int %>% mutate(ref_EML=ifelse(ref_EML=="EPI_Amnion","Epi",ref_EML)) %>% mutate(ref_EML=ifelse(ref_EML=="PriS_Amnion","PriS",ref_EML))  %>% mutate(ref_EML=ifelse(ref_EML=="EPI.PrE.INT","Epi",ref_EML))  
  x.int2 <- x.int %>% mutate(ref_EML=ifelse(ref_EML=="EPI_Amnion","Amnion",ref_EML)) %>% mutate(ref_EML=ifelse(ref_EML=="PriS_Amnion","Amnion",ref_EML))  %>% mutate(ref_EML=ifelse(ref_EML=="EPI.PrE.INT","Endoderm",ref_EML))  
  
  x.unique.int1 <-  x.unique %>% bind_rows(x.int1) %>% group_by(query_cell,query_EML,query_pj,ref_pj,ref_EML,n_hits) %>% summarise(nc_hits=sum(nc_hits),min_pval=min(min_pval),majority_frac=sum(majority_frac))
  x.unique.int2 <-  x.unique %>% bind_rows(x.int2) %>% group_by(query_cell,query_EML,query_pj,ref_pj,ref_EML,n_hits) %>% summarise(nc_hits=sum(nc_hits),min_pval=min(min_pval),majority_frac=sum(majority_frac))
  
  
  x.mod <- x.unique.int1 %>% bind_rows(x.unique.int2) %>% group_by(query_cell,query_EML,query_pj,ref_pj,ref_EML) %>% top_n(1,nc_hits)  %>% ungroup() %>% unique() %>% bind_rows(x.basic)
  
  return(x.mod) 
}
FunForMatAnno <- function(x) {
  x[x %in% c("EPI","Epi")] ="Epiblast"
  x[x %in% c("EarlyTE","medium_TE","Pre.ST","ST","TB.late","TB.medium3","TB.early","TB.medium1","TB.medium2","early_TE","EVT","late_TE","Pre.EVT","TB.apoptosis")] ="TE"
  x[x %in% c("PE","PrE")] = "Endoderm"
  x[x %in% c("B1_B2")] = "EarlyBlastocyst"
  x[x %in% c("IB_Epi","EB_ELC","nicolBla_ELC")]="ELC"
  x[x %in% c("IB_TE","EB_TLC","nicolBla_TLC")]="TLC"
  x[x %in% c("IB_PE","EB_HLC","nicolBla_HLC")]="HLC"
  x[x %in% c("IB_IM1","IB_IM3","IB_NR","IB_IM2","IB_uc","EB_U6","EB_UI13","EB_U5","EB_U10","Unknown","Trans","nicolBla_nc")]="Undefined"
  x[x %in% c("ambiguous","NoSigHits")] <- "Uncertain"
  return(x)
}

APUPJ <- function(temp,PJ) {
  temp.plot <- list()
  for ( n in unique(temp$EML[temp$pj==PJ])) {
    temp.sel.cell <- temp %>% filter(pj==PJ & EML==n) %>% pull(cell)
    temp.plot[[n]] <-ggplot()+geom_point(temp %>% filter(!cell %in% temp.sel.cell),mapping=aes(x=X_umap1,y=X_umap2),size=0.5,col="grey") + geom_point(temp %>% filter(cell %in% temp.sel.cell),mapping=aes(x=X_umap1,y=X_umap2),size=1,col="red")+theme_classic()+NoAxes()+NoLegend()+ggtitle(n)
    
  }
  print(cowplot::plot_grid(plotlist=temp.plot))
}

FunCSUMAP <- function(data.temp.UMAP,label1,label2) {
  temp.plot <-list()
  if ("Amnion(human)" %in% data.temp.UMAP$rename_EML) {
    temp.plot$m1 <-   ggplot()+geom_point(data.temp.UMAP %>% filter(rename_EML!="Amnion(human)") %>% filter(spec=="human"),mapping=aes(x=UMAP_1,y=UMAP_2,color=rename_EML,shape=spec),size=0.5)+geom_point(data.temp.UMAP %>% filter(rename_EML!="Amnion(NHP)") %>% filter(spec=="NHP"),mapping=aes(x=UMAP_1,y=UMAP_2,color=rename_EML,shape=spec),size=0.5)+geom_point(data.temp.UMAP %>% filter(rename_EML %in% c("Amnion(human)")) ,mapping=aes(x=UMAP_1,y=UMAP_2),shape=21,size=3,fill=cross.lineage.col.set["Amnion(human)"],color="grey",)+geom_point(data.temp.UMAP %>% filter(rename_EML %in% c("Amnion(NHP)")) ,mapping=aes(x=UMAP_1,y=UMAP_2),shape=24,size=3,fill=cross.lineage.col.set["Amnion(NHP)"],color="grey",alpha=0.75)+geom_text(data=data.temp.UMAP %>% group_by(rename_EML) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3 ) + theme_classic()+NoAxes()+scale_color_manual(values=cross.lineage.col.set)+scale_shape_manual(values = species.shape.set)+NoAxes()+ggtitle(paste(label1,label2,sep="+"))+theme(plot.title = element_text(hjust=0.5))+NoLegend()
    umap.xlim=ggplot_build(temp.plot$m1)$layout$panel_scales_x[[1]]$range$range
    umap.ylim=ggplot_build(temp.plot$m1)$layout$panel_scales_y[[1]]$range$range
    temp.plot$m2 <-    ggplot()+geom_point(data.temp.UMAP %>% filter(rename_EML!="Amnion(human)") %>% filter(spec=="human"),mapping=aes(x=UMAP_1,y=UMAP_2,color=rename_EML,shape=spec),size=0.5)+geom_point(data.temp.UMAP %>% filter(rename_EML!="Amnion(NHP)") %>% filter(spec=="NHP"),mapping=aes(x=UMAP_1,y=UMAP_2,color=rename_EML,shape=spec),size=0.5)+geom_point(data.temp.UMAP %>% filter(rename_EML %in% c("Amnion(human)")) ,mapping=aes(x=UMAP_1,y=UMAP_2),shape=21,size=3,fill=cross.lineage.col.set["Amnion(human)"],color="grey",)+geom_point(data.temp.UMAP %>% filter(rename_EML %in% c("Amnion(NHP)")) ,mapping=aes(x=UMAP_1,y=UMAP_2),shape=24,size=3,fill=cross.lineage.col.set["Amnion(NHP)"],color="grey",alpha=0.75)+theme_classic()+NoAxes()+scale_color_manual(values=cross.lineage.col.set)+scale_shape_manual(values = species.shape.set)+NoAxes()+NoLegend()+ggtitle(paste(label1,label2,sep="+"))+theme(plot.title = element_text(hjust=0.5))+xlim(umap.xlim)+ylim(umap.ylim)
    temp.plot$m3 <-  ggplot()+geom_text(data=data.temp.UMAP %>% group_by(rename_EML) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3 ) +theme_classic()+NoAxes()+scale_color_manual(values=cross.lineage.col.set)+scale_shape_manual(values = species.shape.set)+NoAxes()+NoLegend()+ggtitle(paste(label1,label2,sep="+"))+theme(plot.title = element_text(hjust=0.5))+xlim(umap.xlim)+ylim(umap.ylim)
    
  }else{
    temp.plot$m1 <-   ggplot()+geom_point(data.temp.UMAP  %>% filter(spec=="human"),mapping=aes(x=UMAP_1,y=UMAP_2,color=rename_EML,shape=spec),size=0.5)+geom_point(data.temp.UMAP %>% filter(rename_EML!="Amnion(NHP)") %>% filter(spec=="NHP"),mapping=aes(x=UMAP_1,y=UMAP_2,color=rename_EML,shape=spec),size=0.5) +geom_point(data.temp.UMAP %>% filter(rename_EML %in% c("Amnion(NHP)")) ,mapping=aes(x=UMAP_1,y=UMAP_2),shape=24,size=3,fill=cross.lineage.col.set["Amnion(NHP)"],color="grey")+geom_text(data=data.temp.UMAP %>% group_by(rename_EML) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3 )+theme_classic()+NoAxes()+scale_color_manual(values=cross.lineage.col.set)+scale_shape_manual(values = species.shape.set)+NoAxes()+ggtitle(paste(label1,label2,sep="+"))+theme(plot.title = element_text(hjust=0.5))+NoLegend()
    umap.xlim=ggplot_build(temp.plot$m1)$layout$panel_scales_x[[1]]$range$range
    umap.ylim=ggplot_build(temp.plot$m1)$layout$panel_scales_y[[1]]$range$range
    temp.plot$m2 <-    ggplot()+geom_point(data.temp.UMAP %>% filter(spec=="human"),mapping=aes(x=UMAP_1,y=UMAP_2,color=rename_EML,shape=spec),size=0.5)+geom_point(data.temp.UMAP %>% filter(rename_EML!="Amnion(NHP)") %>% filter(spec=="NHP"),mapping=aes(x=UMAP_1,y=UMAP_2,color=rename_EML,shape=spec),size=0.5)+geom_point(data.temp.UMAP %>% filter(rename_EML %in% c("Amnion(NHP)")) ,mapping=aes(x=UMAP_1,y=UMAP_2),shape=24,size=3,fill=cross.lineage.col.set["Amnion(NHP)"],color="grey")+theme_classic()+NoAxes()+scale_color_manual(values=cross.lineage.col.set)+scale_shape_manual(values = species.shape.set)+NoAxes()+NoLegend()+ggtitle(paste(label1,label2,sep="+"))+theme(plot.title = element_text(hjust=0.5))+xlim(umap.xlim)+ylim(umap.ylim)
    temp.plot$m3 <-  ggplot()+geom_text(data=data.temp.UMAP %>% group_by(rename_EML) %>% summarise(UMAP_1=median(UMAP_1),UMAP_2=median(UMAP_2)),mapping=aes(x=UMAP_1,y=UMAP_2,label=rename_EML),size=3 ) +theme_classic()+NoAxes()+scale_color_manual(values=cross.lineage.col.set)+scale_shape_manual(values = species.shape.set)+NoAxes()+NoLegend()+ggtitle(paste(label1,label2,sep="+"))+theme(plot.title = element_text(hjust=0.5))+xlim(umap.xlim)+ylim(umap.ylim)
  }
  temp.plot$m4 <-  ggplot()+geom_point(data.temp.UMAP %>% filter(spec=="human"),mapping=aes(x=UMAP_1,y=UMAP_2,shape=spec,col=spec),size=0.5)+geom_point(data.temp.UMAP %>% filter(spec=="NHP"),mapping=aes(x=UMAP_1,y=UMAP_2,shape=spec,col=spec),size=0.5)+theme_classic()+NoAxes()+scale_color_manual(values=c("human"="red","NHP"="grey"))+scale_shape_manual(values = species.shape.set)+NoAxes()+NoLegend()+ggtitle(paste(label1,label2,sep="+"))+theme(plot.title = element_text(hjust=0.5))+xlim(umap.xlim)+ylim(umap.ylim)
  
  
  #ggrepel::geom_text_repel
  return(temp.plot)
}


FunDEG <- function(data,cl,c1,c2,psd.counts=1) {
  data.sub=subset(data,cells=c(cl[[c1]],cl[[c2]]))
  ident=as.factor(c(rep("c1",length(cl[[c1]])),rep("c2",length(cl[[c2]]))))
  names(ident) =c(cl[[c1]],cl[[c2]])
  Idents(data.sub)=ident[colnames(data.sub)]
  suppressMessages(G1G2.DEG <- FindMarkers(data.sub,ident.1="c1",ident.2="c2",verbose=F,test.use="roc",assay="RNA",slot="data",pseudocount.use=psd.counts) %>%tibble::rownames_to_column("gene") %>%  tbl_df())
  return(G1G2.DEG)
}


FunRF_FindAllMarkers_para <- function(data,ct=0.3,assay="RNA",slot="data") {
  temp.id <- unique(as.vector(Idents(data)))
  temp.gene <- list()
  for (n in c(1:(length(temp.id)-1))) {
    for (m in c(2:length(temp.id))) {
      if (m >n) {
        temp.gene[[paste(m,n,sep="_")]] <- paste(m,n,sep="_")
      }
    }
  }
  temp.c <- names(temp.gene)
  temp.gene <- foreach (x=temp.gene,w=temp.c, .combine=c) %dopar% {
    rv=list()
    rv[[w]]=FindMarkers(data,ident.1=temp.id[as.numeric(unlist(strsplit(x,"_") )[2])],ident.2=temp.id[as.numeric(unlist(strsplit(x,"_") )[1])],test.use = "roc",verbose=F,assay=assay,slot=slot) %>% tibble::rownames_to_column("gene") %>% tbl_df() %>%  mutate(BG=ifelse(myAUC>=0.5,temp.id[as.numeric(unlist(strsplit(x,"_") )[2])],temp.id[as.numeric(unlist(strsplit(x,"_") )[1])])) %>%  mutate(EG=ifelse(myAUC>=0.5,temp.id[as.numeric(unlist(strsplit(x,"_") )[1])],temp.id[as.numeric(unlist(strsplit(x,"_") )[2])]))
    rv
  }
  
  temp.gene <- do.call("bind_rows",temp.gene)
  temp.out <- list()
  temp.out$detail <- temp.gene
  temp.out$stat <- temp.gene  %>% filter(power > ct)  %>% group_by(gene,BG) %>% summarise(EGset=paste(EG,collapse=","),nEGset=n()) %>% group_by(gene,EGset,nEGset,) %>% summarise(BGset=paste(BG,collapse=","),nBGset=n())
  temp.out$sig <- temp.gene %>% filter(power > ct) %>% group_by(gene,BG) %>% summarise(n=n(),power=mean(power)) %>% inner_join(temp.out$stat %>% filter(nBGset==1) %>% ungroup() %>% select(BGset,nEGset) %>% unique()%>% group_by(BGset) %>% top_n(1,nEGset) %>% ungroup()%>% rename(n=nEGset,BG=BGset),by=c("BG","n")) %>% mutate(n=(length(temp.id)-n)) %>% arrange(desc(power)) %>% ungroup() %>% mutate(set=BG) %>%  select(-BG) %>% mutate(NP="pos")%>% arrange(set,NP,desc(power))
  return(
    temp.out
  )
}

APPJ <- function(data.temp,PJ) {
  temp.plot <- list()
  for ( n in unique(data.temp@meta.data$EML[data.temp@meta.data$pj==PJ])) {
    temp.plot[[n]] <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj==PJ & data.temp@meta.data$EML==n])+theme_void()+NoLegend()+ggtitle(n)+theme(plot.title = element_text(hjust=0.5,face="bold"))
  }
  print(cowplot::plot_grid(plotlist=temp.plot))
}

APnGuo <- function(data.temp) {
  temp.plot <- list()
  #temp.plot$raw <- DimPlot(data.temp,label=T)+theme_void()+NoLegend()
  #temp.plot$EML <- DimPlot(data.temp,group.by="EML",label=T)+theme_void()+NoLegend()
  for ( n in unique(data.temp@meta.data$EML[data.temp@meta.data$pj=="nBGuo"])) {
    temp.plot[[n]] <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="nBGuo" & data.temp@meta.data$EML==n])+theme_void()+NoLegend()+ggtitle(n)+theme(plot.title = element_text(hjust=0.5,face="bold"))
  }
  print(cowplot::plot_grid(plotlist=temp.plot))
}
AP4 <- function(data.temp) {
  temp.plot <- list()
  temp.plot$TE <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$EML %in% c("TE","STB","EVT","CTB") ])+theme_void()+NoLegend()+ggtitle("TE")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$AMLC <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$EML %in% c("Tsw-AMLC","AMLC") ])+theme_void()+NoLegend()+ggtitle("AMLC")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$D2_TE1 <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$EML %in% c("IB_TE") ])+theme_void()+NoLegend()+ggtitle("IB_TE")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  #temp.plot$D2_TE2 <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$EML %in% c("EB_TLC") ])+theme_void()+NoLegend()+ggtitle("EB_TLC")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  #temp.plot$D2_Epi <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$EML %in% c("IB_Epi","EB_ELC") ])+theme_void()+NoLegend()+ggtitle("EpiLike")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  #temp.plot$D2_PE <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$EML %in% c("IB_PE","EB_HLC") ])+theme_void()+NoLegend()+ggtitle("PELike")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$EM <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$cellType=="EM" & data.temp@meta.data$pj !="JPF2019" ])+theme_void()+NoLegend()+ggtitle("EM only")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$IBD2 <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="IBD2"])+theme_void()+NoLegend()+ggtitle("IBD2")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$EBD2 <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="EBD2"])+theme_void()+NoLegend()+ggtitle("EBD2")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$nBGuo <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="nBGuo"])+theme_void()+NoLegend()+ggtitle("nBGuo")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  print(cowplot::plot_grid(plotlist=temp.plot))
}

AP3 <- function(data.temp) {
  temp.plot <- list()
  temp.plot$E10X <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$seqType=="10X" & data.temp@meta.data$EML=="EPSC" ])+theme_void()+NoLegend()+ggtitle("EPSC(10X)")
  temp.plot$ES <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$seqType=="smt2" & data.temp@meta.data$EML=="EPSC" ])+theme_void()+NoLegend()+ggtitle("EPSC(smt2)")
  temp.plot$N10X <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$seqType=="10X" & data.temp@meta.data$EML=="Naive" ])+theme_void()+NoLegend()+ggtitle("Naive(10X)")
  temp.plot$NS <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$seqType=="smt2" & data.temp@meta.data$EML=="Naive" ])+theme_void()+NoLegend()+ggtitle("Naive(smt2)")
  temp.plot$PS <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$seqType=="smt2" & data.temp@meta.data$EML=="Primed" ])+theme_void()+NoLegend()+ggtitle("Primed(smt2)")
  temp.plot$P10X <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$seqType=="10X" & data.temp@meta.data$EML=="Primed" ])+theme_void()+NoLegend()+ggtitle("Primed(10X)")
   temp.plot$TswAMLC <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$EML=="Tsw-AMLC" ])+theme_void()+NoLegend()+ggtitle("Tsw-AMLC")+theme(plot.title = element_text(hjust=0.5,face="bold"))
   temp.plot$IB_TE <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$EML %in% c("IB_TE","EB_TLC") ])+theme_void()+NoLegend()+ggtitle("TELike")+theme(plot.title = element_text(hjust=0.5,face="bold"))
   temp.plot$IB_TE <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$EML %in% c("IB_Epi","EB_ELC") ])+theme_void()+NoLegend()+ggtitle("EpiLike")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$EM <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$cellType=="EM" & data.temp@meta.data$pj !="JPF2019" ])+theme_void()+NoLegend()+ggtitle("EM only")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$IBD2 <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="IBD2"])+theme_void()+NoLegend()+ggtitle("IBD2")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$EBD2 <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="EBD2"])+theme_void()+NoLegend()+ggtitle("EBD2")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  print(cowplot::plot_grid(plotlist=temp.plot))
}
AP <- function(data.temp) {
  temp.plot <- list()
  #temp.plot$raw <- DimPlot(data.temp,label=T)+theme_void()+NoLegend()
  temp.plot$EML <- DimPlot(data.temp,group.by="EML",label=T)+theme_void()+NoLegend()
  print(cowplot::plot_grid(plotlist=temp.plot))
  
}

AP2 <- function(data.temp) {
  temp.plot <- list()
  temp.plot$E10X <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$seqType=="10X" & data.temp@meta.data$EML=="EPSC" ])+theme_void()+NoLegend()+ggtitle("EPSC(10X)")
  temp.plot$ES <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$seqType=="smt2" & data.temp@meta.data$EML=="EPSC" ])+theme_void()+NoLegend()+ggtitle("EPSC(smt2)")
  temp.plot$N10X <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$seqType=="10X" & data.temp@meta.data$EML=="Naive" ])+theme_void()+NoLegend()+ggtitle("Naive(10X)")
  temp.plot$NS <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$seqType=="smt2" & data.temp@meta.data$EML=="Naive" ])+theme_void()+NoLegend()+ggtitle("Naive(smt2)")
  temp.plot$PS <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$seqType=="smt2" & data.temp@meta.data$EML=="Primed" ])+theme_void()+NoLegend()+ggtitle("Primed(smt2)")
  temp.plot$P10X <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$seqType=="10X" & data.temp@meta.data$EML=="Primed" ])+theme_void()+NoLegend()+ggtitle("Primed(10X)")
  temp.plot$BAP12 <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$EML %in% c("BAP_sub1","BAP_sub2") ])+theme_void()+NoLegend()+ggtitle("BAP_sub12")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$TswAMLC <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$EML=="Tsw-AMLC" ])+theme_void()+NoLegend()+ggtitle("Tsw-AMLC")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$EM <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$cellType=="EM" & data.temp@meta.data$pj !="JPF2019" ])+theme_void()+NoLegend()+ggtitle("EM only")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$OkaTsc <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="TscOka"])+theme_void()+NoLegend()+ggtitle("TscOka")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$BAP <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="BAPdiff"])+theme_void()+NoLegend()+ggtitle("BAPdiff")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  print(cowplot::plot_grid(plotlist=temp.plot))
}


APEpi <- function(data.temp) {
  temp.plot <- list()
  temp.plot$SPHEPI.early <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="SPH2016" & data.temp@meta.data$EML=="EPI"& data.temp@meta.data$devTime %in% c("E5") ])+theme_void()+NoLegend()+ggtitle("SPH_epi(D5)")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$SPHEPI.late <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="SPH2016" & data.temp@meta.data$EML=="EPI"& data.temp@meta.data$devTime %in% c("E6","E7") ])+theme_void()+NoLegend()+ggtitle("SPH_epi(D6,D7)")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$D3EPI.early <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="D3post" & data.temp@meta.data$EML=="EPI" &data.temp@meta.data$devTime %in% c("E6","E7")])+theme_void()+NoLegend()+ggtitle("3D_epi(D6,D7)")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$D3EPI.late <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="D3post" & data.temp@meta.data$EML=="EPI" & (!data.temp@meta.data$devTime %in% c("E6","E7"))])+theme_void()+NoLegend()+ggtitle("3D_epi(D8-D14)")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$SRTEPI <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="SRT10X" & data.temp@meta.data$EML=="EPI" ])+theme_void()+NoLegend()+ggtitle("SRT_epi")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$CSEPI <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="CS7" & data.temp@meta.data$EML=="Epiblast" ])+theme_void()+NoLegend()+ggtitle("CS_epi")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  print(cowplot::plot_grid(plotlist=temp.plot))
}


APPE <- function(data.temp) {
  temp.plot <- list()
  temp.plot$SPHTE <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="SPH2016" & data.temp@meta.data$EML=="PE" ])+theme_void()+NoLegend()+ggtitle("SPH_PE")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$D3PE.early <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="D3post" & data.temp@meta.data$EML=="PrE" &data.temp@meta.data$devTime %in% c("E6","E7")])+theme_void()+NoLegend()+ggtitle("3D_PE(D6,D7)")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$D3PE.late <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="D3post" & data.temp@meta.data$EML=="PrE" & (!data.temp@meta.data$devTime %in% c("E6","E7"))])+theme_void()+NoLegend()+ggtitle("3D_PE(post-imp)")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$SRTPE <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="SRT10X" & data.temp@meta.data$EML=="PE" ])+theme_void()+NoLegend()+ggtitle("SRT_PE")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$CSPE <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="CS7" & data.temp@meta.data$EML=="Endoderm" ])+theme_void()+NoLegend()+ggtitle("CS_Endoderm")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  print(cowplot::plot_grid(plotlist=temp.plot))
}

FunMaSF <- function(temp,n) {
  set.seed(456)
  return(head(temp[sample(nrow(temp)),],n))
}
