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

AP4 <- function(data.temp) {
  temp.plot <- list()
  temp.plot$TE <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$EML %in% c("TE","STB","EVT","CTB") ])+theme_void()+NoLegend()+ggtitle("TE")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$TswAMLC <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$EML=="Tsw-AMLC" ])+theme_void()+NoLegend()+ggtitle("Tsw-AMLC")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$D2_TE1 <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$EML %in% c("IB_TE") ])+theme_void()+NoLegend()+ggtitle("IB_TE")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$D2_TE2 <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$EML %in% c("EB_TLC") ])+theme_void()+NoLegend()+ggtitle("EB_TLC")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$D2_Epi <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$EML %in% c("IB_Epi","EB_ELC") ])+theme_void()+NoLegend()+ggtitle("EpiLike")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$D2_PE <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$EML %in% c("IB_PE","EB_HLC") ])+theme_void()+NoLegend()+ggtitle("PELike")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$EM <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$cellType=="EM" & data.temp@meta.data$pj !="JPF2019" ])+theme_void()+NoLegend()+ggtitle("EM only")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$IBD2 <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="IBD2"])+theme_void()+NoLegend()+ggtitle("IBD2")+theme(plot.title = element_text(hjust=0.5,face="bold"))
  temp.plot$EBD2 <- DimPlot(data.temp,cells.highlight=colnames(data.temp)[data.temp@meta.data$pj=="EBD2"])+theme_void()+NoLegend()+ggtitle("EBD2")+theme(plot.title = element_text(hjust=0.5,face="bold"))
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
