#' To combine all dataset
# R3.6
rm(list=ls())
rewrite=FALSE

condaENV <- "/home/chenzh/miniconda3/envs/R3.6"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(edgeR))
suppressMessages(library(ggplot2))
suppressMessages(library(cowplot))
suppressMessages(library(pheatmap))
suppressMessages(library(umap))
suppressMessages(library(Seurat))
suppressMessages(library(readxl))

setwd("/home/chenzh/My_project/CheckBlastoids")

#' loading previous annotation of JP
load("tmp_data/Prev.old.M.fromJP.Rdata",verbose = T)
#save(old.M.all,file=)
# from load("tmp_data/AllData.clean.Robj",verbose = T) old.M.all=M.all


counts <- list()
rpkms <- list()
tpms <- list()
metas <- list()


#' deal with Sophie2016
pj="SPH2016"
cp.count <- read.delim("~/My_project/CheckBlastoids/tmp_data/CP_data/merge.rsem_counts.csv",stringsAsFactors = F,head=T,row.names = 1,sep=",")  %>% tibble::rownames_to_column("Gene")%>% tbl_df()
temp <- read.delim("/home/chenzh/My_project/ScRNA_analysis/data/Published/E-MTAB-3929.sdrf.txt",stringsAsFactors = F,head=T) 
rownames(temp)=temp$Comment.ENA_RUN.
colnames(cp.count) <- c("Gene",temp[colnames(cp.count)[-1],] %>% pull("Source.Name"))
counts[[pj]] <- cp.count             

rpkms[[pj]] <- read.delim("~/My_project/CheckBlastoids/tmp_data/CP_data/merge.rsem_rpkm.csv",stringsAsFactors = F,head=T,row.names = 1,sep=",")  %>% tibble::rownames_to_column("Gene")%>% tbl_df()   
colnames(rpkms[[pj]]) <- c("Gene",temp[colnames(rpkms[[pj]])[-1],] %>% pull("Source.Name"))

tpms[[pj]] <- read.delim("~/My_project/CheckBlastoids/tmp_data/CP_data/merge.rsem_tpm.csv",stringsAsFactors = F,head=T,row.names = 1,sep=",")  %>% tibble::rownames_to_column("Gene")%>% tbl_df()     
colnames(tpms[[pj]]) <- c("Gene",temp[colnames(tpms[[pj]])[-1],] %>% pull("Source.Name"))


metas[[pj]] <- old.M.all %>% filter(cell %in% colnames(counts[[pj]]))
metas[[pj]] <- metas[[pj]] %>% mutate(libsize=colSums(counts[[pj]] %>% select(-Gene) )[metas[[pj]]$cell]) %>% mutate(nGene=colSums((counts[[pj]]  %>% select(-Gene)) >0)[metas[[pj]]$cell])   %>% mutate(pj=pj) #%>% mutate(cell=paste0(pj,"_",cell))

#(3) Cell lineage  from another published paper #Integrated analysis of single-cell embryo data yields a unified transcriptome signature for the human pre-implantation epiblas
cell_lineage_data <- read_xlsx("doc/stirparo2018_tableS4.xlsx", sheet = 1)
cell_lineage_data <- cell_lineage_data[cell_lineage_data$Study == "Petropoulos et al., 2016 (ERP012552)", ] %>% rename('new_lineage'='Revised lineage (this study)') %>% mutate(new_lineage=gsub("epiblast","EPI",new_lineage)) %>% mutate(new_lineage=gsub("intermediate","INT",new_lineage))%>% mutate(new_lineage=gsub("primitive_endoderm","PE",new_lineage))%>% mutate(new_lineage=gsub("Inner cell mass","ICM",new_lineage)) %>% mutate(new_lineage=gsub("trophectoderm","TE",new_lineage))%>% mutate(new_lineage=gsub("undefined","EM_NA",new_lineage))
 # 1,481 cells
cell_lineage_data$Cell <- gsub("_", ".", cell_lineage_data$Cell)
#metas$SPH2016 <- metas$SPH2016 %>% left_join(cell_lineage_data %>% rename(cell=Cell) %>% select(cell,new_lineage),by="cell") %>% mutate(EML=new_lineage) %>% mutate(SID=paste(devTime,EML,sep="_")) %>% select(-new_lineage)
rm(cell_lineage_data)




# GSE136447 (SS2) #3D postimplantation
pj="D3post"

temp <-  read.csv("~/My_project/CheckBlastoids/data/GSE136447/SraRunInfo.csv",stringsAsFactors = F,head=T) %>% tbl_df() %>% select(Run,SampleName) %>% left_join(read.delim("~/My_project/CheckBlastoids/data/GSE136447/GSE136447.tsv",stringsAsFactors = F,head=T) %>% tbl_df() %>% select(Accession,Title) %>% rename(SampleName=Accession),by="SampleName") %>% mutate(cell=gsub("Embryo_","",Title)) %>% left_join(read.delim("~/My_project/CheckBlastoids/data/GSE136447/41586_2019_1875_MOESM10_ESM.tsv",stringsAsFactors = F,head=T) %>% tbl_df() %>% rename(cell=SampleID),by="cell") %>% tibble::column_to_rownames("Run")


counts[[pj]] <- read.delim("~/My_project/CheckBlastoids/tmp_data/GSE136447/merge.rsem_counts.csv",stringsAsFactors = F,head=T,row.names = 1,sep=",")  %>% tibble::rownames_to_column("Gene")%>% tbl_df()           
colnames(counts[[pj]]) <- c("Gene",temp[colnames(counts[[pj]])[-1],] %>% pull("cell"))


metas[[pj]] <- temp %>% tbl_df()  %>% select(-SampleName,-Title)%>% rename(devTime=Day)%>% mutate(devTime=gsub("^D","E",devTime)) %>% rename(subCT=EmbryoID,EML=Group) %>% mutate(SID=paste(devTime,EML,sep="_")) %>% mutate(batch="batch7",seqType="smt2",cellType="EM")

metas[[pj]] <-  metas[[pj]] %>% mutate(libsize=colSums(counts[[pj]] %>% select(-Gene) )[metas[[pj]]$cell]) %>% mutate(nGene=colSums((counts[[pj]]  %>% select(-Gene)) >0)[metas[[pj]]$cell], pj=pj)                                                       

rpkms[[pj]] <- read.delim("~/My_project/CheckBlastoids/tmp_data/GSE136447/merge.rsem_rpkm.csv",stringsAsFactors = F,head=T,row.names = 1,sep=",")  %>% tibble::rownames_to_column("Gene")%>% tbl_df()   
colnames(rpkms[[pj]]) <- c("Gene",temp[colnames(rpkms[[pj]])[-1],] %>% pull("cell"))

tpms[[pj]] <- read.delim("~/My_project/CheckBlastoids/tmp_data/GSE136447/merge.rsem_tpm.csv",stringsAsFactors = F,head=T,row.names = 1,sep=",")  %>% tibble::rownames_to_column("Gene")%>% tbl_df()     
colnames(tpms[[pj]]) <- c("Gene",temp[colnames(tpms[[pj]])[-1],] %>% pull("cell"))



# GSE66507 (SS2) #blakeley
pj="Blakeley"

counts[[pj]] <- read.delim("~/My_project/CheckBlastoids/tmp_data/Blakeley_Data/merge.rsem_counts.csv",stringsAsFactors = F,head=T,row.names = 1,sep=",")  %>% tibble::rownames_to_column("Gene")%>% tbl_df()   

temp=colnames(counts[[pj]])[-1]
temp.em <- paste0("Bla_",strsplit(temp,"\\.")  %>% lapply(function(x) x[1]) %>% unlist())
temp.EML <- strsplit(temp,"\\d")  %>% lapply(function(x) x[1]) %>% unlist()
metas[[pj]] <- data.frame(cell=paste0("Bla_",temp),devTime="E7",subCT=temp.em,EML=temp.EML,SID=paste0("E7_",temp.EML),batch="batchBla",seqType="smt2",cellType="EM") %>% tbl_df() %>% mutate_all(as.vector)
colnames(counts[[pj]]) <- c("Gene",paste0("Bla_",temp))

metas[[pj]] <-  metas[[pj]] %>% mutate(libsize=colSums(counts[[pj]] %>% select(-Gene) )[metas[[pj]]$cell]) %>% mutate(nGene=colSums((counts[[pj]]  %>% select(-Gene)) >0)[metas[[pj]]$cell], pj=pj)                                                       

# GSE171820 (SS2) #Austin Smith
pj="nBGuo"#GSE171820

counts[[pj]] <- read.delim("~/My_project/CheckBlastoids/tmp_data/nBlastoids//merge.rsem_counts.csv",stringsAsFactors = F,head=T,row.names = 1,sep=",")  %>% tibble::rownames_to_column("Gene")%>% tbl_df()   
temp <- read.delim("~/My_project/CheckBlastoids/data/GSE171820/GSE171820.anno.txt",stringsAsFactors = F,sep="\t",head=F) %>% tbl_df() %>% mutate(V3=ifelse(V3=="Early Trophectoderm","EarlyTE",V3))%>% mutate(V3=ifelse(V3=="Inner Cell Mass","ICM",V3))%>% mutate(V3=ifelse(V3=="Inner Cell Mass-Trophectoderm Transition","ICM-TE_trans",V3))%>% mutate(V3=ifelse(V3=="Primitive Endoderm","PE",V3)) %>% mutate(V3=ifelse(V3=="Trophectoderm","TE",V3)) %>% mutate(V3=ifelse(V3=="Transitioning","Trans",V3)) %>% mutate(V3=ifelse(V2=="Blastoid" & V3=="EarlyTE", "EarlyTLC",V3)) %>% mutate(V3=ifelse(V2=="Blastoid" & V3=="Epiblast", "ELC",V3)) %>% mutate(V3=ifelse(V2=="Blastoid" & V3=="ICM", "ICMLC",V3)) %>% mutate(V3=ifelse(V2=="Blastoid" & V3=="PE", "HLC",V3)) %>% mutate(V3=ifelse(V2=="Blastoid" & V3=="TE", "TLC",V3))

metas[[pj]] <- data.frame(cell=paste0("nBG_",temp$V1),devTime=temp$V5,subCT=strsplit(temp$V1,"\\.") %>% lapply(function(x){paste(x[1],x[2],sep=".")}) %>% unlist(),EML=temp$V3,SID=paste(temp$V5,temp$V3,sep="_"),batch="batchnB",seqType="smt2",cellType=temp$V2) %>% tbl_df() %>% mutate_all(as.vector)
colnames(counts[[pj]]) <- c("Gene",paste0("nBG_",colnames(counts[[pj]])[-1]))

metas[[pj]] <-  metas[[pj]] %>% mutate(libsize=colSums(counts[[pj]] %>% select(-Gene) )[metas[[pj]]$cell]) %>% mutate(nGene=colSums((counts[[pj]]  %>% select(-Gene)) >0)[metas[[pj]]$cell], pj=pj)  %>% mutate(cellType=ifelse(cellType=="Blastocyst","EM",cellType))                                              


pj="nicolBla"
counts[[pj]] <- read.delim("/home/chenzh/My_project/CheckBlastoids/data/GSE177689/CR_RefSeq/merge.rsem_counts.csv",stringsAsFactors = F,head=T,row.names = 1,sep=",")  %>% tibble::rownames_to_column("Gene")%>% tbl_df()       
metas[[pj]] <-  read.delim("data/GSE177689/blastoid.H9.okae_bts5__scRNAseq.unfiltered__metadata.tsv",stringsAsFactors = F,head=T,sep="\t") %>% tbl_df() %>% mutate(time = case_when(grepl("24h",samplename) ~ "24h",grepl("60h",samplename) ~ "60h", grepl("96h",samplename) ~ "96h",TRUE ~ gsub("-.*","",samplename))) %>% tbl_df() %>% filter(nFeature_RNA > 2000 & percent.mt < 12.5) %>% rename(nGene=nFeature_RNA,devTime=time,cell=sampleid) %>% mutate(batch="nicolBla",seqType="smt2",cell=paste("nicolBla",cell,sep="_")) %>% mutate(EML=devTime) %>% select(cell,nGene,batch,EML,seqType,devTime) %>% filter(!cell %in% c("nicolBla_150680")) %>% mutate(cellType=ifelse(devTime %in% c("24h","60h","96h"),"blastoids",EML), SID=ifelse(devTime %in% c("24h","60h","96h"),paste0("blastoids_",devTime),paste0("nicolBla_",EML)))
# unable to download nicolBla_150680.sra #%>% inner_join(read.delim("data/GSE177689/GSE177689_series_matrix.trans.txt",stringsAsFactors = F,head=F,sep="\t") %>% tbl_df() %>% setNames(c("EML","cell","GSM")),by="cell") # the filtering from the paper
metas[[pj]] <-  metas[[pj]] %>% mutate(libsize=colSums(counts[[pj]] %>% select(-Gene) )[metas[[pj]]$cell]) %>% mutate(nGene=colSums((counts[[pj]]  %>% select(-Gene)) >0)[metas[[pj]]$cell], pj=pj,subCT=cellType)  %>% select(cell,devTime,subCT,EML,SID,batch,seqType,cellType,libsize,nGene,pj)


pj="JPF2019"
temp.list <- list()
temp.list[["H9AMN"]] <- Read10X(data.dir = "tmp_data/Amnion_EPSC/H9_Amnion/outs/filtered_feature_bc_matrix")
temp.list[["Pos48ELS"]] <- Read10X(data.dir = "tmp_data/Amnion_EPSC/Posterior48h_1/outs/filtered_feature_bc_matrix")

GI=rownames(temp.list[["H9AMN"]])
SI=names(temp.list)
temp <- NULL
for ( n in names(temp.list)) {
  colnames(temp.list[[n]]) <- paste(n,"10X",colnames(temp.list[[n]]),sep="_") 
  temp.list[[n]] <- temp.list[[n]] %>% as.data.frame() %>% tbl_df() #%>% tibble::rownames_to_column("Gene")
}
counts[[pj]] <- data.frame(Gene=GI) %>% tbl_df() %>% mutate(Gene=as.vector(Gene)) %>% bind_cols(do.call("bind_cols",temp.list)) ## each sample has the same rownames
rm(temp.list)

temp.list <- list()
CI=colnames(counts[[pj]] )

n="H9AMN"
temp.list[[n]] <- data.frame(SID=paste(n,"10X",sep="_"),batch=pj,seqType="10X",cellType="H9",devTime=n,EML=n,cell=CI[grepl(n,CI)],subCT="None") %>% tbl_df()
n="Pos48ELS"
temp.list[[n]] <- data.frame(SID=paste(n,"10X",sep="_"),batch=pj,seqType="10X",cellType="Pos48",devTime=n,EML=n,cell=CI[grepl(n,CI)],subCT="None") %>% tbl_df()


metas[[pj]] <-  bind_rows(do.call("bind_rows",temp.list))
metas[[pj]] <- metas[[pj]]  %>% mutate(libsize=colSums(counts[[pj]] %>% select(-Gene) )[metas[[pj]]$cell]) %>% mutate(nGene=colSums((counts[[pj]]  %>% select(-Gene)) >0)[metas[[pj]]$cell], pj=pj)  


# https://www.biorxiv.org/content/10.1101/2020.07.21.213512v1.full
pj="CS7" 

counts[[pj]] <- read.delim("~/My_project/CheckBlastoids/tmp_data/human_gastrula/merge.rsem_counts.csv",stringsAsFactors = F,head=T,row.names = 1,sep=",")  %>% tibble::rownames_to_column("Gene")%>% tbl_df()       

metas[[pj]] <- readRDS("data/Gastrulation/annot_umap.rds") %>% tbl_df() %>% rename(cell=cell_name,EML=cluster_id,subCT=sub_cluster) %>% mutate(EML=gsub(" ","_",EML)) %>%select(cell,EML,subCT) %>% mutate(cell=gsub("SS.","",cell),SID=paste("CS7",EML,sep="_"),batch="batchCS",seqType="smt2",cellType="EM",devTime="CS7") %>% select(SID,batch,seqType,cellType,devTime,EML,cell,subCT)

metas[[pj]] <-  metas[[pj]] %>% mutate(libsize=colSums(counts[[pj]] %>% select(-Gene))[metas[[pj]]$cell]) %>% mutate(nGene=colSums((counts[[pj]]  %>% select(-Gene)) >0)[metas[[pj]]$cell], pj=pj)                                                       


pj="IBD2"
temp.list <- list()
temp.list[["IBD2"]] <- Read10X(data.dir = "tmp_data/IBlastoids/out/IBD2/outs/filtered_feature_bc_matrix")


GI=rownames(temp.list[["IBD2"]])
SI=names(temp.list)
temp <- NULL
for ( n in names(temp.list)) {
  colnames(temp.list[[n]]) <- paste(n,"10X",colnames(temp.list[[n]]),sep="_") 
  temp.list[[n]] <- temp.list[[n]] %>% as.data.frame() %>% tbl_df() #%>% tibble::rownames_to_column("Gene")
}
counts[[pj]] <- data.frame(Gene=GI) %>% tbl_df() %>% mutate(Gene=as.vector(Gene)) %>% bind_cols(do.call("bind_cols",temp.list)) ## each sample has the same rownames
rm(temp.list)

temp.list <- list()
CI=colnames(counts[[pj]] )
for ( n in SI) {
  temp.list[[n]] <- data.frame(SID=paste(n,"10X",sep="_"),batch="batchIBD2",seqType="10X",cellType=n,devTime=n,EML=n,cell=CI[grepl(n,CI)],subCT="None") %>% tbl_df()
}

metas[[pj]] <-  bind_rows(do.call("bind_rows",temp.list))
metas[[pj]] <- metas[[pj]]  %>% mutate(libsize=colSums(counts[[pj]] %>% select(-Gene) )[metas[[pj]]$cell]) %>% mutate(nGene=colSums((counts[[pj]]  %>% select(-Gene)) >0)[metas[[pj]]$cell], pj=pj)  

pj="EBD2"
temp.list <- list()
temp.list[["LW36"]] <- Read10X(data.dir = "tmp_data/EBlastoids/out/LW36/outs/filtered_feature_bc_matrix")
temp.list[["LW60"]] <- Read10X(data.dir = "tmp_data/EBlastoids/out/LW60/outs/filtered_feature_bc_matrix")
temp.list[["LW61"]] <- Read10X(data.dir = "tmp_data/EBlastoids/out/LW61/outs/filtered_feature_bc_matrix")


GI=rownames(temp.list[["LW61"]])
SI=names(temp.list)
temp <- NULL
for ( n in names(temp.list)) {
  colnames(temp.list[[n]]) <- paste(n,"10X",colnames(temp.list[[n]]),sep="_") 
  temp.list[[n]] <- temp.list[[n]] %>% as.data.frame() %>% tbl_df() #%>% tibble::rownames_to_column("Gene")
}
counts[[pj]] <- data.frame(Gene=GI) %>% tbl_df() %>% mutate(Gene=as.vector(Gene)) %>% bind_cols(do.call("bind_cols",temp.list)) ## each sample has the same rownames
rm(temp.list)

temp.list <- list()
CI=colnames(counts[[pj]] )
for ( n in SI) {
  temp.list[[n]] <- data.frame(SID=paste(n,"10X",sep="_"),batch="batchEBD2",seqType="10X",cellType=n,devTime=n,EML=n,cell=CI[grepl(n,CI)],subCT="None") %>% tbl_df()
}

metas[[pj]] <-  bind_rows(do.call("bind_rows",temp.list))
metas[[pj]] <- metas[[pj]]  %>% mutate(libsize=colSums(counts[[pj]] %>% select(-Gene) )[metas[[pj]]$cell]) %>% mutate(nGene=colSums((counts[[pj]]  %>% select(-Gene)) >0)[metas[[pj]]$cell], pj=pj)  





# Do the QC based on the MT.percent and nGene
gtf.anno <- read.delim( "~/Genome/Human/refdata-cellranger-GRCh38-3.0.0/genes/gene.gtf.anno",stringsAsFactors = F,row.names = 5,head=F)
mt.gene <- gtf.anno$V6[gtf.anno$V1=="MT"]
ribo.gene <- c(gtf.anno$V6[grepl("^RPS",gtf.anno$V6)],gtf.anno$V6[grepl("^RPL",gtf.anno$V6)],gtf.anno$V6[grepl("^MRPL",gtf.anno$V6)],gtf.anno$V6[grepl("^MRPS",gtf.anno$V6)])




dup.gene <- names(table(gtf.anno$V6))[table(gtf.anno$V6 ) > 1]
GI <- unique(gtf.anno$V6)
for (n in c("SPH2016","D3post","Blakeley","nBGuo","nicolBla" )) {
  counts[[n]]$Gene <- gtf.anno[counts[[n]]$Gene,]$V6
  counts[[n]] <- (counts[[n]] %>% filter(Gene %in% dup.gene ) %>% gather(cell,counts,-Gene) %>% group_by(Gene,cell) %>% summarise(counts=sum(counts)) %>% spread(cell,counts) %>% select(colnames(counts[[n]])) %>% bind_rows(counts[[n]] %>% filter(!Gene %in% dup.gene ) ) %>% tibble::column_to_rownames("Gene"))[GI,]
  metas[[n]] <- metas[[n]] %>% left_join(data.frame(cell=colnames(counts[[n]]),mt.perc=colSums(counts[[n]][mt.gene,])/colSums(counts[[n]])) %>% tbl_df() ,by="cell")
}
for (n in c("JPF2019","IBD2","EBD2" )) {
  counts[[n]] <- (counts[[n]] %>% filter(Gene %in% c(dup.gene,paste0(dup.gene,".1")))  %>% mutate(Gene=sub("\\.1","",Gene)) %>% gather(cell,counts,-Gene)  %>% group_by(Gene,cell) %>% summarise(counts=sum(counts)) %>% spread(cell,counts) %>% select(colnames(counts[[n]])) %>% bind_rows(counts[[n]] %>%  filter(!Gene %in% c(dup.gene,paste0(dup.gene,".1"))))  %>% tibble::column_to_rownames("Gene"))[GI,]
  metas[[n]] <- metas[[n]] %>% left_join(data.frame(cell=colnames(counts[[n]]),mt.perc=colSums(counts[[n]][mt.gene,])/colSums(counts[[n]])) %>% tbl_df() ,by="cell")
}
for (n in c("CS7" )) {
  counts[[n]]$Gene <- gtf.anno[counts[[n]]$Gene,]$V6
  counts[[n]] <- (counts[[n]] %>% filter(Gene %in% dup.gene ) %>% gather(cell,counts,-Gene) %>% group_by(Gene,cell) %>% summarise(counts=sum(counts)) %>% spread(cell,counts) %>% select(colnames(counts[[n]])) %>% bind_rows(counts[[n]] %>% filter(!Gene %in% dup.gene ) ) %>% tibble::column_to_rownames("Gene"))[GI,]
  metas[[n]] <- metas[[n]] %>% left_join(data.frame(cell=colnames(counts[[n]]),mt.perc=colSums(counts[[n]][mt.gene,])/colSums(counts[[n]])) %>% tbl_df() ,by="cell")
}


counts.all <-do.call("bind_cols",counts)
meta.all <- do.call("bind_rows",metas)


#' saving files
TD="Nov14_2021"
if (!file.exists(paste0("tmp_data/",TD,"/all.counts.meta.Rdata")) | rewrite) {
  print("save output")
  save(counts.all,meta.all,file=paste0("tmp_data/",TD,"/all.counts.meta.Rdata"))
  save(gtf.anno,mt.gene,ribo.gene,file="tmp_data/gene.meta.Rdata")
}

