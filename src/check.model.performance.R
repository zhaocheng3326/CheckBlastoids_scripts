#' ---
#' title: "model prediction"
#' output: 
#'  html_document:
#'    code_folding: hide
#'    
#' ---
#' 
#' 
# R4.0

rm(list=ls())
rewrite=FALSE

condaENV <- "/home/chenzh/miniconda3/envs/R4.0"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)


suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(Seurat))
suppressMessages(library(ggalluvial))
suppressMessages(library(readxl))
#source("~/PC/SnkM/SgCell.R")
rename <- dplyr::rename
#argv=commandArgs(TRUE)

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


pv1 <- 0.015 
pv2 <- 0.05

# selected by 
#query_ref %>% filter(ref_EML=="Amnion") %>% mutate(query_EML=FunForMatAnno(query_EML)) %>% filter(query_EML=="TE") %>% group_by(query_cell) %>% top_n(1,-1*pval)%>% arrange(pval)  %>% head(300) %>% tail(1)


# loading meta data
ref.pj <-  c("D3post","CS7","SPH2016")
meta.filter<- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% mutate(EML=ifelse(EML %in% c("Ectoderm"),"Amnion",EML))

meta.filter <- meta.filter

ref.meta <- readRDS(paste0("tmp_data/",TD,"/ref.meta.rds")) %>% select(ref_cell,ref_EML,ref_pj) #%>% mutate(ref_EML=ifelse(ref_EML=="EPI_PriS","Epi",ref_EML))
ref.meta.detail <- ref.meta %>% left_join(meta.filter %>% select(cell,devTime) %>% rename(ref_cell=cell,ref_devTime=devTime),by="ref_cell") %>% mutate(ref_devTime=ifelse(ref_devTime %in% c("E3","E4","E5","E6","E7"),"Pre",ifelse(ref_devTime %in% c("E8","E9","E10","E11","E12","E13","E14"), "Post","CS7"))) %>% mutate(ref_EML.detail=paste(ref_EML,ref_devTime,sep="_")) %>% select(-ref_devTime)

query.meta <- meta.filter %>% filter(!cell %in% ref.meta$ref_cell) %>% rename(query_cell=cell,query_EML=EML,query_pj=pj) %>% bind_rows(ref.meta  %>% rename(query_cell=ref_cell,query_EML=ref_EML,query_pj=ref_pj)) %>% select(query_cell,query_EML,query_pj)

#' loading NHP meta
load(paste0("tmp_data/",TD,"/NHP.Rdata"),verbose=T)
query.meta <- query.meta %>% bind_rows(NHP.meta  %>% rename(query_cell=cell,query_EML=EML,query_pj=pj)  %>% select(query_cell,query_EML,query_pj) )
rm(GI.list);rm(NHP.counts.mod);rm(NHP.meta)

#' meta data from Castel et al., 2020 meta data
CastelBulk.meta <- read.delim("data/Castel_et_al_2020-master/Results/sampleAnnot.mergedTechReplicate.tsv",stringsAsFactors = F,head=T) %>% tbl_df() %>% rename(cell=Name,EML=Group) %>% select(cell,EML) %>% mutate(pj="CastelBulk")
query.meta <- query.meta %>% bind_rows(CastelBulk.meta  %>% rename(query_cell=cell,query_EML=EML,query_pj=pj)  %>% select(query_cell,query_EML,query_pj) )

#' meta data from zhou et al.,
zhou2019.meta <- read.delim("data/Castel_et_al_2020-master/Dataset/ZhouPetro/sampleAnnot.tsv",,stringsAsFactors = F,head=T) %>% tbl_df() %>% filter(Dataset=="Zhou") %>% rename(cell=Name,EML=finalClusters) %>% mutate(batch="batch_zhou2019",pj="zhou2019")
query.meta <- query.meta %>% bind_rows(zhou2019.meta  %>% rename(query_cell=cell,query_EML=EML,query_pj=pj)  %>% select(query_cell,query_EML,query_pj) )




#' decide to  use 2500
PDVG=2000;REFVG="2500"
pred_dir=paste0("tmp_data/",TD,"/predict_dir/output/PDVG",PDVG,"_REFVG",REFVG)
pred_files <- list.files(pred_dir,pattern="_pred.out")
query_ref.raw <- lapply( pred_files,function(x) {read.delim(paste(pred_dir,x,sep="/")) %>% tbl_df()}) %>% do.call("bind_rows",.) %>% left_join(ref.meta,by="ref_cell") %>% left_join(query.meta,by="query_cell") %>% filter(pval <pv2)  
query_ref <- query_ref.raw %>% filter(pval <pv2)  

#' ## for the most significant ones
query_ref %>% filter(ref_EML=="TE" | ref_EML=="Amnion") %>% group_by(query_cell,ref_EML) %>% summarise(pval=min(pval)) %>% ungroup() %>% spread(ref_EML,pval) %>% replace(.,is.na(.),0.1) %>% inner_join(query.meta,by="query_cell")  %>% ggplot+geom_point(mapping=aes(x=-log10(Amnion),y=-log10(TE),col=query_pj))+geom_vline(xintercept=-log10(0.015))+geom_hline(yintercept=-log10(0.015))+xlim(0.9,4)+ylim(0.9,4)+facet_wrap(~query_pj)

query_ref.sig <- query_ref %>% filter(pval <pv1)
query_ref.sig.int.type <- query_ref.sig  %>% filter(ref_EML %in% c("EPI_Amnion","EPI.PrE.INT","PriS_Amnion")) %>% select(query_cell,ref_EML,ref_pj) %>% unique()

query_ref.sig.int.type <- query_ref.sig.int.type  %>% mutate(ref_EML=ifelse(ref_EML=="EPI_Amnion","Epi",ref_EML)) %>% mutate(ref_EML=ifelse(ref_EML=="PriS_Amnion","PriS",ref_EML))  %>% mutate(ref_EML=ifelse(ref_EML=="EPI.PrE.INT","Epi",ref_EML)) %>% bind_rows(query_ref.sig.int.type %>% mutate(ref_EML=ifelse(ref_EML=="EPI_Amnion","Amnion",ref_EML)) %>% mutate(ref_EML=ifelse(ref_EML=="PriS_Amnion","Amnion",ref_EML))  %>% mutate(ref_EML=ifelse(ref_EML=="EPI.PrE.INT","Endoderm",ref_EML)))

pSig_frac.cutoff <- 0.5 
pSig_nc_hits.cutoff <- 2


# prediction results for each ref dataset
query_ref.sig.anno  <- query_ref.sig %>% group_by(query_cell,query_EML,query_pj,ref_pj,ref_EML) %>%  summarise(nc_hits=n_distinct(ref_cell),mean_pval=mean(pval),min_pval=min(pval)) %>% inner_join(query_ref.sig %>% group_by(query_cell,query_EML,query_pj,ref_pj) %>% summarise(n_hits=n_distinct(ref_cell)),by=c("query_cell","query_EML","query_pj","ref_pj"))%>% ungroup() %>% mutate(majority_frac=nc_hits/n_hits) %>% FunHandleINT()%>% filter(nc_hits>= pSig_nc_hits.cutoff & majority_frac > pSig_frac.cutoff )
# for all reference dataset,merge all results  selected by the pvalues
query_ref.sig.anno.merge <- query_ref.sig.anno %>% group_by(query_cell,query_EML,query_pj)%>% top_n(1,-1*mean_pval)%>% top_n(1,majority_frac) %>% top_n(1,nc_hits)  %>% ungroup() %>% mutate(Type=paste0("Pval<",pv1))     

#' for the  predictions less than pv2 (but pv1 supported)
query_ref.loose.part1 <- query_ref %>% filter((!query_cell %in% c(query_ref.sig.anno.merge$query_cell)) & query_cell %in% query_ref.sig$query_cell)%>% filter(pval < pv2)

#' for the predictions only less than pv2
query_ref.loose.part2 <- query_ref %>% filter((!query_cell %in% query_ref.sig$query_cell))%>% filter(pval < pv2)

pL1_nc_hits.cutoff <- 2
pL1_frac.cutoff <- 3/5-0.001
pL2_nc_hits.cutoff <- 2
pL2_frac.cutoff <- 4/5-0.001

query_ref.loose.part1.anno <- query_ref.loose.part1 %>% group_by(query_cell,query_EML,query_pj,ref_pj,ref_EML) %>%  summarise(nc_hits=n_distinct(ref_cell),mean_pval=mean(pval),min_pval=mean(pval)) %>% inner_join(query_ref.loose.part1  %>% group_by(query_cell,query_EML,query_pj,ref_pj) %>% summarise(n_hits=n_distinct(ref_cell)),by=c("query_cell","query_EML","query_pj","ref_pj"))%>% ungroup() %>% mutate(majority_frac=nc_hits/n_hits) %>% filter(majority_frac > pL1_frac.cutoff  & nc_hits >= pL1_nc_hits.cutoff )
#' only consider the ones supported by at least one case (pvalue < pv1)
query_ref.loose.part1.anno.support <- query_ref.loose.part1.anno  %>% semi_join( query_ref.sig %>% select(query_cell,ref_EML) %>% unique(),by=c("query_cell","ref_EML")) %>% unique()


query_ref.loose.part1.anno.merge <- query_ref.loose.part1.anno.support %>% group_by(query_cell,query_EML,query_pj) %>% top_n(1,-1*mean_pval) %>% top_n(1,majority_frac)  %>% top_n(1,nc_hits)  %>% ungroup() %>% mutate(Type=paste0("Pval<",pv1, " & Pval <",pv2))       
query_ref.loose.part2.anno <- query_ref.loose.part2 %>% group_by(query_cell,query_EML,query_pj,ref_pj,ref_EML) %>%  summarise(nc_hits=n_distinct(ref_cell),mean_pval=mean(pval),min_pval=min(pval)) %>% inner_join(query_ref.loose.part2 %>% group_by(query_cell,query_EML,query_pj,ref_pj) %>% summarise(n_hits=n_distinct(ref_cell)),by=c("query_cell","query_EML","query_pj","ref_pj"))%>% ungroup() %>% mutate(majority_frac=nc_hits/n_hits) %>% FunHandleINT()  %>% filter(majority_frac > pL2_frac.cutoff  & nc_hits >= pL2_nc_hits.cutoff )# %>% FunHandleINT() 


temp.query_ref.loose.part2.sel.query_cell <-  c(query_ref.loose.part2.anno  %>% select(query_cell,ref_EML) %>% unique()  %>% group_by(query_cell) %>% summarise(nC=n()) %>% filter(nC ==1) %>% pull(query_cell) )
query_ref.loose.part2.anno.merge <- query_ref.loose.part2.anno %>% filter(query_cell %in% temp.query_ref.loose.part2.sel.query_cell ) %>% group_by(query_cell,query_EML,query_pj) %>% top_n(1,-1*mean_pval) %>% top_n(1,majority_frac) %>% ungroup()%>% mutate(Type=paste0("Pval<",pv1, " & Pval <",pv2))  


query_ref.anno.out <- query_ref.sig.anno.merge %>% bind_rows(query_ref.loose.part1.anno.merge) %>% bind_rows(query_ref.loose.part2.anno.merge)
query_ref.anno.out <- query_ref.anno.out  %>% bind_rows(query_ref %>% select(query_cell, query_EML,query_pj) %>% filter(!query_cell %in% query_ref.anno.out$query_cell)  %>% filter(query_cell %in% query_ref$query_cell ) %>% mutate(ref_EML="ambiguous",ref_pj="ambiguous")%>% unique()) %>% arrange(query_cell, query_EML,query_pj)
query_ref.anno.out  <- query_ref.anno.out  %>% bind_rows(query.meta %>% select(query_cell, query_EML,query_pj) %>% filter(!query_cell %in% query_ref.anno.out$query_cell) %>% mutate(ref_EML="NoSigHits",ref_pj="NoSigHits")%>% unique()) %>% arrange(query_cell, query_EML,query_pj)

#+ fig.width=9,fig.height=9
#pdf("temp.pdf",7,7)

for (p in unique(query_ref.anno.out$query_pj)) {
  temp.stat <- query_ref.anno.out %>% filter(query_pj==p) %>% group_by(ref_EML,query_EML) %>% summarise(nCell=n_distinct(query_cell)) %>% spread(ref_EML,nCell)  %>% replace(.,is.na(.),0) %>% tibble::column_to_rownames("query_EML") %>% apply(1,function(x){round(x)}) %>% t()
  pheatmap::pheatmap(temp.stat,scale="row",cluster_rows=F,cluster_cols=F,display_numbers=temp.stat,main=paste("Prediction of",p),legend = FALSE)
  
  temp.query <- query_ref.anno.out %>% filter(query_pj==p) %>% filter(!ref_EML %in% c("ambiguous","NoSigHits")) %>% mutate(query_EML=FunForMatAnno(query_EML),ref_EML=FunForMatAnno(ref_EML)) %>% group_by(ref_EML,query_EML) %>% summarise(nCell=n_distinct(query_cell))  %>% to_lodes_form(axes=1:2,id="Cohort") %>% mutate(x=factor(x,c("query_EML","ref_EML"),ordered = T)) %>% mutate(stratum=as.vector(stratum)) %>% mutate(stratum=as.factor(stratum))
  print(
    ggplot(temp.query,aes(x = x, stratum = stratum, alluvium = Cohort, y = nCell,fill = stratum, label = stratum)) + scale_x_discrete(expand = c(.1, .1)) +geom_flow() +geom_stratum(alpha = .5) +geom_text(stat = "stratum", size = 3) +theme(legend.position = "none") +ggtitle("")+theme_void()+NoLegend()+ggtitle(p)+theme(plot.title = element_text(hjust=0.5))
  )
}
#dev.off()

if (!file.exists(paste0("tmp_data/",TD,"/prediction.model.performance.rds")) | rewrite) {
  print("save output")
  saveRDS(query_ref.anno.out,file= (paste0("tmp_data/",TD,"/prediction.model.performance.rds")))
  saveRDS(query_ref.raw,file= (paste0("tmp_data/",TD,"/query_ref.raw.rds")))
}



print(
  query_ref.anno.out  %>%  select(query_cell,ref_EML)  %>% mutate(ref_EML=ifelse(ref_EML %in% c("ambiguous","NoSigHits"),ref_EML,"Annotated" )) %>% rename(cell=query_cell,Type=ref_EML) %>% inner_join(meta.filter,by="cell") %>% ggplot+geom_histogram(mapping=aes(x=nGene,fill=Type))+geom_vline(xintercept=2000)+ facet_wrap(~Type)
)
print(
  query_ref.anno.out %>% filter(query_pj=="EBD2")  %>%  select(query_cell,ref_EML)  %>% mutate(ref_EML=ifelse(ref_EML %in% c("ambiguous","NoSigHits"),ref_EML,"Annotated" )) %>% rename(cell=query_cell,Type=ref_EML) %>% inner_join(meta.filter,by="cell") %>% ggplot+geom_histogram(mapping=aes(x=nGene,fill=Type))+geom_vline(xintercept=2000)+ facet_wrap(~Type)+ggtitle("EBD2")+FunTitle()
)
print(
  query_ref.anno.out %>% filter(query_pj=="IBD2")  %>%  select(query_cell,ref_EML)  %>% mutate(ref_EML=ifelse(ref_EML %in% c("ambiguous","NoSigHits"),ref_EML,"Annotated" )) %>% rename(cell=query_cell,Type=ref_EML) %>% inner_join(meta.filter,by="cell") %>% ggplot+geom_histogram(mapping=aes(x=nGene,fill=Type))+geom_vline(xintercept=2000)+ facet_wrap(~Type)+ggtitle("IBD2")+FunTitle()
)



