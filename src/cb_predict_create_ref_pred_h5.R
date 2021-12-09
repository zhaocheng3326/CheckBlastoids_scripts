#' title: "create h5 files and training model"
# R4.0
rm(list=ls())
re.run <- FALSE

condaENV <- "/home/chenzh/miniconda3/envs/R4.0"
LBpath <- paste0(condaENV ,"/lib/R/library")
.libPaths(LBpath)

suppressMessages(library(dplyr))
suppressMessages(library(tidyr))
suppressMessages(library(ggplot2))
suppressMessages(library(scran))
suppressMessages(library(batchelor))
suppressMessages(library(Seurat))
suppressMessages(library(harmony))
#suppressMessages(library(scPred))
suppressMessages(library(SeuratDisk))
suppressMessages(library(SeuratData))
suppressMessages(library(rhdf5))
suppressMessages(library(SeuratWrappers))
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

#' loading data
load("tmp_data/gene.meta.Rdata",verbose=T)
counts.filter <- readRDS(paste0("tmp_data/",TD,"/counts.filter.rds"))

meta.filter<- readRDS(paste0("tmp_data/",TD,"/meta.filter.rds")) %>% mutate(EML=ifelse(EML %in% c("Ectoderm"),"Amnion",EML)) 


temp.M <- meta.filter %>% filter( pj %in% c("D3post","CS7","SPH2016"))
#' get the expressed genes 
expG.set <- list()
for (b in c("D3post","CS7","SPH2016")) { 
  temp.cell <- temp.M%>% filter(pj==b) %>% pull(cell)
  expG.set[[b]] <- rownames(counts.filter )[rowSums(counts.filter[,temp.cell] >=1) >=5]
}
temp.sel.expG <-unlist(expG.set) %>% unique() %>% as.vector()

data.merge <- CreateSeuratObject(counts.filter[temp.sel.expG,c(temp.M$cell)], meta.data = (temp.M %>% tibble::column_to_rownames("cell"))) %>% NormalizeData(verbose = FALSE) 
data.spt <- SplitObject(data.merge, split.by = "pj") 


if (re.run) {
  for (nGene in c(2500)) {
    dir.create(paste0("tmp_data/",TD,"/cb_predict_ref_VG",nGene))
    data.spt.temp <- data.spt %>% lapply(function(x){x=FindVariableFeatures(x,nfeatures=nGene,verbose=F)})
    for ( n in names(data.spt)) {
      SaveH5Seurat(data.spt.temp[[n]], filename = paste0("tmp_data/",TD,"/cb_predict_ref_VG",nGene,"/ref_",n,".h5Seurat"))
      Convert(paste0("tmp_data/",TD,"/cb_predict_ref_VG",nGene,"/ref_",n,".h5Seurat"), dest = "h5ad")
    }
  }
  
  system("python src/cb_predict_ref_create.py ~/My_project/CheckBlastoids/tmp_data/Nov14_2021/cb_predict_ref_VG2500/")
  #system("python src/cb_predict_ref_create.py ~/My_project/CheckBlastoids/tmp_data/Nov14_2021/cb_predict_ref_VG1000/")
  #system("python src/cb_predict_ref_create.py ~/My_project/CheckBlastoids/tmp_data/Nov14_2021/cb_predict_ref_VG1500/")
  #system("python src/cb_predict_ref_create.py ~/My_project/CheckBlastoids/tmp_data/Nov14_2021/cb_predict_ref_VG2000/")
  #
  #system("python src/cb_predict_ref_create.py ~/My_project/CheckBlastoids/tmp_data/Nov14_2021/cb_predict_ref_VG3000/")
  #system("python src/cb_predict_ref_create.py ~/My_project/CheckBlastoids/tmp_data/Nov14_2021/cb_predict_ref_VG3500/")
  #system("python src/cb_predict_ref_create.py ~/My_project/CheckBlastoids/tmp_data/Nov14_2021/cb_predict_ref_VG4000/")
  #system("python src/cb_predict_ref_create.py ~/My_project/CheckBlastoids/tmp_data/Nov14_2021/cb_predict_ref_VG4500/")
  #system("python src/cb_predict_ref_create.py ~/My_project/CheckBlastoids/tmp_data/Nov14_2021/cb_predict_ref_VG5000/")
  #system("python src/cb_predict_ref_create.py ~/My_project/CheckBlastoids/tmp_data/Nov14_2021/cb_predict_ref_VG5500/")
  #system("python src/cb_predict_ref_create.py ~/My_project/CheckBlastoids/tmp_data/Nov14_2021/cb_predict_ref_VG6000/")
  #system("python src/cb_predict_ref_create.py ~/My_project/CheckBlastoids/tmp_data/Nov14_2021/cb_predict_ref_VGSPG2500/")
  #system("python src/cb_predict_ref_create.py ~/My_project/CheckBlastoids/tmp_data/Nov14_2021/cb_predict_ref_VGSPG3000/")
  
}


#' output predicted data
dir.create(paste0("tmp_data/",TD,"/predict_dir"))

# reprocessed human dataset
for (b in unique(meta.filter$pj)) {
  print(b)
  temp.counts <- counts.filter[,meta.filter %>% filter(pj==b) %>% pull(cell)]
  gz1 <- gzfile(paste0("tmp_data/",TD,"/predict_dir/",b,"_cell.counts.gz"), "w")
  write.csv(temp.counts  %>% as.matrix() %>% as.data.frame()%>% tibble::rownames_to_column("Gene"), gz1,quote=F,row.names=F)
  close(gz1)
}
# NHP dataset
load(paste0("tmp_data/",TD,"/NHP.Rdata"),verbose=T)
gz1 <- gzfile(paste0("tmp_data/",TD,"/predict_dir/","NHP","_cell.counts.gz"), "w")
write.csv(NHP.counts.mod %>% tibble::rownames_to_column("Gene"), gz1,quote=F,row.names=F)
close(gz1)

#bulk dataset for Castel_et_al_2020
temp.counts <- read.delim("data/Castel_et_al_2020-master/Results/exprRaw.mergedTechReplicate.tsv",stringsAsFactors = F,head=T,row.names = 1)
gz1 <- gzfile(paste0("tmp_data/",TD,"/predict_dir/","CastelBulk","_cell.counts.gz"), "w")
write.csv(temp.counts%>% tibble::rownames_to_column("Gene"), gz1,quote=F,row.names=F)
close(gz1)

#' zhou et al 2019  https://www.nature.com/articles/s41586-019-1500-0
zhou2019.meta <- read.delim("data/Castel_et_al_2020-master/Dataset/ZhouPetro/sampleAnnot.tsv",,stringsAsFactors = F,head=T) %>% tbl_df() %>% filter(Dataset=="Zhou") %>% rename(cell=Name,EML=finalClusters) %>% mutate(batch="batch_zhou2019",pj="zhou2019")
temp.counts <- read.delim("data/Castel_et_al_2020-master/Dataset/Zhou2019/exprDatRawZhou2019.tsv",stringsAsFactors = F,head=T,row.names = 1)[,zhou2019.meta$cell]
gz1 <- gzfile(paste0("tmp_data/",TD,"/predict_dir/","zhou2019","_cell.counts.gz"), "w")
write.csv(temp.counts%>% tibble::rownames_to_column("Gene"), gz1,quote=F,row.names=F)
close(gz1)






#' generate pred h5ad file
for (b in c("Blakeley","CS7","D3post","EBD2","IBD2","JPF2019","nBGuo","NHP","SPH2016","CastelBulk","zhou2019",,"nicolBla")) {
  temp.counts <- read.delim(gzfile(paste0("tmp_data/",TD,"/predict_dir/",b,"_cell.counts.gz")),header = T,row.names = 1,sep=",")
  for (nVG in c(2000)) {
    print(paste(b,nVG))
    data.temp <- CreateSeuratObject(temp.counts) %>% NormalizeData(verbose = FALSE) %>% FindVariableFeatures(nfeatures=nVG,verbose=F) 
    SaveH5Seurat(data.temp , filename = paste0("tmp_data/",TD,"/predict_dir/",b,"_",nVG,"_cell.counts.gz.h5Seurat"),verbose=F,overwrite=T)
    Convert(paste0("tmp_data/",TD,"/predict_dir/",b,"_",nVG,"_cell.counts.gz.h5Seurat"), dest = "h5ad",verbose=F,overwrite=T)
  }
}


re.run <- FALSE
if (re.run) {
  system(
    '
 for REFVG in 2500 #2000  3000 3500 4000 4500 5000 5500 SPG2500 SPG3000 #6000
  do
	for PDVG in 2000 
	do
		for pj in Blakeley CS7 D3post EBD2 IBD2 JPF2019 nBGuo NHP SPH2016 CastelBulk zhou2019 nicolBla
		do 
			echo "$REFVG $PDVG $pj"
			mkdir -p ~/My_project/CheckBlastoids/tmp_data/Nov14_2021/predict_dir/output/PDVG${PDVG}_REFVG${REFVG}
			python ~/My_project/CheckBlastoids/src/cb.calculator.py -i ~/My_project/CheckBlastoids/tmp_data/Nov14_2021/predict_dir/${pj}_${PDVG}_cell.counts.gz.h5ad -c ~/My_project/CheckBlastoids/tmp_data/Nov14_2021/cb_predict_ref_VG${REFVG}/ -o ~/My_project/CheckBlastoids/tmp_data/Nov14_2021/predict_dir/output/PDVG${PDVG}_REFVG${REFVG}/${pj}_pred.out
		done
	done
done
    '
  )
}

