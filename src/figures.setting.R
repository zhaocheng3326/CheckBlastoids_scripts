#' figure setting

zs.limit <- 2.5
heat.col <- colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100)
#heat2.col <- colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#CC4678FF","#F89441FF","#F89441FF","#F0F921FF","#F0F921FF","#F0F921FF","#F0F921FF"))(100)
rename <- dplyr::rename
select<- dplyr::select
filter <- dplyr::filter
#' setting 
lineage.col.set <-c("Amnion"="#E41A1C" ,"AMLC"="#E41A1C" , "ELC"="#4DAF4A" , "Epiblast"="#4DAF4A","Endoderm"="#984EA3","HLC"="#984EA3","TE"= "#377EB8", "TLC"= "#51C0CC","MeLC"= "#FCCDE5","Mesoderm"= "#FCCDE5", "Prelineage"= "#B35806","ICM"="#FFFFB3","ICM-TE_trans"="#FDB863", "PriS"="#DB8E00",  "Undef"= "grey50", "Unknown"="grey33","ExE_Mes"="#F781BF","hES"="#00FF92","hPGCLC"="#CD8BBC","Stem-blastoids"="#F8766D","Iblastoids"="#00BA38","Embryonic"="royalblue3","Stem-blastoids"="grey33","Iblastoids"="grey33","EarlyBlastocyst"="#FDB863","8C_Morula"="#B35806","Mes"="#FCCDE5","EPI_Amnion"="#FDB863","naive_H9"="#4DFA4A","primed_H9"="#00FF92","okae_bts5"="black")

lineage.bind.col.set <-c("Amnion"="#E41A1C" , "ELC"="#CCEBC5" , "Epiblast"="#4DAF4A","Endoderm"="#984EA3","HLC"="#A66DCC","TE"= "#377EB8", "TLC"= "#51C0CC","Mesoderm"= "#FCCDE5",  "PriS"="#DB8E00",  "PriS_Amnion"="#B35806","Mes"="#FCCDE5","EPI_Amnion"="#4EFF0D","EPI_PriS"="#7C2136","PriS_Amnion"="#FCF49B")

cross.lineage.col.set <-c("Amnion(human)"="#E41A1C" ,"Amnion(NHP)"="#FC000D" , "E-Amnion"="#FED9A6","ELC"="#CCEBC5" , "Epiblast"="#4DAF4A","Endoderm"="#984EA3","PE"="#984EA3","HLC"="#A66DCC","TE"= "#377EB8", "TLC"= "#51C0CC","Mesoderm"= "#FCCDE5", "EPI-AM-trans"= "yellow","ICM"="#FFFFB3","ICM-TE_trans"="#FDB863", "PriS"="#DB8E00",  "Undef"= "grey50","Unknown"="grey33", "ExE_Mes"="#F781BF","ExE_Mech"="#B35806","naive_H9"="#4DFA4A","primed_H9"="#00FF92","okae_bts5"="black")

#' setting 
recheck.lineage.col.set <-c("Amnion"="#E41A1C" ,"Epi"="#4DAF4A","Epi(Prev-PSA-EPI)"="#4DAF4A","Endoderm"="#984EA3","Endoderm(Prev-PSA-EPI)"="#984EA3","TE"= "#377EB8", "TE(Prev-ICM)"= "#377EB8","TE(Prev-PSA-EPI)"= "#377EB8","Mesoderm"= "#FCCDE5", "Mes(Prev-Epi)"="#FCCDE5", "Mes(Prev-Epi)"="#FCCDE5","Mes(Prev-PSA-EPI)"="#FCCDE5", "Epiblast"="#4DAF4A","PriS"="#DB8E00","EPI_Amnion"="#FDB863","EPI_PriS"="#7C2136","PriS_Amnion"="#FCF49B")

species.shape.set <- c(human=16,NHP=17)


heat.col <- colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100)

vln.col.set <- c("TE"="#DA8F00","CTB"= "#86AC00","EVT"="#00C094","STB"= "#00B6EB" ,"Amnion"="#F8766D","TLC(SB)"="#DA8F00","TLC(IB)"="#FF6B96","AMLC"= "#F8766D","Tsw-AMLC"= "#F8766D","ELC"="#00BADE" ,"TLC"= "#64B200","HLC"="#B385FF","MeLC"= "#EF67EB")

#cluster.col.set <- c("C0"="#00C094","C1"="#E38900","C2"="#99A800","C3"="#00BFC4","C4"="#00BC56","C5"="#F8766D","C6"="#53B400","C7"="#C49A00","C8"="#06A4FF","C9"="#A58AFF","C10"="#DF70F8","C11"="#FB61D7","C12"="#FF66A8","Others"="grey66")

EML.shape.set <- c('8C_Morula'=0,'EarlyBlastocyst'=1,'Epiblast'=2,'Endoderm'=3,'PriS'=4,'Mes'=5,'ExE_Mes'=6,'Amnion'=17,'TE'=15)



pj.shape.set <- c(CS7=21,D3post=22,SPH2016=24,nBGuo=23,Blakeley=9, zhou2019=13)
pj.shape.set.solid <- c(CS7=16,D3post=15,SPH2016=17,nBGuo=18)
label.pj <- list()
label.pj$SPH2016="Petropoulos et al., 2016"
label.pj$D3post="Xiang et al., 2020"
label.pj$CS7="Tyser et al., 2021"
label.pj$JPF2019="Zheng et al., 2019"
label.pj$nBGuo="Yanagida et al., 2021"

pj.col <- c("SPH2016"="#FF61C6","D3post"="#7CA800","D3"="#7CA800","D3"="#7CA800","CS7"="#C79000","Blakeley"="#F2766D","zhou2019"="lightblue","IBD2"="#00B9BE","EBD2"="#00B867","nBGuo"="#C17CFF","JPF2019"="#00A3FF",NHP="purple",nicolBla="royalblue")




#' 
#' heat.col <- colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100)
#' rename <- dplyr::rename
#' select<- dplyr::select
#' filter <- dplyr::filter
#' #' setting 
#' lineage.col.set <-c("Amnion"="#F8766D" ,"AMLC"="#F8766D" , "ELC"="#00BADE" , "Epiblast"="#00BADE","Endoderm"="#B385FF","HLC"="#B385FF","TE"= "#64B200", "TLC"= "#64B200","MeLC"= "#EF67EB","Mesoderm"= "#EF67EB", "Prelineage"= "#7B554E", "PriS"="#DB8E00",  "Undef"= "grey50", "ExE_Mes"="#C4423E","hES"="#51C0CC","hPGCLC"="#CD8BBC","Stem-blastoids"="#F8766D","Iblastoids"="#00BA38","Embryonic"="royalblue3","ICM"="grey50","Stem-blastoids"="grey33","Iblastoids"="grey33")
#' 
#' heat.col <- colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100)
#' 
#' vln.col.set <- c("TE"="#DA8F00","CTB"= "#86AC00","EVT"="#00C094","STB"= "#00B6EB" ,"Amnion"="#F8766D","TLC(SB)"="#DA8F00","TLC(IB)"="#FF6B96","AMLC"= "#F8766D","Tsw-AMLC"= "#F8766D","ELC"="#00BADE" ,"TLC"= "#64B200","HLC"="#B385FF","MeLC"= "#EF67EB")
#' 
#' cluster.col.set <- c("C0"="#00C094","C1"="#E38900","C2"="#99A800","C3"="#00BFC4","C4"="#00BC56","C5"="#F8766D","C6"="#53B400","C7"="#C49A00","C8"="#06A4FF","C9"="#A58AFF","C10"="#DF70F8","C11"="#FB61D7","C12"="#FF66A8","Others"="grey66")
#' 
#' pj.col.set <- c(CS7="#F8766D",D3post="#B79F00",EBD2="#00BA38",IBD2="#00BFC4",JPF2019="#619CFF",SPH2016="#F564E3")
#' pj.shape.set <- c(CS7=21,D3post=22,SPH2016=24)
#' pj.shape.set.solid <- c(CS7=16,D3post=15,SPH2016=17)
#' label.pj <- list()
#' label.pj$SPH2016="Petropoulos et al., 2016"
#' label.pj$D3post="Xiang et al., 2020"
#' label.pj$CS7="Tyser et al., 2021"
#' label.pj$JPF2019="Zheng et al., 2019"
