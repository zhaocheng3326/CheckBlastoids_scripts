#' figure setting

heat.col <- colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100)
rename <- dplyr::rename
select<- dplyr::select
filter <- dplyr::filter
#' setting 
lineage.col.set <-c("Amnion"="#F8766D" ,"AMLC"="#F8766D" , "ELC"="#00BADE" , "Epiblast"="#00BADE","Endoderm"="#B385FF","HLC"="#B385FF","TE"= "#64B200", "TLC"= "#64B200","MeLC"= "#EF67EB","Mesoderm"= "#EF67EB", "Prelineage"= "#7B554E", "PriS"="#DB8E00",  "Undef"= "grey50", "ExE_Mes"="#C4423E","hES"="#51C0CC","hPGCLC"="#CD8BBC","Stem-blastoids"="#F8766D","Iblastoids"="#00BA38","Embryonic"="royalblue3","ICM"="grey50","Stem-blastoids"="grey33","Iblastoids"="grey33")

heat.col <- colorRampPalette(c("#0D0887FF","#0D0887FF","#0D0887FF","#0D0887FF","#7E03A8FF","#CC4678FF","#F89441FF","#F0F921FF","#F0F921FF"))(100)

vln.col.set <- c("TE"="#DA8F00","CTB"= "#86AC00","EVT"="#00C094","STB"= "#00B6EB" ,"Amnion"="#F8766D","TLC(SB)"="#DA8F00","TLC(IB)"="#FF6B96","AMLC"= "#F8766D","Tsw-AMLC"= "#F8766D","ELC"="#00BADE" ,"TLC"= "#64B200","HLC"="#B385FF","MeLC"= "#EF67EB")

cluster.col.set <- c("C0"="#00C094","C1"="#E38900","C2"="#99A800","C3"="#00BFC4","C4"="#00BC56","C5"="#F8766D","C6"="#53B400","C7"="#C49A00","C8"="#06A4FF","C9"="#A58AFF","C10"="#DF70F8","C11"="#FB61D7","C12"="#FF66A8","Others"="grey66")

pj.col.set <- c(CS7="#F8766D",D3post="#B79F00",EBD2="#00BA38",IBD2="#00BFC4",JPF2019="#619CFF",SPH2016="#F564E3")
pj.shape.set <- c(CS7=21,D3post=22,SPH2016=24)
pj.shape.set.solid <- c(CS7=16,D3post=15,SPH2016=17)
label.pj <- list()
label.pj$SPH2016="Petropoulos et al., 2016"
label.pj$D3post="Xiang et al., 2020"
label.pj$CS7="Tyser et al., 2021"
label.pj$JPF2019="Zheng et al., 2019"
