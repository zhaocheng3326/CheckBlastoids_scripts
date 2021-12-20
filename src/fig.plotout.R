#' R4.0
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
suppressMessages(library(ggalluvial))

# working directory
DIR <- "~/My_project/CheckBlastoids"
setwd(DIR)

#' Loading R functions
#source("~/PC/R_code/functions.R")
#source("~/PC/SnkM/SgCell.R")
source("src/local.quick.fun.R")
source("src/figures.setting.R")


#' loading options
TD="Nov14_2021"
MM="Pub"

load(paste0("tmp_data/",TD,"/figplot.Rdata"),verbose=T)



pdf("temp.fig1A.pdf",15,18)
(
  plot_spacer()+plot.results$umap.em+plot.results$umap.em.legend+plot_spacer()+plot_layout(widths = c(0.6,2,1,0.6))

)/(
  ( plot.results$umap.IBD2|plot.results$umap.EBD2|plot.results$umap.nicolBla1|plot.results$umap.nBGuo)
  
)/(
    plot.results$NHP.IBD2.mnn|plot.results$NHP.EBD2.mnn|plot.results$NHP.nicolBla.mnn|plot.results$NHP.nBGuo.mnn
  ) +plot_layout(heights = c(2.5,1.75,1.75))

(
  plot_spacer()+plot.results.jpeg$umap.em+plot.results$umap.em.legend+plot_spacer()+plot_layout(widths = c(0.6,2,1,0.6))
  
)/(
  ( plot.results.jpeg$umap.IBD2|plot.results.jpeg$umap.EBD2|plot.results.jpeg$umap.nicolBla1|plot.results.jpeg$umap.nBGuo)
  
)/(
  plot.results.jpeg$NHP.IBD2.mnn|plot.results.jpeg$NHP.EBD2.mnn|plot.results.jpeg$NHP.nicolBla.mnn|plot.results.jpeg$NHP.nBGuo.mnn
) +plot_layout(heights = c(2.5,1.75,1.75))

dev.off()



pdf("temp.fig1A_legend.pdf")
plot.results$umap.em.legend
plot.results$CS.legend
plot.results$qr.bind.legend
dev.off()



pdf("temp.fig1C_heatmap.pdf",7,7)
plot.results$MK.heatmap.EM
plot.results$MK.heatmap.blastoids
plot.results$MK.heatmap.blastoids.meanlogExp
dev.off()
pdf("temp.fig1d_pca.pdf",7,4)
plot.results$psd.human.pca
dev.off()


temp.sel.gene <- c("ISL1","GABRP","IGFBP5","GATA2","GCM1","BIN2")
pdf("temp.fig1E.AmnVSTE.Vin.pdf",15,10)
print(
  plot_grid(plotlist=plot.results[paste0("MK.Vln.",temp.sel.gene)],nrow=2,ncol=3)
)
print(
  plot_grid(plotlist=plot.results.jpeg[paste0("MK.Vln.",temp.sel.gene)],nrow=2,ncol=3)
)
dev.off()

pdf("temp.fig1G.model.perf.pdf",7,5)
plot.results$model.perf
dev.off()

pdf("temp.fig1G.blastoids.model.pdf",10,10,family="Arial")
cowplot::plot_grid(plot.results$IBD2.model.perf,plot.results$EBD2.model.per,plot.results$nicolBla.model.perf,plot.results$nBGuo.model.perf)
dev.off()



pdf("temp.supfig1A.pdf",12,5,family="Arial")
plot.results$new.anno.IBD2 +NoLegend()|plot.results$new.anno.EBD2+NoLegend() |plot.results$new.anno.nicolBla+NoLegend() |plot.results$new.anno.nBGuo+NoLegend() 
dev.off()




pdf("temp.sup1A.pdf",20,6)
plot.results$umap.nicolBla2|plot.results$umap.JPF|plot.results$umap.NOJPF.em|plot.results$umap.NOJPF.IBD2
plot.results.jpeg$umap.nicolBla2|plot.results.jpeg$umap.JPF|plot.results.jpeg$umap.NOJPF.em|plot.results.jpeg$umap.NOJPF.IBD2
plot.results.txt$umap.nicolBla2|plot.results.txt$umap.JPF|plot.results.txt$umap.NOJPF.em|plot.results.txt$umap.NOJPF.IBD2
dev.off()

pdf("temp.sup1B.pdf",15,6)
plot.results$NHP.CS7.mnn|plot.results$NHP.CS7.cca|plot.results$NHP.CS7.harmony
plot.results.jpeg$NHP.CS7.mnn|plot.results.jpeg$NHP.CS7.cca|plot.results.jpeg$NHP.CS7.harmony
plot.results.txt$NHP.CS7.mnn|plot.results.txt$NHP.CS7.cca|plot.results.txt$NHP.CS7.harmony
plot.results$NHP.CS7.mnn.spec|plot.results$NHP.CS7.cca.spec|plot.results$NHP.CS7.harmony.spec
dev.off()

pdf("temp.sup1C.pdf",20,6)
plot.results$NHP.IBD2.mnn.spec|plot.results$NHP.EBD2.mnn.spec|plot.results$NHP.nicolBla.mnn.spec|plot.results$NHP.nBGuo.mnn.spec
plot.results$NHP.IBD2.harmony|plot.results$NHP.EBD2.harmony|plot.results$NHP.nicolBla.harmony|plot.results$NHP.nBGuo.harmony
plot.results.jpeg$NHP.IBD2.harmony|plot.results.jpeg$NHP.EBD2.harmony|plot.results.jpeg$NHP.nicolBla.harmony|plot.results.jpeg$NHP.nBGuo.harmony
plot.results.txt$NHP.IBD2.harmony|plot.results.txt$NHP.EBD2.harmony|plot.results.txt$NHP.nicolBla.harmony|plot.results.txt$NHP.nBGuo.harmony
plot.results$NHP.IBD2.harmony.spec|plot.results$NHP.EBD2.harmony.spec|plot.results$NHP.nicolBla.harmony.spec|plot.results$NHP.nBGuo.harmony.spec

plot.results$NHP.IBD2.cca|plot.results$NHP.EBD2.cca|plot.results$NHP.nicolBla.cca|plot.results$NHP.nBGuo.cca
plot.results.jpeg$NHP.IBD2.cca|plot.results.jpeg$NHP.EBD2.cca|plot.results.jpeg$NHP.nicolBla.cca|plot.results.jpeg$NHP.nBGuo.cca
plot.results.txt$NHP.IBD2.cca|plot.results.txt$NHP.EBD2.cca|plot.results.txt$NHP.nicolBla.cca|plot.results.txt$NHP.nBGuo.cca
plot.results$NHP.IBD2.cca.spec|plot.results$NHP.EBD2.cca.spec|plot.results$NHP.nicolBla.cca.spec|plot.results$NHP.nBGuo.cca.spec

dev.off()

pdf("temp.supFig3.reanno.pdf",12,9,family="Arial")
(plot.results$ReCheck.D3post.MK.Vln.NANOG+plot.results$ReCheck.D3post.MK.Vln.BMP2)/(plot.results$ReCheck.D3post.MK.Vln.GATA2+plot.results$ReCheck.D3post.MK.Vln.PMP22) |(plot.results[["CS7.recheck.old.anno"]] +plot.results[["CS7.recheck.new.anno"]])/plot.results$CS7.recheck.mk.INT.nolegend
plot.results$CS7.recheck.mk.INT
dev.off()

pdf("temp.sup1D.pdf")
plot.results$NHP.MK.heatmap.EM
dev.off()
pdf("temp.sup.pca.pdf",7,5.5)
plot.results$psd.NHP.human.pca
dev.off()

pdf("temp.sup.cluster.pdf",7,7)
plot.results$umap.cluster 
plot.results.jpeg$umap.cluster 
dev.off()


temp.sel.gene <- c("LIX1","TMEM88","RGS5","PMP22","VIM")#
pdf("temp.sup.MesMk.Vin.pdf",15,9)
print(
  plot_grid(plotlist=plot.results[c(paste0("Mes.MK.Vln.IBD2.",temp.sel.gene),paste0("Mes.MK.Vln.EBD2.",temp.sel.gene),paste0("Mes.MK.Vln.nicolBla.",temp.sel.gene))],nrow=3,ncol=5)
)
dev.off()
# pdf("D2.AmnVSTE.Vin.pdf",15,5)
# print(
#   plot_grid(plotlist=plot.results[paste0("MK.Vln.",temp.sel.gene)],nrow=2,ncol=3)
# )
# print(
#   plot_grid(plotlist=plot.results.jpeg[paste0("MK.Vln.",temp.sel.gene)],nrow=2,ncol=3)
# )
# dev.off()
#dev.off()
# pdf("D2.UMAP.plot.text.pdf",9,6)
# (
#   plot_spacer()+plot.results$umap.em.tiff1+plot.results$umap.em_legend+plot_spacer()+plot_layout(widths = c(0.6,2,1,0.6))
#   
# )/(
#   (Plot.umap$JPF.tiff1|Plot.umap$EBD2.tiff1|Plot.umap$IBD2.tiff1|Plot.umap$nBGuoids.tiff1)
#   
# )+plot_layout(heights = c(2.5,1.5))
# dev.off()
# pdf("D2.UMAP.plot.point.pdf",9,6)
# (
#   plot_spacer()+plot.results$umap.em.tiff2+plot.results$umap.em_legend+plot_spacer()+plot_layout(widths = c(0.6,2,1,0.6))
#   
# )/(
#   (Plot.umap$JPF.tiff2|Plot.umap$EBD2.tiff2|Plot.umap$IBD2.tiff2|Plot.umap$nBGuoids.tiff2)
#   
# )+plot_layout(heights = c(2.5,1.5))
# dev.off()



data.ob.umap %>% filter(seurat_clusters=="C13") %>% filter(pj %in% c("IBD2","EBD2","nBGuo","nicolBla")) %>% filter(subCT %in% c("None","blastoids","EB_ELC","Gata3_Day3.1","Gata3_Day3.3","Gata3_Day3.9")) %>% pull(pj) %>% table()



temp.c1 <-  data.ob.umap %>% filter(cluster_EML=="AMLC" & pj=="IBD2") %>% pull(cell)
temp.c2 <-  data.ob.umap %>% filter(cluster_EML=="AMLC" & pj=="JPF2019") %>% pull(cell)


temp.exp <-lognormExp.mBN[sub.amnion.mk,c(temp.c1,temp.c2)] 
temp.sel.exp <- t(apply(temp.exp,1,scale))
colnames(temp.sel.exp) <- colnames(temp.exp)
rownames(temp.sel.exp) <- rownames(temp.exp)
temp.sel.exp[temp.sel.exp>zs.limit] <- zs.limit
temp.sel.exp[temp.sel.exp<  (-1*zs.limit)] <- -1*zs.limit

temp.sel.exp%>% pheatmap(scale="row",show_colnames=F,cluster_cols=F,cluster_rows=F,gaps_col=(length(temp.c1)))
data.frame(IBD2_AMLC=rowMeans(expm1(lognormExp.mBN[sub.amnion.mk,c(temp.c1)]) ), JPF2019_AMLC=(rowMeans(expm1(lognormExp.mBN[sub.amnion.mk,c(temp.c2)]))))
data.frame(IBD2_AMLC=rowMeans(expm1(lognormExp.mBN[sub.amnion.mk,c(temp.c1)]) ), JPF2019_AMLC=(rowMeans(expm1(lognormExp.mBN[sub.amnion.mk,c(temp.c2)])))) %>% pheatmap(scale="row",cluster_rows=F)