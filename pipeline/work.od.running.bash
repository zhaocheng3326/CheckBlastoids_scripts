#!/bin/bash
DIR=/home/chenzh/My_project/test_smncRNA
TMPD=$DIR/tmp_data
SRC=$DIR/src
DATA=$DIR/data
RE=$DIR/results
DOC=$DIR/doc
BIN=$DIR/bin

BDOC=$DIR/big_doc

cd $DIR



~/miniconda3/envs/R3.6/bin/Rscript -e 'source("src/Def_PASE.R")'
~/miniconda3/envs/R3.6/bin/Rscript -e 'rmarkdown::render("src/Def_PASE.R", output_dir = "html_report")'
~/miniconda3/envs/R3.6/bin/Rscript -e 'rmarkdown::render("src/Def_Eblastoids.R", output_dir = "html_report")'
~/miniconda3/envs/R3.6/bin/Rscript -e 'rmarkdown::render("src/Def_Iblastoids.R", output_dir = "html_report")'
~/miniconda3/envs/R4.0/bin/Rscript -e 'rmarkdown::render("src/Def_nicolasBlastoids.R", output_dir = "html_report")'
~/miniconda3/envs/R3.6/bin/Rscript -e 'source("src/QC.downsampling.R")'
~/miniconda3/envs/R3.6/bin/Rscript -e 'source("src/scran.norm.R")'
~/miniconda3/envs/R4.0/bin/Rscript -e 'source("src/Rec.nonHumanPrimate.R")'
~/miniconda3/envs/R3.6/bin/Rscript -e 'rmarkdown::render("src/Whole.humanDataSet.IT.R", output_dir = "html_report")'
~/miniconda3/envs/R3.6/bin/Rscript -e 'rmarkdown::render("src/withoutAMLC.humanDataSet.IT.R", output_dir = "html_report")'
~/miniconda3/envs/R4.0/bin/Rscript -e 'rmarkdown::render("src/CrossSpecies.NHP.IT.R", output_dir = "html_report")'
~/miniconda3/envs/R4.0/bin/Rscript -e 'source("src/singleCell.TEvsAmnion.Marker.R")'
~/miniconda3/envs/R4.0/bin/Rscript -e 'source("src/create.pseudo.bulk.R")'
~/miniconda3/envs/R3.6/bin/Rscript -e 'rmarkdown::render("src/pseudo.bulk.PCA.R", output_dir = "html_report")'
~/miniconda3/envs/R4.0/bin/Rscript -e 'rmarkdown::render("src/Check.CS7.subType.R", output_dir = "html_report")'
~/miniconda3/envs/R4.0/bin/Rscript -e 'source("src/UpDate.current.anno.R")'
~/miniconda3/envs/R4.0/bin/Rscript -e 'source("src/cb_predict_create_ref_pred_h5.R")'
~/miniconda3/envs/R4.0/bin/Rscript -e 'rmarkdown::render("src/check.model.performance.R", output_dir = "html_report")'
~/miniconda3/envs/R4.0/bin/Rscript -e 'rmarkdown::render("src/Figure.generate.R", output_dir = "html_report")'


#rsync -v -PRatrl ~/My_project/JP_project/./tmp_data/Prev.old.M.fromJP.Rdata /san/Sophie_ongoing/Cheng_tmpData/My_project/CheckBlastoids/
#rsync -v -PRatr ~/My_project/JP_project/./tmp_data/CP_data/merge.rsem_counts.csv /san/Sophie_ongoing/Cheng_tmpData/My_project/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./tmp_data/CP_data/merge.rsem_rpkm.csv /san/Sophie_ongoing/Cheng_tmpData/My_project/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./tmp_data/CP_data/merge.rsem_tpm.csv /san/Sophie_ongoing/Cheng_tmpData/My_project/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./tmp_data/GSE136447/merge.rsem_counts.csv /san/Sophie_ongoing/Cheng_tmpData/My_project/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./tmp_data/GSE136447/merge.rsem_rpkm.csv /san/Sophie_ongoing/Cheng_tmpData/My_project/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./tmp_data/GSE136447/merge.rsem_tpm.csv /san/Sophie_ongoing/Cheng_tmpData/My_project/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./tmp_data/Blakeley_Data/merge.rsem_counts.csv /san/Sophie_ongoing/Cheng_tmpData/My_project/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./tmp_data/nBlastoids//merge.rsem_counts.csv /san/Sophie_ongoing/Cheng_tmpData/My_project/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./tmp_data/IBlastoids/out/IBD2/outs/filtered_feature_bc_matrix /san/Sophie_ongoing/Cheng_tmpData/My_project/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./tmp_data/EBlastoids/out/LW36/outs/filtered_feature_bc_matrix /san/Sophie_ongoing/Cheng_tmpData/My_project/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./tmp_data/EBlastoids/out/LW60/outs/filtered_feature_bc_matrix /san/Sophie_ongoing/Cheng_tmpData/My_project/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./tmp_data/EBlastoids/out/LW61/outs/filtered_feature_bc_matrix /san/Sophie_ongoing/Cheng_tmpData/My_project/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./tmp_data/Amnion_EPSC/H9_Amnion/outs/filtered_feature_bc_matrix /san/Sophie_ongoing/Cheng_tmpData/My_project/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./tmp_data/Amnion_EPSC/Posterior48h_1/outs/filtered_feature_bc_matrix /san/Sophie_ongoing/Cheng_tmpData/My_project/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./tmp_data/human_gastrula/merge.rsem_counts.csv /san/Sophie_ongoing/Cheng_tmpData/My_project/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./doc/stirparo2018_tableS4.xlsx ~/My_project/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./data/GSE171820/GSE171820.anno.txt /lanner-lab/sophie-lab-backup/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./data/GSE177689/CR_RefSeq/merge.rsem_counts.csv /lanner-lab/sophie-lab-backup/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./data/GSE177689/blastoid.H9.okae_bts5__scRNAseq.unfiltered__metadata.tsv /lanner-lab/sophie-lab-backup/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./data/Gastrulation/annot_umap.rds /lanner-lab/sophie-lab-backup/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./data/GSE136447/SraRunInfo.csv /lanner-lab/sophie-lab-backup/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./data/GSE136447/GSE136447.tsv /lanner-lab/sophie-lab-backup/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./data/GSE136447/41586_2019_1875_MOESM10_ESM.tsv /lanner-lab/sophie-lab-backup/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./data/GSE134571/GSE134571_Posterior48h_H9_Amnion_Merged.rds /lanner-lab/sophie-lab-backup/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./data/stirparo2018_tableS4.xlsx /lanner-lab/sophie-lab-backup/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./data/Primate_10X_Amn/raw_data/mart_export_idtranslate_fun.txt /lanner-lab/sophie-lab-backup/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./data/Meistermann_et_al/1-s2.0-S1934590921001855-mmc2.xlsx /lanner-lab/sophie-lab-backup/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./data/Castel_et_al_2020-master/Dataset/ZhouPetro/sampleAnnot.tsv /lanner-lab/sophie-lab-backup/CheckBlastoids/
#rsync -v -PRatrl ~/My_project/JP_project/./data/Castel_et_al_2020-master/Results/sampleAnnot.mergedTechReplicate.tsv /lanner-lab/sophie-lab-backup/CheckBlastoids/


