# -*- snakemake -*-
rule merge_raw_data:
	input:
		"Combine.all.rawData.R" 
	output:
		"Def_PASE.R",
		"Def_Eblastoids.R",
		"Def_Iblastoids.R",
		"Def_nicolasBlastoids.R"

rule QC:
	input:
		"Def_PASE.R",
		"Def_Eblastoids.R",
		"Def_Iblastoids.R",
		"Def_nicolasBlastoids.R"

	output:
		"QC.downsampling.R"

rule normalization:
	input:
		"QC.downsampling.R"
	output:
		"scran.norm.R"

rule  recover_NHP_data:
	input:
		"Combine.all.rawData.R"
	output:
		"Rec.nonHumanPrimate.R" 

rule multi_data_integration:
	input:
		"QC.downsampling.R",
		"scran.norm.R"
	output:
		"Whole.humanDataSet.IT.R"

rule multi_data_integration_without_PASE:
	input:
		"QC.downsampling.R",
		"scran.norm.R"
	output:
		"withoutAMLC.humanDataSet.IT.R"

rule cross_speciese_integration:
	input:
		"QC.downsampling.R",
		"Rec.nonHumanPrimate.R"
	output:
		"CrossSpecies.NHP.IT.R"

rule detect_Amnion_VS_TE_DEG:
	input:
		"scran.norm.R"
	output:
		"singleCell.TEvsAmnion.Marker.R"

rule create_pseudu_bulk:
	input:
		"QC.downsampling.R",
		"Rec.nonHumanPrimate.R"
	output:
		"create.pseudo.bulk.R"

rule pseudo_bulk_pca:
	input:
		"create.pseudo.bulk.R"
	output:
		"pseudo.bulk.PCA.R"

rule check_em_reference_annotation:
	input:
		"Whole.humanDataSet.IT.R"
	output:
		"Check.CS7.subType.R"#

rule update_em_reference_annotation:
	input:
		"Check.CS7.subType.R"
	output:
		"UpDate.current.anno.R"

rule generate_h5_file:
	input:
		"QC.downsampling.R",
		"Rec.nonHumanPrimate.R"
	output:
		"cb_predict_create_ref_pred_h5.R"

rule check_model_performance:
	input:
		"UpDate.current.anno.R",
		"cb_predict_create_ref_pred_h5.R"
	output:
		"check.model.performance.R"

rule generate_figures_ob:
	input:
		"Whole.humanDataSet.IT.R",
		"withoutAMLC.humanDataSet.IT.R",
		"CrossSpecies.NHP.IT.R",
		"singleCell.TEvsAmnion.Marker.R",
 		"Check.CS7.subType.R",
 		"check.model.performance.R",
		"figures.setting.R"
	output:
		"Figure.generate.R"

rule plot_figures:
	input:
		"Figure.generate.R"
	output:
		"fig.plotout.R"


#snakemake plot_figures --rulegraph -npF -s ../pipeline/analysis.pipeline.smk| dot -Tpng  > ../pipeline.png






