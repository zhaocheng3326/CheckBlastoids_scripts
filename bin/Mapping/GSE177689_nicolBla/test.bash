1.lts_workflows_sm_scrnaseq/rules/utils/dbutils_fasta_to_genbank.smk:rule dbutils_fasta_to_genbank:
 input: doc/ref/ERCC.fasta
 output: doc/ref/ERCC.genbank
#Note use the biopython to transfer the fasta to gen.bank

2. rules/utils/dbutils_make_transcript_annot_gtf.smk: rule dbutils_make_transcript_annot_gtf:
    input: doc/ref/Homo_sapiens.GRCh38.89.gtf, doc/ref/ERCC.genbank
    output: doc/ref/Homo_sapiens.GRCh38.89-ERCC.gtf
# Transfer the ERCC.geneBank as gtf information into the Homo_sapiens.GRCh38.89.gtf file

3.rule ucsc/ucsc_gtfToGenepred:
    input: doc/ref/Homo_sapiens.GRCh38.89-ERCC.gtf
    output: doc/ref/Homo_sapiens.GRCh38.89-ERCC.genePred
# order: gtfToGenePred -genePredExt -ignoreGroupsWithoutExons doc/ref/Homo_sapiens.GRCh38.89-ERCC.gtf doc/ref/Homo_sapiens.GRCh38.89-ERCC.genePred
## Half of the file  that the col6(coding start region)== col7 (coding end region) ## cause those genes don't have the UTR region

4.ucsc/gtf_to_bed12:
    input: doc/ref/Homo_sapiens.GRCh38.89-ERCC.gtf
    output: doc/ref/Homo_sapiens.GRCh38.89-ERCC.bed12
# order perl scripts/gtf2bed.pl Homo_sapiens.GRCh38.89-ERCC.gtf > Homo_sapiens.GRCh38.89-ERCC.bed12


5.localrule ucsc/ucsc_genepred_to_refFlat:
    input: doc/ref/Homo_sapiens.GRCh38.89-ERCC.genePred
    output: doc/ref/Homo_sapiens.GRCh38.89-ERCC.refFlat
#order awk '{printf("%s\t%s\n", $1, $0)}'  doc/ref/Homo_sapiens.GRCh38.89-ERCC.genePred > doc/ref/Homo_sapiens.GRCh38.89-ERCC.refFlat
## Half of the file,  col7(coding start region)== col8 (coding end region)

6. rule star_index:
    input: doc/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa, doc/ref/ERCC.fasta, doc/ref/Homo_sapiens.GRCh38.89-ERCC.gtf
    output: doc/ref/star_index/
STAR --runThreadN 20 --runMode genomeGenerate --genomeDir ./doc/ref/star_index --genomeFastaFiles doc/ref/Homo_sapiens.GRCh38.dna.primary_assembly.fa doc/ref/ERCC.fasta --genomeSAindexNbases 14 --sjdbGTFfile doc/ref/Homo_sapiens.GRCh38.89-ERCC.gtf --sjdbOverhang 42 --outFileNamePrefix $(dirname doc/ref/star_index/Genome)/ > doc/ref/star_index/Genome.log
### --sjdbOverhang 42 (should be the max(length)-1

7.rule star_align_se:
    input: data/samples/ERR1042769/ERR1042769.fastq.gz, doc/ref/star_index/Genome, doc/ref/star_index/SA, doc/ref/star_index/SAindex, doc/ref/star_index/genomeParameters.txt, doc/ref/star_index/chrLength.txt, doc/ref/star_index/chrName.txt, doc/ref/star_index/chrNameLength.txt, doc/ref/star_index/chrStart.txt, doc/ref/star_index/geneInfo.tab, doc/ref/star_index/sjdbInfo.txt, doc/ref/star_index/sjdbList.fromGTF.out.tab, doc/ref/star_index/sjdbList.out.tab, doc/ref/star_index/transcriptInfo.tab
    output: data/samples/ERR1042769/ERR1042769.Aligned.out.bam, data/samples/ERR1042769/ERR1042769.Log.final.out

# ORDER STAR --runThreadN 5 --genomeDir ./doc/ref/star_index --readFilesIn data/samples/ERR1042769/ERR1042769.fastq.gz --outSAMstrandField intronMotif --readFilesCommand zcat --outSAMtype BAM Unsorted  --outFileNamePrefix data/samples/ERR1042769/ERR1042769.

8.rule bamtools/bamtools_filter_unique:
    input: data/samples/ERR1042578/ERR1042578.Aligned.out.bam
    output: data/samples/ERR1042578/ERR1042578.Aligned.out_unique.bam
# ORDER: bamtools filter -in data/samples/ERR1042769/ERR1042769.Aligned.out.bam -out data/samples/ERR1042769/ERR1042769.Aligned.out_unique.bam -mapQuality ">=255"
## The mapping quality MAPQ (column 5) is 255 for uniquely mapping reads, and int(-10*log10(1-1/Nmap)) for multi-mapping reads. This scheme is same as the one used by TopHat and is compatible
with Cufflinks

9. rule picard/picard_merge_sam:
    input: data/samples/ERR1042769/ERR1042769.Aligned.out_unique.bam
    output: data/samples/ERR1042769/ERR1042769.merge.bam

## Not very sure deatils ,seems the order of reads in sam files is changed.

10.rule rpkmforgenes_from_bam:
    input: doc/ref/Homo_sapiens.GRCh38.89-ERCC.refFlat, data/samples/ERR1042769/ERR1042769.merge.bam
    output: data/samples/ERR1042769/ERR1042769.merge.rpkmforgenes

# order rpkmforgenes.py -readcount -fulltranscript -mRNAnorm -rmnameoverlap -bothendsceil   -a doc/ref/Homo_sapiens.GRCh38.89-ERCC.refFlat -bamu -i data/samples/ERR1042769/ERR1042769.merge.bam -o data/samples/ERR1042769/ERR1042769.merge.rpkmforgenes &> /dev/null
##-readcount :to add the number of reads to the output -fulltranscript to not remove 3'UTRs (default) -mRNAnorm:to normalize by the number of reads matching mRNA exons (default) -rmnameoverlap to ignore regions shared my multiple genes -bothendsceil to set -bothends but round the read count upward


11.  rule picard/rule picard_sort_sam:
input: data/samples/ERR1042578/ERR1042578.merge.bam
output: data/samples/ERR1042578/ERR1042578.merge.sort.bam, data/samples/ERR1042578/ERR1042578.merge.sort.bai, data/samples/ERR1042578/ERR1042578.merge.sort.bam.bai

# order picard -Xmx8g -Djava.io.tmpdir=/tmp  SortSam I=data/samples/ERR1042578/ERR1042578.merge.bam O=data/samples/ERR1042578/ERR1042578.merge.sort.bam SORT_ORDER=coordinate VALIDATION_STRINGENCY=SILENT CREATE_INDEX=TRUE;
#ln data/samples/ERR1042578/ERR1042578.merge.sort.bai data/samples/ERR1042578/ERR1042578.merge.sort.bam.bai
# sort reads by position

12. rule rseqc/rseqc_read_distribution:
    input: data/samples/ERR1042578/ERR1042578.merge.bam, doc/ref/Homo_sapiens.GRCh38.89-ERCC.bed12
    output: data/samples/ERR1042578/ERR1042578.merge_rseqc/read_distribution.txt

# order read_distribution.py  -i data/samples/ERR1042578/ERR1042578.merge.bam -r doc/ref/Homo_sapiens.GRCh38.89-ERCC.bed12 &> data/samples/ERR1042578/ERR1042578.merge_rseqc/read_distribution.txt
## upstrem, intron, 5k and soon

13.rule featurecounts/featurecounts:
    input: doc/ref/Homo_sapiens.GRCh38.89-ERCC.gtf, data/samples/ERR1042578/ERR1042578.merge.sort.bam
    output: data/samples/ERR1042578/ERR1042578.merge.sort.biotypes
#order: featureCounts -g gene_biotype -a doc/ref/Homo_sapiens.GRCh38.89-ERCC.gtf -o temp_fc data/samples/ERR1042578/ERR1042578.merge.sort.bam > /dev/null
#       cut -f 1,7 temp_fc > data/samples/ERR1042578/ERR1042578.merge.sort.biotypes
#A read is said to overlap a feature if at least one read base is found to overlap the feature. For paired-end data, a fragment (or template) is said to overlap a feature if any of the two reads from that fragment is found to overlap the feature.does not count reads overlapping with more than one feature use -o to change.
# use gene_biotype as features

14. localrule multiqc/multiqc.smk biotypes_to_multiqc:
    input: data/samples/ERR1042578/ERR1042578.merge.sort.biotypes
    output: data/samples/ERR1042578/ERR1042578.merge.sort.biotypes_mqc.txt
### combine rules in multiqc.setting.smk

15. rule rseqc_geneBody_coverage:
    input: data/samples/ERR1042578/ERR1042578.merge.sort.bam, data/samples/ERR1042578/ERR1042578.merge.sort.bam.bai, doc/ref/hg38.HouseKeepingGenes.bed
    output: data/samples/ERR1042578/ERR1042578.merge_rseqc/geneBody_coverage.geneBodyCoverage.txt
# order geneBody_coverage.py  -i data/samples/ERR1042578/ERR1042578.merge.sort.bam -o $(dirname data/samples/ERR1042578/ERR1042578.merge_rseqc/geneBody_coverage.geneBodyCoverage.txt)/geneBody_coverage -r doc/ref/hg38.HouseKeepingGenes.bed 2> /dev/null calculate 100 percentage in gene location

16. localrule multiqc_file_list: multiqc/multiqc.smk
    input: use all samles'geneBody_coverage.geneBodyCoverage.txt read_distribution.txt ERR1042406.merge.sort.biotypes_mqc.txt ERR1042406.Log.final.out
    output: multiqc/multiqc_report.html.input.all.txt

## 1.Reads distribution in 100 percentage of GeneBody.  2. Upstream 5k, Downstream 5k,  utr regions of reads distribution 3. Reads distribution in biotypes distribution 4. Final out is the mapping stat.out


17.rule multiqc:
    input: multiqc/multiqc_report.html.input.all.txt
    output: multiqc/multiqc_report.all.html, and all files in multiqc/multiqc_data_all (mqc_featureCounts_biotype_plot_1.txt mqc_star_alignment_plot_1.txt multiqc_rseqc_read_distribution.txt mqc_rseqc_gene_body_coverage_plot_Counts.txt multiqc_data.json multiqc_sources.txt mqc_rseqc_gene_body_coverage_plot_Percentages.txt multiqc_general_stats.txt multiqc_star.txt mqc_rseqc_read_distribution_plot_1.txt) ### actually they all made by multiqc
    jobid: 1
    wildcards: subset=all


        # This is to avoid problems if running parallel jobs. Shadow would be better
        # but requires absolute paths in the input file.
###
        cut -f2 multiqc/multiqc_report.html.input.all.txt > $TMPDIR/file_list
        multiqc -f -d -dd 1 --flat -l $TMPDIR/file_list -o $TMPDIR
        mv $TMPDIR/multiqc_report.html multiqc/multiqc_report.all.html


        # This is to avoid problems if running parallel jobs. Shadow would be better
        # but requires absolute paths in the input file.

        rm -rf multiqc/multiqc_data_all
        TMPDIR=`mktemp -d tmp.XXXXXX`
        cut -f2 multiqc/multiqc_report.html.input.all.txt > $TMPDIR/file_list
        multiqc -f -d -dd 1 --flat -l $TMPDIR/file_list -o $TMPDIR
        mv $TMPDIR/multiqc_report.html multiqc/multiqc_report.all.html
        mv $TMPDIR/multiqc_data multiqc/multiqc_data_all
####


18#. rule rule make_rpkmforgenes_count_matrix:
    input: data/samples/ERR1042406/ERR1042406.merge.rpkmforgenes
    output: results/merge.rpkmforgenes_counts.csv, results/merge.rpkmforgenes_rpkm.csv, results/merge.rpkmforgenes_genes.csv
order:python rules/utils/scripts/merge_rpkmforgenes_data.py -i data/samples/ERR1042406/ERR1042406.merge.rpkmforgenes -o results/merge -f .merge.rpkmforgenes


##Alternative
rules/utils/custom_qc.smk
input: results/merge.rpkmforgenes_rpkm.csv
output: results/merge.rpkmforgenes_rpkm.qc_stats.csv
order:
python lts_workflows_sm_scrnaseq/rules/utils/scripts/gene_expression_qc_stat.py -i merge.rpkmforgenes_rpkm.csv -o merge.rpkmforgenes_rpkm.qc_stats.csv

## count_gene      count_spike     gene_correlation_max    gene_detection  spike_correlation_max   spike_detection Use ERCC as default

#localrule rules/utils/custom_qc.smk merge_all_qc
order
input: results/merge.rpkmforgenes_rpkm.qc_stats.csv, multiqc/multiqc_report.all.html, multiqc/multiqc_data_all
output: results/qc_summary.csv
python /home/soppet/miniconda3/envs/work/lib/python3.6/site-packages/lts_workflows_sm_scrnaseq/rules/utils/scripts/parse_all_qc_stats.py -i results/merge.rpkmforgenes_rpkm.qc_stats.csv --multiqc multiqc/multiqc_data_all -o results/qc_summary.csv -c temp_config.yml

#localrule rules/utils/custom_qc.smk make_qc_report
order:
input: results/qc_summary.csv, results/merge.rpkmforgenes_rpkm.csv, results/merge.rpkmforgenes_counts.csv
output: results/qc_report.html, results/qc_report.settings.yaml, results/qc_report.filtered_cells.txt, results/qc_report.conda_env.yaml, results/qc_report.samplefile.csv, results/qc_report.sceset.RData

## use Rmarkdown get html . core code("utils/scripts/qc_summary_report.Rmd",) output_dir="results",output_file="qc_report.html",params=list(workdir=workdir,qc_summary="results/qc_summary.csv",rpkm_file="results/merge.rpkmforgenes_rpkm.csv",counts_file="results/merge.rpkmforgenes_counts.csv",config_file="results/qc_report.settings.yaml",sample_file="/home/soppet/My_project/scrnaseq_course/tutorial_analysis/doc/Meta_data_human_embryo.csv",output_sceset="results/qc_report.sceset.RData",output_filtered="results/qc_report.filtered_cells.txt")
#
echo 'workdir <- getwd();rmarkdown::render("qc_summary_report.Rmd",output_dir="results",output_file="qc_report.html",params=list(workdir=workdir,qc_summary="results/qc_summary.csv",rpkm_file="results/merge.rpkmforgenes_rpkm.csv",counts_file="results/merge.rpkmforgenes_counts.csv",config_file="results/qc_report.settings.yaml",sample_file="/home/soppet/My_project/scrnaseq_course/tutorial_analysis/doc/Meta_data_human_embryo.csv",output_sceset="results/qc_report.sceset.RData",output_filtered="results/qc_report.filtered_cells.txt"))'| R --vanilla







####################


cd /home/chenzh/Genome/Human/CR_tutorial

cellranger mkfastq --id=tutorial_walk_through --run=cellranger-tiny-bcl-1.2.0  --csv=cellranger-tiny-bcl-simple-1.2.0.csv
#/home/chenzh/Genome/Human/CR_tutorial/tutorial_walk_through/outs/fastq_path/H35KCBCXY

cd ~/Genome/Human/CR_tutorial/run_cellranger_count/
cellranger count --id=run_count_1kpbmcs --fastqs=/home/chenzh/Genome/Human/CR_tutorial/run_cellranger_count/pbmc_1k_v3_fastqs/ --sample=pbmc_1k_v3 --transcriptome=/home/chenzh/Genome/Human/refdata-cellranger-GRCh38-3.0.0 --localcores 20  --localmem 32

cd ~/Genome/Human/CR_tutorial/run_cellranger_aggr/

ellranger reanalyze --id=AGG123_reanalysis  --matrix=AGG123/outs/filtered_feature_bc_matrix.h5  --params=AGG123_reanalysis.csv

devtools::install_github("xuzhougeng/scCATCH")
#https://github.com/dynverse/dynbenchmark.git
#library(celltalker)
#https://arc85.github.io/celltalker/articles/celltalker.html
