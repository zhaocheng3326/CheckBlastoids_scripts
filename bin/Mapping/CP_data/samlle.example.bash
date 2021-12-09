#configfile: "All_dataset.config.yaml"

#workdir: "/home/soppet/My_project/ScRNA_analysis/"


cp -rs ~/My_project/JP_project/data/CP_data/ ./



#java -jar /home/czhao/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 10 -phred33 -trimlog  $TMPD/clean_data/${ARRAY[$SGE_TASK_ID]}/trim.log  $DATA/raw_data/${ARRAY[$SGE_TASK_ID]}_R1.fastq $DATA/raw_data/${ARRAY[$SGE_TASK_ID]}_R2.fastq -baseout $TMPD/clean_data/${ARRAY[$SGE_TASK_ID]}/${ARRAY[$SGE_TASK_ID]}.fastq ILLUMINACLIP:$ADA:2:30:10 LEADING:3 SLIDINGWINDOW:4:30 MINLEN:36 2> $TMPD/clean_data/${ARRAY[$SGE_TASK_ID]}/Trim.sta


RD=/home/chenzh/My_project/JP_project/tmp_data/CP_data/
SA=ERR1041404
GD_star=/home/chenzh/Genome/Human/refdata-cellranger-GRCh38-3.0.0/star
GF=/home/chenzh/Genome/Human/refdata-cellranger-GRCh38-3.0.0/fasta/genome
GTF=/home/chenzh/Genome/Human/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf
RS_GF=/home/chenzh/Genome/Human/refdata-cellranger-GRCh38-3.0.0/rsem
core=20

#STAR --runThreadN 4 --runMode genomeGenerate --genomeDir ${GD_star} --genomeFastaFiles ${GF} --genomeSAindexNbases 14 --sjdbGTFfile ${GTF} --sjdbOverhang 42 --outFileNamePrefix $(dirname doc/ref/star_index/Genome)/ > doc/ref/star_index/Genome.log
mkdir -p ${RD}/${SA}
STAR --runThreadN ${core} --genomeDir ${GD_star} --readFilesIn ${RD}/${SA}.fastq.gz --quantMode TranscriptomeSAM --outSAMstrandField intronMotif --readFilesCommand zcat --outSAMtype BAM Unsorted  --outFileNamePrefix ${RD}/${SA}/${SA}.





rsem-calculate-expression --alignments -p ${core} ${RD}/${SA}/${SA}.Aligned.toTranscriptome.out.bam ${RS_GF} ${RD}/${SA}/${SA}.rsem -q


python ~/My_project/JP_project/src/merge_rsem_data.py -i ERR1041403/ERR1041403.rsem.genes.results ERR1041404/ERR1041404.rsem.genes.results ERR1041405/ERR1041405.rsem.genes.results ERR1041406/ERR1041406.rsem.genes.results -o ~/My_project/JP_project/tmp_data/CP_data/merge -f .rsem.genes.results