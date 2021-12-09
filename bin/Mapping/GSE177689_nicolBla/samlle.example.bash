#configfile: "All_dataset.config.yaml"

#workdir: "/home/soppet/My_project/ScRNA_analysis/"


cp -rs ~/My_project/JP_project/data/GSE136447/ ./



#java -jar /home/czhao/tools/Trimmomatic-0.36/trimmomatic-0.36.jar PE -threads 10 -phred33 -trimlog  $TMPD/clean_data/${ARRAY[$SGE_TASK_ID]}/trim.log  $DATA/raw_data/${ARRAY[$SGE_TASK_ID]}_R1.fastq $DATA/raw_data/${ARRAY[$SGE_TASK_ID]}_R2.fastq -baseout $TMPD/clean_data/${ARRAY[$SGE_TASK_ID]}/${ARRAY[$SGE_TASK_ID]}.fastq ILLUMINACLIP:$ADA:2:30:10 LEADING:3 SLIDINGWINDOW:4:30 MINLEN:36 2> $TMPD/clean_data/${ARRAY[$SGE_TASK_ID]}/Trim.sta


RD=/home/chenzh/My_project/JP_project/tmp_data/GSE136447/
SA=SRR10039951
GD_star=/home/chenzh/Genome/Human/refdata-cellranger-GRCh38-3.0.0/star
GF=/home/chenzh/Genome/Human/refdata-cellranger-GRCh38-3.0.0/fasta/genome
GTF=/home/chenzh/Genome/Human/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf
RS_GF=/home/chenzh/Genome/Human/refdata-cellranger-GRCh38-3.0.0/rsem
core=20

fastq-dump --split-3  ${RD}/${SA}.sra --gzip  -O ${RD}
# fastqc check
mkdir -p ${RD}/${SA}
STAR --runThreadN ${core} --genomeDir ${GD_star} --readFilesIn ${RD}/${SA}_1.fastq.gz ${RD}/${SA}_2.fastq.gz  --quantMode TranscriptomeSAM --outSAMstrandField intronMotif --readFilesCommand zcat --outSAMtype BAM Unsorted  --outFileNamePrefix ${RD}/${SA}/${SA}.

rsem-calculate-expression --no-bam-output --alignments --paired-end --single-cell-prior -p ${core} ${RD}/${SA}/${SA}.Aligned.toTranscriptome.out.bam ${RS_GF} ${RD}/${SA}/${SA}.rsem -q






for I in {1..10}
do
  mkdir -p sub_${I}
  cp cp.quant.map.smk  sub_${I}
  echo "RD: /home/zhaoch/snic2020-15-102/nobackup/CP_data/sub_${I}" >sub_${I}/cp.yaml
  echo "SampleName_file: /home/zhaoch/My_project/JP_project/bin/CP_data/sub_${I}/Sample.name.file" >> sub_${I}/cp.yaml
  tail -n 4 cp.yaml >> sub_${I}/cp.yaml
  mv x0${I}Sample.name sub_${I}/Sample.name.file
done


python ~/My_project/JP_project/src/merge_rsem_data.py -i ERR1041403/ERR1041403.rsem.genes.results ERR1041404/ERR1041404.rsem.genes.results ERR1041405/ERR1041405.rsem.genes.results ERR1041406/ERR1041406.rsem.genes.results -o ~/My_project/JP_project/tmp_data/CP_data/merge -f .rsem.genes.results