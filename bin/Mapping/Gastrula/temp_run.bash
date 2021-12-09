
CRAM_ref=~/My_project/JP_project/data/Gastrulation/Mus_musculus.GRCm38.68.dna.toplevel.fa
RD=/home/chenzh/My_project/JP_project/data/Gastrulation/raw_data/TEMP
SA=sc7786109
GD_star=/home/chenzh/Genome/Human/refdata-cellranger-GRCh38-3.0.0/star
GF=/home/chenzh/Genome/Human/refdata-cellranger-GRCh38-3.0.0/fasta/genome
GTF=/home/chenzh/Genome/Human/refdata-cellranger-GRCh38-3.0.0/genes/genes.gtf
RS_GF=/home/chenzh/Genome/Human/refdata-cellranger-GRCh38-3.0.0/rsem
core=20


samtools sort  -n SS-${SA}.cram  -o ${SA}.sorted.cram --reference  $CRAM_ref
samtools fastq  -1 ${SA}_R1.fastq -2 ${SA}_R2.fastq ${SA}.sorted.cram --reference  $CRAM_ref


mkdir -p ${RD}/${SA}
STAR --runThreadN ${core} --genomeDir ${GD_star} --readFilesIn ${RD}/${SA}_R1.fastq ${RD}/${SA}_R2.fastq  --quantMode TranscriptomeSAM --outSAMstrandField intronMotif  --outSAMtype BAM Unsorted  --outFileNamePrefix ${RD}/${SA}/${SA}.

rsem-calculate-expression --no-bam-output --alignments --paired-end --single-cell-prior -p ${core} ${RD}/${SA}/${SA}.Aligned.toTranscriptome.out.bam ${RS_GF} ${RD}/${SA}/${SA}.rsem -q