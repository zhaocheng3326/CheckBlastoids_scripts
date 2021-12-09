# -*- snakemake -*-
#import os.path
#mkdir -p ${RD}/${SA}
Samples=open(config['SampleName_file'],"r").read().splitlines()
rule sort_cram:
  input:
    cram=config['RD']+"SS-{sample}.cram",
    cram_ref=config['CRAM_ref']
  output:
    cramSort=temp(config['RD']+"{sample}.sorted.cram"),
  shell:
    "samtools sort -n {input.cram}  -o {output.cramSort} --reference  {input.cram_ref}"

rule cram_fq:
  input:
    cramSort=temp(config['RD']+"{sample}.sorted.cram"),
  output:
    fq1=temp(config['RD']+"{sample}_R1.fastq"),
    fq2=temp(config['RD']+"{sample}_R2.fastq"),
  shell:
    "samtools fastq  -1 {output.fq1} -2 {output.fq2} {input.cramSort}"

rule STAR_mapping:
  input:
    fq1=config['RD']+"{sample}_R1.fastq",
    fq2=config['RD']+"{sample}_R2.fastq",
    GTF=config['GTF'],
    GD_star_w=config['GD_star']+"sjdbList.fromGTF.out.tab",
  output:
    bam=config['RD']+"{sample}/{sample}.Aligned.out.bam",
    tbam=temp(config['RD']+"{sample}/{sample}.Aligned.toTranscriptome.out.bam"),
    mlog=config['RD']+"{sample}/{sample}.Log.final.out",
  #threads:4
  shell:
    "STAR --runThreadN 4 --outFilterMultimapNmax 1 --genomeDir $(dirname {input.GD_star_w}) --readFilesIn {input.fq1} {input.fq2} --quantMode TranscriptomeSAM --outSAMstrandField intronMotif  --outSAMtype BAM Unsorted  --outFileNamePrefix $(dirname {output.bam})/$(basename $(dirname {output.bam}))."

rule rsem_exp:
  input:
    bam=config['RD']+"{sample}/{sample}.Aligned.toTranscriptome.out.bam",
    RS_GF=config['RS_GF']+".grp"
  output:
    gcount=config['RD']+"{sample}/{sample}.rsem.genes.results",
  #threads:4
  shell:
    "rsem-calculate-expression --no-bam-output --paired-end --single-cell-prior --alignments -p 4 {input.bam} $(dirname {input.RS_GF})/rsem $(dirname {output.gcount})/$(basename $(dirname {output.gcount})).rsem -q"

rule merge:
  input:
    gcount=expand(config['RD']+"{sample}/{sample}.rsem.genes.results", sample=Samples),
    msrc=config['merge_src']
  output:
    mcount=config['RD']+"merge.rsem_counts.csv",
    mtpm=config['RD']+"merge.rsem_tpm.csv",
  shell:
    "python3.6 {input.msrc} -i {input.gcount} -o $(dirname {output.mcount})/merge -f .rsem.genes.results"
