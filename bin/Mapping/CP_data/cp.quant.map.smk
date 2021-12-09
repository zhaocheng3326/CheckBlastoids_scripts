# -*- snakemake -*-
#import os.path
#mkdir -p ${RD}/${SA}
Samples=open(config['SampleName_file'],"r").read().splitlines()
rule STAR_mapping:
  input:
    fq=config['RD']+"{sample}.fastq.gz",
    GTF=config['GTF'],
    GD_star_w=config['GD_star']+"sjdbList.fromGTF.out.tab",
  output:
    bam=temp(config['RD']+"{sample}/{sample}.Aligned.out.bam"),
    tbam=temp(config['RD']+"{sample}/{sample}.Aligned.toTranscriptome.out.bam"),
    mlog=config['RD']+"{sample}/{sample}.Log.final.out",
  #threads:4
  shell:
     "STAR --runThreadN 4 --outFilterMultimapNmax 1 --genomeDir $(dirname {input.GD_star_w}) --readFilesIn {input.fq} --quantMode TranscriptomeSAM --outSAMstrandField intronMotif --readFilesCommand zcat --outSAMtype BAM Unsorted  --outFileNamePrefix $(dirname {output.bam})/$(basename $(dirname {output.bam}))."

rule rsem_exp:
  input:
    bam=config['RD']+"{sample}/{sample}.Aligned.toTranscriptome.out.bam",
    RS_GF=config['RS_GF']+".grp"
  output:
    gcount=config['RD']+"{sample}/{sample}.rsem.genes.results",
  #threads:4
  shell:
    "rsem-calculate-expression --no-bam-output --single-cell-prior --alignments -p 4 {input.bam} $(dirname {input.RS_GF})/rsem $(dirname {output.gcount})/$(basename $(dirname {output.gcount})).rsem -q"

rule merge:
  input:
    gcount=expand(config['RD']+"{sample}/{sample}.rsem.genes.results", sample=Samples),
    msrc=config['merge_src']
  output:
    mcount=config['RD']+"merge.rsem_counts.csv",
    mtpm=config['RD']+"merge.rsem_tpm.csv",
  shell:
    "python {input.msrc} -i {input.gcount} -o $(dirname {output.mcount})/merge -f .rsem.genes.results"