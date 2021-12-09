# -*- snakemake -*-
#import os.path
#mkdir -p ${RD}/${SA}
Samples=dict()
with open(config['SampleName_file'],"r") as SampleName_file:
  for line in SampleName_file:
    sa=line.strip("\n").split("\t")[0]
    fa=line.strip("\n").split("\t")[1]
    Samples[sa]=Samples.get(sa,list())
    Samples[sa].append(fa)
#Samples=open(config['SampleName_file'],"r").read().splitlines()
rule merge_files:
  input:
    fq1=lambda wildcards:[config['RD']+ "/" + fa + "_1.fastq.gz" for fa in Samples[wildcards.sample]],
    fq2=lambda wildcards:[config['RD']+ "/" + fa + "_2.fastq.gz" for fa in Samples[wildcards.sample]],

  output:
    merge_fq1=temp(config['RD']+"{sample}_R1.fastq.gz"),
    merge_fq2=temp(config['RD']+"{sample}_R2.fastq.gz"),

  run:
    if (len(input.fq1) >1):
      input_file_fq1=" ".join(input.fq1)
      input_file_fq2=" ".join(input.fq2)
      shell("cat {ipf} > {opf}".format(ipf=input_file_fq1,opf=output.merge_fq1))
      shell("cat {ipf} > {opf}".format(ipf=input_file_fq2,opf=output.merge_fq2))
      print("cat {ipf} > {opf}".format(ipf=input_file_fq1,opf=output.merge_fq1))
      print("cat {ipf} > {opf}".format(ipf=input_file_fq2,opf=output.merge_fq2))

    else:
      input_file_fq1=" ".join(input.fq1)
      input_file_fq2=" ".join(input.fq2)
      shell("ln -s {ipf} {opf}".format(ipf=input_file_fq1,opf=output.merge_fq1))
      shell("ln -s {ipf} {opf}".format(ipf=input_file_fq2,opf=output.merge_fq2))
      print("ln -s {ipf} {opf}".format(ipf=input_file_fq1,opf=output.merge_fq1))
      print("ln -s {ipf} {opf}".format(ipf=input_file_fq2,opf=output.merge_fq2))


rule STAR_mapping:
  input:
    fq1=config['RD']+"{sample}_R1.fastq.gz",
    fq2=config['RD']+"{sample}_R2.fastq.gz",
    GTF=config['GTF'],
    GD_star_w=config['GD_star']+"sjdbList.fromGTF.out.tab",
  output:
    bam=config['RD']+"{sample}/{sample}.Aligned.out.bam",
    tbam=temp(config['RD']+"{sample}/{sample}.Aligned.toTranscriptome.out.bam"),
    mlog=config['RD']+"{sample}/{sample}.Log.final.out",
  #threads:4
  shell:
    "STAR --runThreadN 4 --outFilterMultimapNmax 1 --genomeDir $(dirname {input.GD_star_w}) --readFilesIn {input.fq1} {input.fq2} --quantMode TranscriptomeSAM --outSAMstrandField intronMotif --readFilesCommand zcat  --outSAMtype BAM Unsorted  --outFileNamePrefix $(dirname {output.bam})/$(basename $(dirname {output.bam}))."

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
    gcount=expand(config['RD']+"{sample}/{sample}.rsem.genes.results", sample=Samples.keys()),
    msrc=config['merge_src']
  output:
    mcount=config['RD']+"merge.rsem_counts.csv",
    mtpm=config['RD']+"merge.rsem_tpm.csv",
  shell:
    "python {input.msrc} -i {input.gcount} -o $(dirname {output.mcount})/merge -f .rsem.genes.results"