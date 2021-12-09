# -*- snakemake -*-
#import os.path
#mkdir -p ${RD}/${SA}
BAMS=dict()
CELLS=dict()
with open(config['SampleName_file'],"r") as SampleName_file:
    for line in SampleName_file:
        bam=line.strip("\n").split("\t")[0]
        if bam == "raw_file":
            continue
        cell=line.strip("\n").split("\t")[6]
        em=line.strip("\n").split("\t")[4]
        BAMS[cell]=BAMS.get(cell,list())
        BAMS[cell].append(bam)
        CELLS[em]=CELLS.get(em,list())
        if cell not in CELLS[em]:
            CELLS[em].append(cell)


rule merge_files:
    input:
        lambda wildcards:[config['RD']+"/" + bam for bam in BAMS[wildcards.cell]],
    output:
        merge_bam=temp(config['RD']+"/Allele/{cell}.merge.bam"),
    run:
        if (len(input) >1):
            input_file=" -I=".join(input)
            shell("gatk MergeSamFiles -I={ipf} -O {opf}".format(ipf=input_file,opf=output.merge_bam))
            print("gatk MergeSamFiles -I={ipf} -O {opf}".format(ipf=input_file,opf=output.merge_bam))
        else:
            input_file=" ".join(input)
            shell("ln -s {ipf} {opf}".format(ipf=input_file,opf=output.merge_bam))
            print("ln -s {ipf} {opf}".format(ipf=input_file,opf=output.merge_bam))


rule removeERCCandsort:
    input:
        bam=config['RD']+"/Allele/{cell}.merge.bam",
    output:
        sam=temp(config['RD']+"/Allele/{cell}.merge.sam"),
        tpbam=temp(config['RD']+"/Allele/{cell}.merge.remove.ERCC.sam"),
        bam=temp(config['RD']+"/Allele/{cell}.merge.sort.bam")
    shell:
        "samtools view {input.bam} -h |grep -v 'ERCC-' > {output.sam};samtools view -Sb {output.sam} > {output.tpbam};samtools sort {output.tpbam}  -o {output.bam} -@ 2"

rule AddOrReplaceReadGroups:
    input:
        bam=config['RD']+"/Allele/{cell}.merge.sort.bam"
    output:
        bam=temp(config['RD']+"/Allele/{cell}.merge.sort.addg.bam")
    shell:
        "gatk AddOrReplaceReadGroups -I {input.bam} -O {output.bam} -ID {wildcards.cell} -LB {wildcards.cell} -PL illumina -PU {wildcards.cell} -SM {wildcards.cell}"

rule MarkDuplicates:
    input:
        bam=config['RD']+"/Allele/{cell}.merge.sort.addg.bam"
    output:
        bam=temp(config['RD']+"/Allele/{cell}.merge.sort.dedup.bam"),
        metr=temp(config['RD']+"/Allele/{cell}.merge.sort.dedup.metrics.txt")
    shell:
        "gatk MarkDuplicates --INPUT {input.bam} --OUTPUT {output.bam}  --METRICS_FILE {output.metr} --CREATE_INDEX true  --VALIDATION_STRINGENCY SILENT "

rule SplitNCigarReads:
    input:
        bam=config['RD']+"/Allele/{cell}.merge.sort.dedup.bam",
        ref=config['REF_fa']
    output:
        bam=temp(config['RD']+"/Allele/{cell}.merge.sort.dedup.splitN.bam"),
    shell:
        "gatk SplitNCigarReads -R {input.ref} -I {input.bam}  -O  {output.bam}"

rule BaseRecalibrator:
    input:
        bam=config['RD']+"/Allele/{cell}.merge.sort.dedup.splitN.bam",
        ref=config['REF_fa'],
        DBSNP=config['DBSNP']
    output:
        rept=temp(config['RD']+"/Allele/{cell}.merge.sort.dedup.splitN.BaseR.report"),
    shell:
        "gatk BaseRecalibrator -R {input.ref}  -I {input.bam} --use-original-qualities -O {output.rept} -known-sites  {input.DBSNP} "

rule ApplyBQSR:
    input:
        bam=config['RD']+"/Allele/{cell}.merge.sort.dedup.splitN.bam",
        ref=config['REF_fa'],
        rept=temp(config['RD']+"/Allele/{cell}.merge.sort.dedup.splitN.BaseR.report"),
    output:
        bam=temp(config['RD']+"/Allele/{cell}.merge.sort.dedup.splitN.BaseR.BSSR.bam"),
    shell:
        "gatk ApplyBQSR --add-output-sam-program-record -R {input.ref} -I {input.bam}  --use-original-qualities -O {output.bam} --bqsr-recal-file  {input.rept}"

rule extra_sort:
    input:
        bam=config['RD']+"/Allele/{cell}.merge.sort.dedup.splitN.BaseR.BSSR.bam",
        ref=config['REF_fa'],
    output:
        bam=temp(config['RD']+"/Allele/{cell}.merge.sort.dedup.splitN.BaseR.BSSR.es.bam"),
    shell:
        "gatk ReorderSam  -I {input.bam}  -R {input.ref} -O {output.bam} "

rule merge_bam_by_EM:
    input:
        lambda wildcards:[config['RD']+"/Allele/" + cell +".merge.sort.dedup.splitN.BaseR.BSSR.es.bam" for cell in CELLS[wildcards.em]],
    output:
        merge_bam=config['RD']+"/Allele/{em}.merge.bam",
    run:
        if (len(input) >1):
            input_file=" -I=".join(input)
            shell("gatk MergeSamFiles -I={ipf} -O {opf}".format(ipf=input_file,opf=output.merge_bam))
            print("gatk MergeSamFiles -I={ipf} -O {opf}".format(ipf=input_file,opf=output.merge_bam))
        else:
            input_file=" ".join(input)
            shell("cat {ipf} > {opf}".format(ipf=input_file,opf=output.merge_bam))
            print("cat {ipf} > {opf}".format(ipf=input_file,opf=output.merge_bam))


rule HaplotypeCaller:
    input:
        bam=config['RD']+"/Allele/{em}.merge.bam",
        ref=config['REF_fa'],
        DBSNP=config['DBSNP']
    output:
        vcf=config['RD']+"/Allele/{em}.merge.vcf",
        bai=config['RD']+"/Allele/{em}.merge.bam.bai",
    shell:
        "samtools index {input.bam}; gatk HaplotypeCaller -R {input.ref}  -I {input.bam} -O {output.vcf} --dbsnp {input.DBSNP} -dont-use-soft-clipped-bases"

rule shortVCF:
    input:
        vcf=config['RD']+"/Allele/{em}.merge.vcf",
        vsrc=config['VCF_src']
    output:
        vcf=config['RD']+"/Allele/{em}.merge.ana.vcf",
    shell:
        "python3 {input.vsrc} {input.vcf} > {output.vcf}"

rule temp_run:
    input:
        expand(config['RD']+"/Allele/{em}.merge.ana.vcf", em=CELLS.keys()),

