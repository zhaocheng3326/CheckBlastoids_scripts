for ITEM in E3.1 E3.2 E3.3 E3.4 E3.45 E3.46 E3.47 E3.48 E3.49 E3.50 E3.51 E3.52 E3.53 E4.1 E4.10 E4.11 E4.2 E4.24_5_0_4 E4.3 E4.31_5_1 E4.4 E4.5 E4.6 E4.7 E4.8 E4.9 E4.late.33 E4.late.34 E4.late.35 E5.1 E5.10 E5.11 E5.12 E5.13 E5.14 E5.15 E5.16 E5.2 E5.3 E5.37 E5.38 E5.39 E5.40 E5.41 E5.42 E5.43 E5.5 E5.6 E5.7 E5.8 E5.9 E5.early.31 E5.early.36 E6.1 E6.10 E6.11 E6.12 E6.13 E6.14 E6.15 E6.16 E6.17 E6.18 E6.2 E6.22 E6.3 E6.4 E6.6 E6.7 E6.8 E6.9 E7.10 E7.11 E7.12 E7.13 E7.14 E7.15 E7.16 E7.17 E7.19 E7.2 E7.3 E7.4 E7.5 E7.6 E7.7 E7.8 E7.9
do
  mkdir -p /home/zhaoch/My_project/JP_project/bin/CP_data/allele/sub_yaml/${ITEM}
  cd /home/zhaoch/My_project/JP_project/bin/CP_data/allele/sub_yaml/${ITEM}
  #mv ../../sample.${ITEM}.meta.tsv ./
  cat /home/zhaoch/My_project/JP_project/bin/CP_data/allele/cp.map.allele.smk  > running.smk
  echo "RD: /home/zhaoch/snic2020-16-71/nobackup/CP_data" > running.yaml
  echo "SampleName_file: /home/zhaoch/My_project/JP_project/bin/CP_data/allele/sub_yaml/${ITEM}/sample.${ITEM}.meta.tsv" >> running.yaml
  echo "GTF: /home/zhaoch/snic2020-16-71/nobackup/Genome/GC_refdata-GRCh38/genes-ERCC.gtf" >> running.yaml
  echo "GD_star: /home/zhaoch/snic2020-16-71/nobackup/Genome/GC_refdata-GRCh38/genome.fa_star_index" >> running.yaml
  echo "REF_fa: /home/zhaoch/snic2020-16-71/nobackup/Genome/GC_refdata-GRCh38/genome.fa" >> running.yaml
  echo "DBSNP: /home/zhaoch/snic2020-16-71/nobackup/Genome/GC_refdata-GRCh38/dbSNP.valid.vcf" >> running.yaml
  echo "VCF_src: /home/zhaoch/snic2020-16-71/nobackup/Genome/resource/allel_analysis.py" >> running.yaml

  echo '#!/bin/bash -l' > ${ITEM}.run_slurm.sbatch.bash
  echo "#SBATCH -A snic2020-15-102" >> ${ITEM}.run_slurm.sbatch.bash
  echo "#SBATCH -p core" >> ${ITEM}.run_slurm.sbatch.bash
  echo "#SBATCH -n 18" >> ${ITEM}.run_slurm.sbatch.bash
  echo "#SBATCH -t 5-00:00:00" >> ${ITEM}.run_slurm.sbatch.bash
  echo "#SBATCH -J MAL_CP" >> ${ITEM}.run_slurm.sbatch.bash
  echo "cd /home/zhaoch/My_project/JP_project/bin/CP_data/allele/sub_yaml/${ITEM}" >> ${ITEM}.run_slurm.sbatch.bash
  echo "source ~/miniconda3/bin/activate GATK4" >> ${ITEM}.run_slurm.sbatch.bash
  echo "snakemake  temp_run -k  -s running.smk --configfile running.yaml -p -j 4" >> ${ITEM}.run_slurm.sbatch.bash
done

for ITEM in E3.2 E3.3 E3.4 E3.45 E3.46 E3.47 E3.48 E3.49 E3.50 E3.51 E3.52 E3.53 E4.1 E4.10 E4.11 E4.2 E4.24_5_0_4 E4.3 E4.31_5_1 E4.4 E4.5 E4.6 E4.7 E4.8 E4.9 E4.late.33 E4.late.34 E4.late.35 E5.1 E5.10 E5.11 E5.12 E5.13 E5.14 E5.15 E5.16 E5.2 E5.3 E5.37 E5.38 E5.39 E5.40 E5.41 E5.42 E5.43 E5.5 E5.6 E5.7 E5.8 E5.9 E5.early.31 E5.early.36 E6.1 E6.10 E6.11 E6.12 E6.13 E6.14 E6.15 E6.16 E6.17 E6.18 E6.2 E6.22 E6.3 E6.4 E6.6 E6.7 E6.9 E7.10 E7.11 E7.12 E7.13 E7.14 E7.15 E7.16 E7.17 E7.19 E7.2 E7.3 E7.4 E7.5 E7.6 E7.7 E7.8 E7.9
do
  cd /home/zhaoch/My_project/JP_project/bin/CP_data/allele/sub_yaml/${ITEM}
  sbatch ${ITEM}.run_slurm.sbatch.bash
done	


python ~/My_project/ScRNA_analysis/src/vcf_merge_v2.py E10_EM1.merge.ana.vcf E10_EM2.merge.ana.vcf E10_EM3.merge.ana.vcf E10_EM4.merge.ana.vcf E10_EM5.merge.ana.vcf E3_EM10.merge.ana.vcf E3_EM11.merge.ana.vcf E3_EM1.merge.ana.vcf E3_EM2.merge.ana.vcf E3_EM3.merge.ana.vcf E3_EM4.merge.ana.vcf E3_EM5.merge.ana.vcf E3_EM6.merge.ana.vcf E3_EM7.merge.ana.vcf E3_EM8.merge.ana.vcf E3_EM9.merge.ana.vcf E4_EM1.merge.ana.vcf E4_EM2.merge.ana.vcf E4_EM3.merge.ana.vcf E4_EM4.merge.ana.vcf E4_EM5.merge.ana.vcf E4_EM6.merge.ana.vcf E4_EM7.merge.ana.vcf E4_EM8.merge.ana.vcf E4_EM9.merge.ana.vcf E5_EM10.merge.ana.vcf E5_EM11.merge.ana.vcf E5_EM12.merge.ana.vcf E5_EM13.merge.ana.vcf E5_EM14.merge.ana.vcf E5_EM15.merge.ana.vcf E5_EM16.merge.ana.vcf E5_EM17.merge.ana.vcf E5_EM18.merge.ana.vcf E5_EM19.merge.ana.vcf E5_EM1.merge.ana.vcf E5_EM20.merge.ana.vcf E5_EM21.merge.ana.vcf E5_EM22.merge.ana.vcf E5_EM23.merge.ana.vcf E5_EM24.merge.ana.vcf E5_EM25.merge.ana.vcf E5_EM26.merge.ana.vcf E5_EM27.merge.ana.vcf E5_EM2.merge.ana.vcf E5_EM3.merge.ana.vcf E5_EM4.merge.ana.vcf E5_EM5.merge.ana.vcf E5_EM6.merge.ana.vcf E5_EM7.merge.ana.vcf E5_EM8.merge.ana.vcf E5_EM9.merge.ana.vcf E6_EM10.merge.ana.vcf E6_EM11.merge.ana.vcf E6_EM12.merge.ana.vcf E6_EM13.merge.ana.vcf E6_EM14.merge.ana.vcf E6_EM15.merge.ana.vcf E6_EM16.merge.ana.vcf E6_EM17.merge.ana.vcf E6_EM18.merge.ana.vcf E6_EM19.merge.ana.vcf E6_EM1.merge.ana.vcf E6_EM20.merge.ana.vcf E6_EM21.merge.ana.vcf E6_EM22.merge.ana.vcf E6_EM23.merge.ana.vcf E6_EM24.merge.ana.vcf E6_EM25.merge.ana.vcf E6_EM2.merge.ana.vcf E6_EM3.merge.ana.vcf E6_EM4.merge.ana.vcf E6_EM5.merge.ana.vcf E6_EM6.merge.ana.vcf E6_EM7.merge.ana.vcf E6_EM8.merge.ana.vcf E6_EM9.merge.ana.vcf E7_EM10.merge.ana.vcf E7_EM11.merge.ana.vcf E7_EM12.merge.ana.vcf E7_EM13.merge.ana.vcf E7_EM1.merge.ana.vcf E7_EM2.merge.ana.vcf E7_EM3.merge.ana.vcf E7_EM4.merge.ana.vcf E7_EM5.merge.ana.vcf E7_EM6.merge.ana.vcf E7_EM7.merge.ana.vcf E7_EM8.merge.ana.vcf E7_EM9.merge.ana.vcf | sed -e 's/#CHROM/CHROM/'  >  merge.ana.vcf 
