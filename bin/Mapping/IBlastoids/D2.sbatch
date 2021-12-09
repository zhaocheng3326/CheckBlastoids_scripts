#!/bin/bash -l
 
#SBATCH -A snic2020-15-102
#SBATCH -p core
#SBATCH -n 19
#SBATCH -t 5-00:00:00
#SBATCH -J D2
export PATH=/home/zhaoch/snic2020-16-71/nobackup/cellranger-3.1.0:$PATH
cd /home/zhaoch/snic2020-16-71/nobackup/IBlastoids/out
cellranger count --id=IBD2 --fastqs=/home/zhaoch/snic2020-16-71/nobackup/IBlastoids/RL2157_2020_06_24_20_2175_32F38F1_scRNAseq --sample=RL2157_2020_06_24_20_2175_32F38F1_scRNAseq  --transcriptome=/home/zhaoch/snic2020-16-71/nobackup/Genome/refdata-cellranger-GRCh38-3.0.0 
