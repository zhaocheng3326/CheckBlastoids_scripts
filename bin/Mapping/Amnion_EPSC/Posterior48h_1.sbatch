#!/bin/bash -l
 
#SBATCH -A snic2020-15-102
#SBATCH -p core
#SBATCH -n 19
#SBATCH -t 5-00:00:00
#SBATCH -J MQ_H9AM
export PATH=/home/zhaoch/snic2020-16-71/nobackup/cellranger-3.1.0:$PATH
cd /home/zhaoch/snic2020-16-71/nobackup/Amnion_EPSC/out
cellranger count --id=Posterior48h_1 --fastqs=/home/zhaoch/snic2020-16-71/nobackup/Amnion_EPSC/Posterior48h_1 --sample=Posterior48h_1 --transcriptome=/home/zhaoch/snic2020-16-71/nobackup/Genome/refdata-cellranger-GRCh38-3.0.0 
