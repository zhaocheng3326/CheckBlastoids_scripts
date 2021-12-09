#!/bin/bash -l
 
#SBATCH -A snic2020-15-102
#SBATCH -p core
#SBATCH -n 19
#SBATCH -t 5-00:00:00
#SBATCH -J LW61
export PATH=/home/zhaoch/snic2020-16-71/nobackup/cellranger-3.1.0:$PATH
cd /home/zhaoch/snic2020-16-71/nobackup/EBlastoids/out
cellranger count --id=LW61 --fastqs=/home/zhaoch/snic2020-16-71/nobackup/EBlastoids/LW61 --sample=LW61  --transcriptome=/home/zhaoch/snic2020-16-71/nobackup/Genome/refdata-cellranger-GRCh38-3.0.0 
