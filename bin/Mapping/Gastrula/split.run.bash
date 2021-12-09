for ITEM in {1..15}
do
  cd /home/zhaoch/snic2020-16-71/nobackup/Gastu/human_gastrula
  mkdir -p sub_${ITEM}
  for file in `\ls -d SS-*cram|grep -v sub_`
  do
    cp -rs /home/zhaoch/snic2020-16-71/nobackup/Gastu/human_gastrula/${file} sub_${ITEM}/
  done
done





cd /home/zhaoch/My_project/JP_project/bin/Gastrula
split -l 80 Sample.name.file Split.sn.file_
mv Split.sn.file_aa Split.sn.file_1
mv Split.sn.file_ab Split.sn.file_2
mv Split.sn.file_ac Split.sn.file_3
mv Split.sn.file_ad Split.sn.file_4
mv Split.sn.file_ae Split.sn.file_5
mv Split.sn.file_af Split.sn.file_6
mv Split.sn.file_ag Split.sn.file_7
mv Split.sn.file_ah Split.sn.file_8
mv Split.sn.file_ai Split.sn.file_9
mv Split.sn.file_aj Split.sn.file_10
mv Split.sn.file_ak Split.sn.file_11
mv Split.sn.file_al Split.sn.file_12
mv Split.sn.file_am Split.sn.file_13
mv Split.sn.file_an Split.sn.file_14
mv Split.sn.file_ao Split.sn.file_15

for ITEM in {1..15}
do
  cd /home/zhaoch/My_project/JP_project/bin/Gastrula
  mkdir -p sub_${ITEM}
  cp CS7.quant.map.smk sub_${ITEM}
  mv Split.sn.file_${ITEM} sub_${ITEM}
  echo "RD: /home/zhaoch/snic2020-16-71/nobackup/Gastu/human_gastrula/sub_${ITEM}/" > sub_${ITEM}/CS7.yaml
  echo "SampleName_file: /home/zhaoch/My_project/JP_project/bin/Gastrula/sub_${ITEM}/Split.sn.file_${ITEM}" >> sub_${ITEM}/CS7.yaml
  tail -n 5 CS7.yaml >> sub_${ITEM}/CS7.yaml
  head -n 6 run_slurm.sbatch.bash > sub_${ITEM}/run_slurm.sbatch.bash
  echo "cd ~/My_project/JP_project/bin/Gastrula/sub_${ITEM}" >> sub_${ITEM}/run_slurm.sbatch.bash
  tail -n 3 run_slurm.sbatch.bash >> sub_${ITEM}/run_slurm.sbatch.bash
done

