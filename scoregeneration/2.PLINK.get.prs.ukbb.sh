#!/bin/sh
#$ -cwd
#$ -j y
#$ -l h_vmem=4G
#$ -N plink.prs


plink2=/stanley/huang_lab/home/yruan/Software/plink2


dat=${base}.${gene}

for phi in ${philist[@]}; do
$plink2 \
--bfile //medpop/esp2/yruan/projects/zyu.prs/data/target/ukbb.zyu.testing.chr${chr} \
--memory 4000 \
--score ${work_dir}/posterior/PRS-CS.${dat}_pst_eff_a1_b0.5_phi${phi}_chr${chr}.txt 2 4 6 \
--out  ${work_dir}/prs/${dat}.phi${phi}

done
