#!/bin/sh
#$ -cwd
#$ -j y
#$ -l h_vmem=10G
#$ -l h_rt=2:00:00
#$ -t 1-48
#$ -N PRSice

source /broad/software/scripts/useuse
reuse .gcc-7.3.0

i=$(expr ${SGE_TASK_ID} - 1)

chr=${chrlist[$i]}
gene=${genelist[$i]}

for r2 in 0.1 0.01 0.001; do
for cdis in 250kb 5Mb; do
if [ ! -s ${work_dir}/${outfolder}/${base}.${gene}.rsq${r2}.cd${cdis}.valid ]; then
$PRSice \
--base ${work_dir}/formated.gwas/${base}.${gene} \
--target ${tar}.chr${chr} \
--out ${work_dir}/${outfolder}/${base}.${gene}.rsq${r2}.cd${cdis} \
--ld /medpop/esp2/yruan/projects/zyu.prs/data/target/g1000_eur.allbase \
--fastscore \
--print-snp \
--memory 10Gb \
--stat BETA \
--A1 ALT \
--A2 REF \
--pvalue P \
--binary-target F \
--snp RSID \
--bar-levels 5E-8,1E-5,0.001,0.01,0.1 \
--clump-r2 ${r2} \
--clump-kb ${cdis} \
--score sum \
--missing CENTER \
--no-regress \
--no-full \
--all-score \
--thread 1

fi


if [ ! -s ${work_dir}/${outfolder}/${base}.${gene}.rsq${r2}.cd${cdis}.prsice ]; then
$PRSice \
--base ${work_dir}/formated.gwas/${base}.${gene} \
--target ${tar}.chr${chr} \
--out ${work_dir}/${outfolder}/${base}.${gene}.rsq${r2}.cd${cdis} \
--ld /stanley/huang_lab/home/yruan/DATA/Raw_Data/g1000_${base_pop} \
--fastscore \
--print-snp \
--memory 10Gb \
--stat BETA \
--A1 ALT \
--A2 REF \
--pvalue P \
--binary-target F \
--snp RSID \
--bar-levels 5E-8,1E-5,0.001,0.01,0.1 \
--clump-r2 ${r2} \
--clump-kb ${cdis} \
--score sum \
--missing CENTER \
--no-regress \
--no-full \
--all-score \
--thread 1 \
--extract ${work_dir}/${outfolder}/${base}.${gene}.rsq${r2}.cd${cdis}.valid
fi
done; done;
