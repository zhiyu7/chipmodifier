#!/bin/sh
#$ -cwd
#$ -j y
#$ -l h_vmem=10G
#$ -l h_rt=2:00:00=
#$ -N prs_cs


source /broad/software/scripts/useuse

reuse Anaconda

i=$(expr ${SGE_TASK_ID} - 1 )

PRScs=/medpop/esp2/yruan/tools/PRScs/PRScs.py
ref_dir=/medpop/esp2/yruan/raw.data/ld.ref.prs-csx.1kg/ldblk_1kg_${base_pop}
work_dir=/medpop/esp2/yruan/projects/zyu.prs/data
bim_prefix=/medpop/esp2/yruan/projects/zyu.prs/data/target/ukbb.zyu.testing.chr${chr}

# require phi and base

python $PRScs \
 --ref_dir=${ref_dir} \
 --bim_prefix=${bim_prefix} \
 --sst_file=${work_dir}/formated.gwas/${base}.${gene} \
 --n_gwas=${n_trn} \
 --chrom=$chr \
 --phi=$phi \
 --out_dir=${work_dir}/posterior/PRS-CS.${base}.${gene}
