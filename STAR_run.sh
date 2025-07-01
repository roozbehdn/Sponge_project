#!/bin/sh
#############################
# File Name : STAR_run.sh
#
# Purpose : [???]
#
# Creation Date : 22-02-2024
#
# Last Modified : Thu 22 Feb 2024 11:25:38 AM PST
#
# Created By : Roozbeh Dehghannasiri
#
##############################
INFILE=$1
for x in $(cat ${INFILE})
do
sbatch -p horence,owners,quake,normal --time=12:00:00 --mem=20000 --wrap="/oak/stanford/groups/horence/Roozbeh/software/STAR-2.7.5a/bin/Linux_x86_64/STAR --runThreadN 4 --genomeDir /oak/stanford/groups/horence/Roozbeh/NOMAD_10x/utility_files/non_model_organisms_references/Spongilla_lacustris/STAR_index_files --readFilesIn ${x}  --outFileNamePrefix ${x}_to_NCBI --twopassMode Basic --alignIntronMax 1000000 --limitOutSJcollapsed 3000000 --chimJunctionOverhangMin 10 --chimSegmentReadGapMax 0 --chimOutJunctionFormat 1 --chimSegmentMin 12 --chimScoreJunctionNonGTAG -4 --chimNonchimScoreDropMin 10 --outSAMtype SAM --chimOutType SeparateSAMold --outSAMunmapped None --clip3pAdapterSeq AAAAAAAAA --outSAMattributes NH HI AS nM NM"
done
