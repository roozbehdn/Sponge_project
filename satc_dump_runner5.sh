#!/bin/sh
#############################
# File Name : satc_dump_runner.sh
#
# Purpose : [???]
#
# Creation Date : 09-02-2024
#
# Last Modified : Wed 22 May 2024 03:30:42 PM PDT
#
# Created By : Roozbeh Dehghannasiri
#
##############################

anchLst=/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/Time_course_analysis/sig_ancors_in_normal_or_immune_mapped_to_peter_cdna.txt
outFldr=/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/Time_course_analysis/10X_counts/
satcFldr=/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/10X_runs/Sponge_SRP216435_all_samples_together/result_satc/

for i in {0..127}
do
  #  satcFile=${satcFldr}/result.bin${i}.satc
    satcFile=${satcFldr}/bin${i}.satc
    outFile=${outFldr}/satcOut${i}.tsv
    cmd="/oak/stanford/groups/horence/Roozbeh/software/R-NOMAD/bin/satc_dump --anchor_list ${anchLst} ${satcFldr} ${outFldr}"
#     echo $cmd
    /oak/stanford/groups/horence/Roozbeh/software/R-NOMAD/bin/satc_dump --anchor_list ${anchLst} ${satcFile} ${outFile}
done
