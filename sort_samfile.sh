#!/bin/sh
#############################
# File Name : sort_samfile.sh
#
# Purpose : [???]
#
# Creation Date : 05-04-2024
#
# Last Modified : Mon 08 Apr 2024 02:38:23 PM PDT
#
# Created By : Roozbeh Dehghannasiri
#
##############################

cut -f2 sig_ancors_in_normal_or_immune_mapped_to_peter_cdna_compactors.tsv | grep -v "comp" | makefasta > compactors.fasta

/oak/stanford/groups/horence/Roozbeh/software/minimap2/minimap2 -N 0 -a /oak/stanford/groups/horence/Roozbeh/Sponge_project/Granny_Data/Jacob_assembly/Peters_cDNAs_11-20.fasta compactors.fasta > compactors_minimap.sam


samtools view -bS compactors_minimap.sam > compactors_minimap.bam
samtools sort compactors_minimap.bam  compactors_minimap_sorted
samtools index compactors_minimap_sorted.bam
