#!/bin/sh
#############################
# File Name : kallisto_run.sh
#
# Purpose : [???]
#
# Creation Date : 20-02-2024
#
# Last Modified : Tue 20 Feb 2024 02:54:23 PM PST
#
# Created By : Roozbeh Dehghannasiri
#
##############################
/oak/stanford/groups/horence/Roozbeh/software/kallisto/kallisto quant -i /oak/stanford/groups/horence/Roozbeh/Sponge_project/Granny_Data/Jacob_assembly/kallisto/transcriptom_with_cDNAs_uniq -o output --plaintext --single -l 200 -s 20 /oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Sponge_RNA_immune_response_RQ23179/Sl_D12_CD-2_245_332_S19_L005_R1_001.fastq.gz /oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Sponge_RNA_immune_response_RQ23179/Sl_D12_LPS-3_269_308_S17_L005_R1_001.fastq.gz /oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Sponge_RNA_immune_response_RQ23179/Sl_D8_CD-3_257_320_S18_L005_R1_001.fastq.gz /oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Sponge_RNA_immune_response_RQ23179/Sl_D8_LPS-2_281_296_S16_L005_R1_001.fastq.gz

