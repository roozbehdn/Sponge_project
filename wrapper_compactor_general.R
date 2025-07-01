library(data.table)
library(tictoc)

args <- commandArgs(trailingOnly = TRUE)
directory = args[1]
num_kmer = args[2]
kmer_len = args[3]


print(directory)
setwd(directory)
system("mkdir compactors_granny_mapping_anchors_30iteration")
system("cut -f2 -d \" \" sample_sheet.txt > compactors_granny_mapping_anchors_30iteration/compactor_samplesheet.txt")

system("cp /oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/Time_course_analysis/sig_ancors_in_normal_or_immune_mapped_to_peter_cdna.txt  compactors_granny_mapping_anchors_30iteration/sig_ancors_in_normal_or_immune_mapped_to_peter_cdna.txt")

setwd(paste(directory,"compactors_granny_mapping_anchors_30iteration",sep=""))


system(paste("sbatch -p horence,owners,quake,normal --time=48:00:00 --mem=50000 --wrap=\"/oak/stanford/groups/horence/Roozbeh/software/R-NOMAD/bin/compactors --max_length 100000 --num_kmers ", num_kmer, " --kmer_len ", kmer_len, " --min_extender_specificity 0.5 --lower_bound 1 --epsilon 0.01 --beta 0.5 compactor_samplesheet.txt sig_ancors_in_normal_or_immune_mapped_to_peter_cdna.txt sig_ancors_in_normal_or_immune_mapped_to_peter_cdna_compactors.tsv\"", sep = ""))
