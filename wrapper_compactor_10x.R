library(data.table)
library(tictoc)

donor_tissue_sample_list = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/scripts/test.txt")
directory = "/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/10X_runs/Mouse_lemur_10x/"

for (counter in 1:nrow(donor_tissue_sample_list)){
  results_directory = paste(directory,donor_tissue_sample_list$Donor[counter],"/",donor_tissue_sample_list$Tissue[counter],sep="")
  print(results_directory)
  setwd(results_directory)
  system("mkdir compactors")
  system(paste("cut -f 2 -d \",\" /oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/10X_runs/Mouse_lemur_10x/input_sample_files/",donor_tissue_sample_list$Donor[counter],"*",donor_tissue_sample_list$Tissue[counter],"*S*.txt > compactors/compactor_samplesheet.txt",sep=""))
 # system("awk ' { print \"/scratch/groups/horence/Roozbeh/single_cell_project/data/TSP21_10X/\" $0 } ' compactors/compactor_samplesheet1.txt > compactors/compactor_samplesheet.txt")
 # system("rm compactors/compactor_samplesheet1.txt")
#  system("rm compactors/compactor_samplesheet2.txt")
  system("cut -f1  result.after_correction.all_anchors.tsv | grep -v \"anch\" > compactors/before_correction_anchors.txt")
  setwd(paste(results_directory,"/compactors",sep=""))
  system("sbatch -p horence,owners,quake --mem=100000 --time=48:00:00 --wrap=\"/oak/stanford/groups/horence/Roozbeh/software/R-NOMAD/bin/compactors --num_kmers 2 compactor_samplesheet.txt before_correction_anchors.txt before_correction_anchors_compactors.tsv\"")
  
}
