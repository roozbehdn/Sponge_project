library(data.table)
library(tictoc)

donor_tissue_sample_list = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/scripts/mouse_lemur_10x_samples.txt")
directory = "/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/10X_runs/Mouse_lemur_10x/"

for (counter in 1:nrow(donor_tissue_sample_list)){
  results_directory = paste(directory,donor_tissue_sample_list$Donor[counter],"/",donor_tissue_sample_list$Tissue[counter],sep="")
  setwd(paste(results_directory,"/compactors",sep=""))
  system("source /oak/stanford/groups/horence/george/dog/bin/activate")
  system("sbatch -p owners,quake,normal,horence --mem=20000 --time=48:00:00 --cpus-per-task=42 --wrap=\"python3 /oak/stanford/groups/horence/george/splash_utils/generic_repeat.py after_correction_anchors_compactors.tsv compactors_generic_repeat.tsv compactor\"")
  system("sbatch -p owners,quake,normal,horence --mem=20000 --time=48:00:00 --cpus-per-task=42 --wrap=\"python3 /oak/stanford/groups/horence/george/splash_utils/periodic_repeat.py after_correction_anchors_compactors.tsv compactors_periodic_repeat.tsv compactor\"")
}
