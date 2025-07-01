library(data.table)

donor_tissue_sample_list = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/scripts/test.txt")
################## input parameters for the run_args.py script ######################
directory = "/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/10X_runs/Mouse_lemur_10x/"
bin_path = "/oak/stanford/groups/horence/Roozbeh/software/R-NOMAD/bin"
technology = "10x" # could be either of base/10x/visium
cbc_len = 16
umi_len = 10
anchor_len = 27
target_len = 27
gap_len = 0
poly_ACGT_len = 6
train_fraction = 0.2
cbc_filtering_thr = 2000
#####################################################################################

for (counter in 1:nrow(donor_tissue_sample_list)){
  results_directory = paste(directory,donor_tissue_sample_list$Donor[counter],"/",donor_tissue_sample_list$Tissue[counter],sep="")
  system(paste("mkdir ", directory, donor_tissue_sample_list$Donor[counter], sep = ""))
  system(paste("mkdir ", results_directory, sep = ""))
  input_samples = donor_tissue_sample_list$Sample_list[counter]
  system(paste("cp /oak/stanford/groups/horence/Roozbeh/software/R-NOMAD/example/run_args.py ",results_directory,sep=""))
  system(paste("cp /oak/stanford/groups/horence/Roozbeh/software/R-NOMAD/bin/splash.py ",results_directory,sep=""))
  setwd(results_directory)
  print(paste("sbatch -p owners,quake,normal --time=24:00:00 --mem=50000 --wrap=\"python3 ", results_directory, "/run_args.py -b ", bin_path, " -d ",  technology, " -c ", cbc_len, " -u ", umi_len, " -a ", anchor_len, " -t ", target_len, " -g ", gap_len, " -p ", poly_ACGT_len, " -f ", train_fraction, " -e ", cbc_filtering_thr,  " -i ", input_samples, "\"", sep = ""))
  system(paste("sbatch -p owners,quake,normal --time=24:00:00 --mem=50000 --wrap=\"python3 ", results_directory, "/run_args.py -b ", bin_path, " -d ", technology, " -c ", cbc_len, " -u ", umi_len, " -a ", anchor_len, " -t ", target_len, " -g ", gap_len, " -p ", poly_ACGT_len, " -f ", train_fraction, " -e ", cbc_filtering_thr, " -i ", input_samples, "\"", sep = ""))
}
