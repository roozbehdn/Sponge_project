library(data.table)


################## input parameters for the run_args.py script ######################
directory = "/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_SRP339256/"
bin_path = "/oak/stanford/groups/horence/Roozbeh/software/R-NOMAD/bin"
input_samples = "/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_SRP339256/sample_sheet.txt"
technology = "base"  # could be 10x/visium/base
cbc_len = 16
umi_len = 10
anchor_len = 27
target_len = 27
gap_len = 0
poly_ACGT_len = 6
#anchor_count_threshold = 1
#anchor_samples_threshold = 1
#anchor_sample_counts_threshold = 1
train_fraction = 0.2
cbc_filtering_thr = 3000
postprocessing_json_name = "classification.json"
#####################################################################################

system(paste("cp /oak/stanford/groups/horence/Roozbeh/software/R-NOMAD/example/run_args.py ",directory,sep=""))
system(paste("cp /oak/stanford/groups/horence/Roozbeh/software/R-NOMAD/bin/splash.py ",directory,sep=""))
setwd(directory)
print(paste("sbatch -p horence,owners,quake,normal --time=48:00:00 --mem=50000 --wrap=\"python3 ", directory, "/run_args.py -b ", bin_path, " -d ",  technology, " -c ", cbc_len, " -u ", umi_len, " -a ", anchor_len," -t ", target_len, " -g ", gap_len, " -p ", poly_ACGT_len, " -f ", train_fraction, " -e ", cbc_filtering_thr,  " -i ", input_samples, " -j ", postprocessing_json_name, "\"", sep = ""))
system(paste("sbatch -p horence,owners,quake,normal --time=48:00:00 --mem=50000 --wrap=\"python3 ", directory, "/run_args.py -b ", bin_path, " -d ", technology, " -c ", cbc_len, " -u ", umi_len, " -a ", anchor_len," -t ", target_len, " -g ", gap_len, " -p ", poly_ACGT_len, " -f ", train_fraction, " -e ", cbc_filtering_thr, " -i ", input_samples, " -j ", postprocessing_json_name, "\"", sep = ""))
