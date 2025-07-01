library(data.table)
#dir="/scratch/groups/horence/Roozbeh/non_model_organisms/Sponge_SRP216435/"
#dir="/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Sponge_10X_SRP216435/"
dir="/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Sponge_RNA_immune_response_RQ23179/"
f=list.files(dir)
f=f[which(f %like% "_R1_001.fastq.gz")]
#f=f[which(f %like% "gz")]
#f=f[which(!f %like% "granulo")]
f1=gsub(".fastq.gz","",f)
f1=gsub("_R1_001","",f1)
output_dir = "/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/alignment_to_granulocyte_markers/Sponge_RNA_immune_response_RQ23179/"

#system(paste("mkdir ",output_dir,sep=""))
for (counter in 1:length(f)){
  
  system(paste("sbatch -p horence,owners,quake,normal --time=12:00:00 --mem=10000 --wrap=\"/home/groups/horence/applications/bowtie2-2.2.1/bowtie2 -p 8 --no-unal --no-hd --no-sq -x /oak/stanford/groups/horence/Roozbeh/Sponge_project/Granny_Data/Jacob_assembly/bowtie2_granulocyte_marker_index/granulocyte_marker_bt2 -U ",dir,f[counter]," > ", output_dir, "bowtie_out_",f1[counter],".sam", "\"", sep=""))
}

