library(data.table)
library(stringdist)
library(RColorBrewer)
library(Hmisc)

#counts = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/Time_Course_runs/Sponge_RNA_immune_response_RQ23179/satc_sig_anchor/control_anchors.txt")
#names(counts) = c("sample","anchor","target","count")

counts = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/sig_RNA_developmental_anchor_counts_granny_mapping_new.tsv")



final_dt = data.table() # this data table includes the number of true targets estimated for each dataset
p_seq_err = 0.01
counts = counts[dataset!="bulk_10X_matched" & !dataset%like%"DNA" ] #dataset!="bulk_10X_matched" & !dataset%like%"DNA"
counts[,anchor_counts_per_dataset:=sum(count),by=paste(anchor,dataset)]
counts[,min_anchor_count_per_dataset:=min(anchor_counts_per_dataset)] # this gives the minimum anchor count per dataset across all datasets I use this to downsample the anchor count in datasets that have higher counts
counts[,target_count_per_dataset:=sum(count),by=paste(anchor,target,dataset)]
list.dataset = unique(counts$dataset)
for (counter_dataset in 1:length(list.dataset)){
  print(counter_dataset)
  counts_dataset = counts[dataset==list.dataset[counter_dataset]]
  counts_dataset[,anchor_count:=sum(count),by=list(anchor)]
  counts_dataset[,target_count:=sum(count),by=list(anchor, target)]
  counts_dataset = unique(counts_dataset[,list(target_count,anchor,target, anchor_count, min_anchor_count_per_dataset)])
  
  counts_dataset[,target_frac:=target_count/anchor_count]
  counts_dataset[,target_order:=rank(-target_frac,ties.method = "random"),by=list(anchor)] #ordering targets within each anchor 
  
  list.anchor = unique(counts_dataset$anchor)
  for (counter_anchor in seq(1,length(list.anchor),1)){
    print(counter_anchor)
    counts_anchor = counts_dataset[anchor==list.anchor[counter_anchor]]
    counts_anchor = counts_anchor[order(-target_frac)]
    
    for (counter_iter in 1:20){
      counts_anchor_downsampled = counts_anchor[sample.int(nrow(counts_anchor), unique(counts_anchor$min_anchor_count_per_dataset), replace = TRUE, prob = counts_anchor$target_frac),]
      counts_anchor_downsampled[,anchor_count:=.N,by=list(anchor)]
      counts_anchor_downsampled[,target_count:=.N,by=list(anchor, target)]
      counts_anchor_downsampled = unique(counts_anchor_downsampled[,list(target_count,anchor,target, anchor_count,target_order)])
      counts_anchor_downsampled[,target_frac:=target_count/anchor_count]
      
      counts_anchor_downsampled[,target_order:=rank(-target_frac,ties.method = "random"),by=list(anchor)] #ordering targets within each anchor 
      counts_anchor_downsampled = counts_anchor_downsampled[order(-target_frac)]
    #  counts_anchor_downsampled = copy(counts_anchor)
      counts_anchor_downsampled[,target_frac_99_CI_lower_bound:=binconf(target_count,anchor_count,0.05)[2],by=1:nrow(counts_anchor_downsampled)]
      counts_anchor_downsampled$true_target = 0
      counts_anchor_downsampled$min_HD = 0
      counts_anchor_downsampled$null_error_prob = 0
      counts_anchor_downsampled$true_target[1] = 1
      counts_anchor_downsampled = counts_anchor_downsampled[target_count>1]
      
      if (nrow(counts_anchor_downsampled)>1){
        for (counter_target in 2:nrow(counts_anchor_downsampled)){
          counts_anchor_downsampled$min_HD[counter_target] = min(stringdist(counts_anchor_downsampled$target[counter_target],  counts_anchor_downsampled[true_target==1]$target)) # compute minimum hamming distance between the target of interest and the set of true targets
          counts_anchor_downsampled$null_error_prob[counter_target] = p_seq_err^counts_anchor_downsampled$min_HD[counter_target]*sum(counts_anchor_downsampled[true_target==1]$target_frac)
          if (counts_anchor_downsampled$null_error_prob[counter_target] < counts_anchor_downsampled$target_frac_99_CI_lower_bound[counter_target]){
            counts_anchor_downsampled$true_target[counter_target]=1
          }
        }
      }
      
      final_dt = rbind(final_dt, data.table(list.anchor[counter_anchor],list.dataset[counter_dataset],sum(counts_anchor_downsampled$true_target)))
    }
    
  }
}

names(final_dt) = c("anchor","dataset","num_true_targets")

## now below I want to classify anchors as granny, HD with granny, N_terminus, C_terminus, ....
granny_repeat = "CCAGCCATCAGAACCCCAGGAACCATCTAACCAGCCATCAGAACCCCAGGAACCATCTAACCAGCCATCAGAACCCCAGGAACCATCTAA" # 3 consecutive granny repeats
granny_substrings <- unique(sapply(1:(nchar(granny_repeat) - 26), function(i) substr(granny_repeat, i, i + 26)))

final_dt[anchor%in%granny_substrings,anchor_type:="Granny"]
final_dt[,granny_HD:=min(stringdist(anchor, granny_substrings)), by=anchor]
final_dt[granny_HD<7, anchor_type:=paste("Granny: HD",granny_HD,sep="")]

lysine_rich_cdna2 = "CCAGCCAGTGGAACAGCCTGACCACAGTGGTCCTGGATCTGGACCTAAAAATCATATAGGGAAGAAAAGACCAAAAAATACAAAGCCAAAAAATAAAAATAATAAACAACCTAAAGATAAAAAAGGTCCCAAGAAAGATAAAAAACCACACAAAAAGCCCATCAAAGACCCAAGAAAAAAGCCCATCCATCATGATCATCAAGGTGGTAGCCGTAGAAATGGGGATCGTAGAGGAGGTGGTCGACACAGAGGTGATGGTGGTAATCGAGGTGGTCAAA"
lysine_rich_cdna1 = "CCAGCCAGTGGAACAGCCTGACCGCAGTTGTCCTGGAAAGAAAAGACCAAAAAAAATAAAGCCAAAATATAAAAAAAACAAAGATATAAATGGTCACATGAAATGTAAAAAAACCCACAAAATAACCATCAAACACCCAAGAAAAAAGCCCATCCATCATAATCATCAAGGTGGTAGCCGTAGAAATGGGGATCGTAGAGGAGGTGGTCAACACAGAGGTGATGCTGGTAATCGAGGTGGTCAAATAGATGGTGGTCATCGAATTGGTGGTGGTCAAATTGGTGGTGGTCGTAGAATTAGCTATCACAGTGATGGTAGACGTCGTTATGGTTGA"
lysine_substrings_rich1 <- sapply(1:(nchar(lysine_rich_cdna1) - 26), function(i) substr(lysine_rich_cdna1, i, i + 26))
lysine_substrings_rich2 <- sapply(1:(nchar(lysine_rich_cdna2) - 26), function(i) substr(lysine_rich_cdna2, i, i + 26))
lysine_substrings = unique(c(lysine_substrings_rich1,lysine_substrings_rich2))
final_dt[anchor%in%lysine_substrings,anchor_type:="lysine"]
final_dt[anchor%like%"GGTGGTCGTCGAT",anchor_type:="C_terminus_repeat"]

final_dt[,average_num_target:=mean(num_true_targets),by=anchor]
final_dt = final_dt[order(-average_num_target)]
final_dt$anchor_type = factor(final_dt$anchor_type,levels=c("Granny: HD0", "Granny: HD1", "Granny: HD2", "Granny: HD3", "Granny: HD4", "Granny: HD5", "Granny: HD6", "C_terminus_repeat", "lysine"))
final_dt = setorder(final_dt, anchor_type ,-average_num_target)
final_dt$anchor = factor(final_dt$anchor, levels = unique(final_dt$anchor))
final_dt_orig = copy(final_dt)
final_dt[,average_num_target_per_anchor_dataset:=mean(num_true_targets),by=paste(anchor,dataset)]
final_dt = final_dt[!duplicated(paste(anchor,dataset))]
ggplot(final_dt, aes(x=anchor, y=average_num_target_per_anchor_dataset,  color=anchor_type,shape = dataset)) +     geom_point() + theme_bw()  +scale_color_manual(values=brewer.pal(length(unique(final_dt$anchor_type)), 'Set1')) +scale_shape_manual(values= c(1,3,8,17,6,9)) + ylim(0,15)# 6,9



compactors_RNA_development = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/Time_Course_runs/Sponge_RNA_developmental_timecourse_RQ23078_RC/compactors_granny_mapping_anchors_4iteration/sig_ancors_in_normal_or_immune_mapped_to_peter_cdna_compactors.tsv")
compactors_RNA_immune = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/Time_Course_runs/Sponge_RNA_immune_response_RQ23179_RC/compactors_granny_mapping_anchors_4iteration/sig_ancors_in_normal_or_immune_mapped_to_peter_cdna_compactors.tsv")
compactors_DNA_development = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/Time_Course_runs/Sponge_DNA_developmental_timecourse_RQ23077/compactors_granny_mapping_anchors_4iteration/sig_ancors_in_normal_or_immune_mapped_to_peter_cdna_compactors.tsv")
compactors_DNA_immune = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/Time_Course_runs/Sponge_DNA_immune_response_RQ23180/compactors_granny_mapping_anchors_4iteration/sig_ancors_in_normal_or_immune_mapped_to_peter_cdna_compactors.tsv")

compactors_RNA_development = compactors_RNA_development[num_extended==0]
compactors_RNA_immune = compactors_RNA_immune[num_extended==0]
compactors_DNA_development = compactors_DNA_development[num_extended==0]
compactors_DNA_immune = compactors_DNA_immune[num_extended==0]

compactors_RNA_development[,num_comp_per_anchor:=length(unique(compactor)),by=anchor]
compactors_RNA_immune[,num_comp_per_anchor:=length(unique(compactor)),by=anchor]
compactors_DNA_development[,num_comp_per_anchor:=length(unique(compactor)),by=anchor]
compactors_DNA_immune[,num_comp_per_anchor:=length(unique(compactor)),by=anchor]

final_dt = merge(final_dt,compactors_RNA_development[!duplicated(anchor),list(anchor,num_comp_per_anchor)], all.x=T, all.y=F, by.x="anchor", by.y="anchor")
final_dt[is.na(num_comp_per_anchor),num_comp_per_anchor:=0]
setnames(final_dt,"num_comp_per_anchor","num_comp_RNA_development")
final_dt = merge(final_dt,compactors_RNA_immune[!duplicated(anchor),list(anchor,num_comp_per_anchor)], all.x=T, all.y=F, by.x="anchor", by.y="anchor")
final_dt[is.na(num_comp_per_anchor),num_comp_per_anchor:=0]
setnames(final_dt,"num_comp_per_anchor","num_comp_RNA_immune")
final_dt = merge(final_dt,compactors_DNA_development[!duplicated(anchor),list(anchor,num_comp_per_anchor)], all.x=T, all.y=F, by.x="anchor", by.y="anchor")
final_dt[is.na(num_comp_per_anchor),num_comp_per_anchor:=0]
setnames(final_dt,"num_comp_per_anchor","num_comp_DNA_development")
final_dt = merge(final_dt,compactors_DNA_immune[!duplicated(anchor),list(anchor,num_comp_per_anchor)], all.x=T, all.y=F, by.x="anchor", by.y="anchor")
final_dt[is.na(num_comp_per_anchor),num_comp_per_anchor:=0]
setnames(final_dt,"num_comp_per_anchor","num_comp_DNA_immune")

final_dt_num_compactors = final_dt[!duplicated(anchor),list(anchor,num_comp_RNA_development,num_comp_RNA_immune,num_comp_DNA_development,num_comp_DNA_immune)] 
final_dt_num_compactors = melt(final_dt_num_compactors,id.vars = "anchor")
setnames(final_dt_num_compactors,c("variable","value"),c("dataset","num_compactor"))
final_dt_num_compactors[,dataset:=gsub("num_comp_","",dataset),by=dataset]
final_dt_num_compactors$anchor=factor(final_dt_num_compactors$anchor,levels=levels(final_dt_orig$anchor))
final_dt_num_compactors$dataset = factor(final_dt_num_compactors$dataset)
ggplot(final_dt_num_compactors ,aes(x = anchor, y = num_compactor, color = dataset,group=1)) +  geom_line() + theme_bw()  
ggplot(final_dt_num_compactors[dataset=="DNA_development"] ,aes(x = anchor, y = num_compactor, color = dataset,group=1)) +  geom_line() + theme_bw() + ylim(0,40)
ggplot(final_dt_num_compactors[dataset=="RNA_development"] ,aes(x = anchor, y = num_compactor, color = dataset,group=1)) + geom_line() + theme_bw() + ylim(0,40)
ggplot(final_dt_num_compactors[dataset=="DNA_immune"] ,aes(x = anchor, y = num_compactor, color = dataset,group=1)) + geom_line() + theme_bw() + ylim(0,40)
ggplot(final_dt_num_compactors[dataset=="RNA_immune"] ,aes(x = anchor, y = num_compactor, color = dataset,group=1))  + geom_line() + theme_bw() + ylim(0,40)

########################################################################################################################
### now below I want to look at each data point and immune vs normal development for each anchor #######################
########################################################################################################################

counts = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/sig_RNA_developmental_anchor_counts_granny_mapping.tsv")


final_dt = data.table() # this data table includes the number of true targets estimated for each dataset
p_seq_err = 0.01
counts = counts[dataset%like%"develop"]
counts[,anchor_counts_per_dataset_day:=sum(count),by=paste(anchor,dataset,Day)]
counts[,min_anchor_count_across_dataset_days:=min(anchor_counts_per_dataset_day),by=anchor] # this gives the minimum anchor count across days and datasets for each anchor
list.dataset = unique(counts$dataset)

for (counter_dataset in 1:length(list.dataset)){
  print(counter_dataset)
  counts_dataset = counts[dataset==list.dataset[counter_dataset]]
  
  list.day = unique(counts_dataset$Day)
  for (counter_day in 1:length(list.day)){
    counts_day = counts_dataset[Day==list.day[counter_day]]
    list.anchor = unique(counts_day$anchor)
    for (counter_anchor in seq(1,length(list.anchor),1)){
      print(counter_anchor)
      counts_anchor = counts_day[anchor==list.anchor[counter_anchor]]
      counts_anchor[,target_frac:=count/anchor_counts_per_dataset_day]
      counts_anchor = counts_anchor[order(-target_frac)]
      
      for (counter_iter in 1:20){
        counts_anchor_downsampled = counts_anchor[sample.int(nrow(counts_anchor), unique(counts_anchor$min_anchor_count_across_dataset_days), replace = TRUE, prob = counts_anchor$target_frac),]
        counts_anchor_downsampled[,anchor_count:=.N,by=list(anchor)]
        counts_anchor_downsampled[,target_count:=.N,by=list(anchor, target)]
        counts_anchor_downsampled = unique(counts_anchor_downsampled[,list(target_count,anchor,target, anchor_count)])
        counts_anchor_downsampled[,target_frac:=target_count/anchor_count]
        counts_anchor_downsampled[,target_frac_99_CI_lower_bound:=binconf(target_count,anchor_count,0.05)[2],by=1:nrow(counts_anchor_downsampled)]
        counts_anchor_downsampled[,target_order:=rank(-target_frac,ties.method = "random"),by=list(anchor)] #ordering targets within each anchor 
        counts_anchor_downsampled = counts_anchor_downsampled[order(-target_frac)]
        
        counts_anchor_downsampled$true_target = 0
        counts_anchor_downsampled$min_HD = 0
        counts_anchor_downsampled$null_error_prob = 0
        counts_anchor_downsampled$true_target[1] = 1
        #    counts_anchor_downsampled = counts_anchor_downsampled[target_count>1]
        
        if (nrow(counts_anchor_downsampled)>1){
          for (counter_target in 2:nrow(counts_anchor_downsampled)){
            counts_anchor_downsampled$min_HD[counter_target] = min(stringdist(counts_anchor_downsampled$target[counter_target],  counts_anchor_downsampled[true_target==1]$target)) # compute minimum hamming distance between the target of interest and the set of true targets
            counts_anchor_downsampled$null_error_prob[counter_target] = p_seq_err^counts_anchor_downsampled$min_HD[counter_target]*sum(counts_anchor_downsampled[true_target==1]$target_frac)
            if (counts_anchor_downsampled$null_error_prob[counter_target] < counts_anchor_downsampled$target_frac_99_CI_lower_bound[counter_target]){
              counts_anchor_downsampled$true_target[counter_target]=1
            }
          }
        }
        
        final_dt = rbind(final_dt, data.table(list.anchor[counter_anchor], list.dataset[counter_dataset], list.day[counter_day], sum(counts_anchor_downsampled$true_target)))
      }
    }
  }
}

names(final_dt) = c("anchor","dataset","day","num_true_targets")

granny_repeat = "CCAGCCATCAGAACCCCAGGAACCATCTAACCAGCCATCAGAACCCCAGGAACCATCTAACCAGCCATCAGAACCCCAGGAACCATCTAA" # 3 consecutive granny repeats
granny_substrings <- unique(sapply(1:(nchar(granny_repeat) - 26), function(i) substr(granny_repeat, i, i + 26)))

final_dt[anchor%in%granny_substrings,anchor_type:="Granny"]
final_dt[,granny_HD:=min(stringdist(anchor, granny_substrings)), by=anchor]
final_dt[granny_HD<7, anchor_type:=paste("Granny: HD",granny_HD,sep="")]

lysine_rich_cdna2 = "CCAGCCAGTGGAACAGCCTGACCACAGTGGTCCTGGATCTGGACCTAAAAATCATATAGGGAAGAAAAGACCAAAAAATACAAAGCCAAAAAATAAAAATAATAAACAACCTAAAGATAAAAAAGGTCCCAAGAAAGATAAAAAACCACACAAAAAGCCCATCAAAGACCCAAGAAAAAAGCCCATCCATCATGATCATCAAGGTGGTAGCCGTAGAAATGGGGATCGTAGAGGAGGTGGTCGACACAGAGGTGATGGTGGTAATCGAGGTGGTCAAA"
lysine_rich_cdna1 = "CCAGCCAGTGGAACAGCCTGACCGCAGTTGTCCTGGAAAGAAAAGACCAAAAAAAATAAAGCCAAAATATAAAAAAAACAAAGATATAAATGGTCACATGAAATGTAAAAAAACCCACAAAATAACCATCAAACACCCAAGAAAAAAGCCCATCCATCATAATCATCAAGGTGGTAGCCGTAGAAATGGGGATCGTAGAGGAGGTGGTCAACACAGAGGTGATGCTGGTAATCGAGGTGGTCAAATAGATGGTGGTCATCGAATTGGTGGTGGTCAAATTGGTGGTGGTCGTAGAATTAGCTATCACAGTGATGGTAGACGTCGTTATGGTTGA"
lysine_substrings_rich1 <- sapply(1:(nchar(lysine_rich_cdna1) - 26), function(i) substr(lysine_rich_cdna1, i, i + 26))
lysine_substrings_rich2 <- sapply(1:(nchar(lysine_rich_cdna2) - 26), function(i) substr(lysine_rich_cdna2, i, i + 26))
lysine_substrings = unique(c(lysine_substrings_rich1,lysine_substrings_rich2))
final_dt[anchor%in%lysine_substrings,anchor_type:="lysine"]
final_dt[anchor%like%"GGTGGTCGTCGAT",anchor_type:="C_terminus_repeat"]

final_dt[,average_num_target:=mean(num_true_targets),by=paste(anchor,dataset,day)]
final_dt = final_dt[!duplicated(paste(anchor,dataset,day))]
final_dt$anchor = factor(final_dt$anchor, levels = unique(final_dt$anchor))
final_dt[,average_num_target_per_anchor_dataset_day:=mean(num_true_targets),by=paste(anchor,dataset,day)]
final_dt = final_dt[!duplicated(paste(anchor,dataset,day))]
final_dt$day= factor(final_dt$day,levels=c(5,8,12))
for (counter_anchor in 1:length(list.anchor)){
  final_dt_anchor = final_dt[anchor==list.anchor[counter_anchor]]
  if (unique(final_dt_anchor$anchor_type)%like%"Granny"){
    p = ggplot(final_dt_anchor, aes(x=day, y=average_num_target_per_anchor_dataset_day,  color=dataset,shape=dataset)) +     geom_point() + theme_bw()  +scale_color_manual(values=brewer.pal(4, 'Set1')) + ggtitle(paste(unique(final_dt_anchor$anchor),unique(final_dt_anchor$anchor_type))) + scale_shape_manual(values= c(0,3,2,8)) + ylim(0,15)
    print(p)
  }
}


compactors_RNA_development_D5 = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/Time_Course_runs/Sponge_RNA_developmental_timecourse_RQ23078_RC/compactors_granny_mapping_anchors_4iteration_D5/sig_ancors_in_normal_or_immune_mapped_to_peter_cdna_compactors.tsv")
compactors_RNA_development_D8 = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/Time_Course_runs/Sponge_RNA_developmental_timecourse_RQ23078_RC/compactors_granny_mapping_anchors_4iteration_D8/sig_ancors_in_normal_or_immune_mapped_to_peter_cdna_compactors.tsv")
compactors_RNA_development_D12 = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/Time_Course_runs/Sponge_RNA_developmental_timecourse_RQ23078_RC/compactors_granny_mapping_anchors_4iteration_D12/sig_ancors_in_normal_or_immune_mapped_to_peter_cdna_compactors.tsv")

compactors_DNA_development_D5 = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/Time_Course_runs/Sponge_DNA_developmental_timecourse_RQ23077/compactors_granny_mapping_anchors_4iteration_D5/sig_ancors_in_normal_or_immune_mapped_to_peter_cdna_compactors.tsv")
compactors_DNA_development_D8 = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/Time_Course_runs/Sponge_DNA_developmental_timecourse_RQ23077/compactors_granny_mapping_anchors_4iteration_D8/sig_ancors_in_normal_or_immune_mapped_to_peter_cdna_compactors.tsv")
compactors_DNA_development_D12 = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/Time_Course_runs/Sponge_DNA_developmental_timecourse_RQ23077/compactors_granny_mapping_anchors_4iteration_D12/sig_ancors_in_normal_or_immune_mapped_to_peter_cdna_compactors.tsv")


compactors_RNA_development_D5 = compactors_RNA_development_D5[num_extended==0]
compactors_RNA_development_D8 = compactors_RNA_development_D8[num_extended==0]
compactors_RNA_development_D12 = compactors_RNA_development_D12[num_extended==0]
compactors_DNA_development_D5 = compactors_DNA_development_D5[num_extended==0]
compactors_DNA_development_D8 = compactors_DNA_development_D8[num_extended==0]
compactors_DNA_development_D12 = compactors_DNA_development_D12[num_extended==0]

compactors_RNA_development_D5[,num_comp_per_anchor:=length(unique(compactor)),by=anchor]
compactors_RNA_development_D8[,num_comp_per_anchor:=length(unique(compactor)),by=anchor]
compactors_RNA_development_D12[,num_comp_per_anchor:=length(unique(compactor)),by=anchor]
compactors_DNA_development_D5[,num_comp_per_anchor:=length(unique(compactor)),by=anchor]
compactors_DNA_development_D8[,num_comp_per_anchor:=length(unique(compactor)),by=anchor]
compactors_DNA_development_D12[,num_comp_per_anchor:=length(unique(compactor)),by=anchor]

compactors_RNA_development_D5[,dataset:="RNA_normal_num_compactor"]
compactors_RNA_development_D5[,day:=5]
compactors_RNA_development_D8[,dataset:="RNA_normal_num_compactor"]
compactors_RNA_development_D8[,day:=8]
compactors_RNA_development_D12[,dataset:="RNA_normal_num_compactor"]
compactors_RNA_development_D12[,day:=12]
compactors_DNA_development_D5[,dataset:="DNA_normal_num_compactor"]
compactors_DNA_development_D5[,day:=5]
compactors_DNA_development_D8[,dataset:="DNA_normal_num_compactor"]
compactors_DNA_development_D8[,day:=8]
compactors_DNA_development_D12[,dataset:="DNA_normal_num_compactor"]
compactors_DNA_development_D12[,day:=12]
compactors = rbind(compactors_RNA_development_D5,compactors_RNA_development_D8,compactors_RNA_development_D12,compactors_DNA_development_D5,compactors_DNA_development_D8,compactors_DNA_development_D12)
compactors = compactors[,list(anchor,dataset,day,num_comp_per_anchor)]
compactors = compactors[!duplicated(paste(anchor,dataset,day))]

###########################################################################################################################
############# now below I want to find the number of targets for each anchor based on the cDNA sequences #################
your_string=fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/Granny_Data/Jacob_assembly/test.fasta",header=F)
your_string = your_string$V1[1]

find_27_chars_after_sequence <- function(your_string, given_sequence) {
  # Find the positions of the given sequence in the string
  positions <- str_locate_all(your_string, given_sequence)[[1]][,2]
  
  # Initialize a list to store the immediate 27 characters after each occurrence
  immediate_chars_list <- c()
  
  # Iterate over each occurrence of the given sequence
  for (counter in 1: length(positions)) {
    # Extract the immediate 27 characters after the occurrence
    immediate_chars <-  substr(your_string, positions[counter] + 1, positions[counter] + 27)
    # Append to the list
    immediate_chars_list <- c(immediate_chars_list, immediate_chars)
  }
  
  return(length(unique(immediate_chars_list)))
}

final_dt$anchor=as.character(final_dt$anchor)
final_dt[,num_comp_cDNA_models:=find_27_chars_after_sequence(your_string,anchor),by=anchor]
final_dt$num_comp_cDNA_models = as.numeric(final_dt$num_comp_cDNA_models)
final_dt_num_compactors = final_dt[!duplicated(anchor),list(anchor,num_comp_cDNA_models)]
final_dt_num_compactors = melt(final_dt_num_compactors,id.vars = "anchor")
setnames(final_dt_num_compactors,c("value"),c("num_compactor_cDNA_model"))
final_dt_num_compactors$anchor=factor(final_dt_num_compactors$anchor,levels=levels(final_dt_orig$anchor))
ggplot(final_dt_num_compactors ,aes(x = anchor, y = num_compactor_cDNA_model, group=1)) +  geom_line() + theme_bw() + ylim(0,40)


find_27_chars_after_sequence_isoseq <- function(your_string, given_sequence) {
  # Find the positions of the given sequence in the string
  positions <- str_locate_all(your_string, given_sequence)[[1]][,2]
  
  # Initialize a list to store the immediate 27 characters after each occurrence
  immediate_chars_list <- c()
  
  # Iterate over each occurrence of the given sequence
  for (counter in 1: length(positions)) {
    # Extract the immediate 27 characters after the occurrence
    immediate_chars <-  substr(your_string, positions[counter] + 1, positions[counter] + 27)
    # Append to the list
    immediate_chars_list <- c(immediate_chars_list, immediate_chars)
  }
  
  return(unique(immediate_chars_list))
}

your_string = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/Granny_Data/isoseq_data/granny_deduped_isoseq.fasta",header=F) 
your_string = your_string[!V1%like%">"]
your_string = paste(your_string$V1,collapse = "") 
isoseq_granny_targets = find_27_chars_after_sequence_isoseq(your_string,"GCCATCAGAACCCCAGGAACCATCTAA")
cDNA_granny_targets = find_27_chars_after_sequence_isoseq(your_string,"GCCATCAGAACCCCAGGAACCATCTAA")
###########################################################################################################################
#################### BELOW I WANT TO SHOW THE DISTRIBUTION OF TARGET COUNTS PER ANCHOR FOR DNA-SEQ #########################
counts = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/sig_RNA_developmental_anchor_counts_granny_mapping.tsv")

final_dt = data.table() # this data table includes the number of true targets estimated for each dataset
p_seq_err = 0.01
counts = counts[dataset%like%"DNA" | dataset%like%"RNA_immune" ]
counts[,anchor_counts_per_dataset:=sum(count),by=paste(anchor,dataset)]
counts[,min_anchor_count_per_dataset:=min(anchor_counts_per_dataset)] # this gives the minimum anchor count per dataset across all datasets I use this to downsample the anchor count in datasets that have higher counts
counts[,target_count_per_dataset:=sum(count),by=paste(anchor,target,dataset)]
list.dataset = unique(counts$dataset)


## now below I want to classify anchors as granny, HD with granny, N_terminus, C_terminus, ....
granny_repeat = "CCAGCCATCAGAACCCCAGGAACCATCTAACCAGCCATCAGAACCCCAGGAACCATCTAACCAGCCATCAGAACCCCAGGAACCATCTAA" # 3 consecutive granny repeats
granny_substrings <- unique(sapply(1:(nchar(granny_repeat) - 26), function(i) substr(granny_repeat, i, i + 26)))

counts[anchor%in%granny_substrings,anchor_type:="Granny"]
counts[,granny_HD:=min(stringdist(anchor, granny_substrings)), by=anchor]
counts[granny_HD<7, anchor_type:=paste("Granny: HD",granny_HD,sep="")]

lysine_rich_cdna2 = "CCAGCCAGTGGAACAGCCTGACCACAGTGGTCCTGGATCTGGACCTAAAAATCATATAGGGAAGAAAAGACCAAAAAATACAAAGCCAAAAAATAAAAATAATAAACAACCTAAAGATAAAAAAGGTCCCAAGAAAGATAAAAAACCACACAAAAAGCCCATCAAAGACCCAAGAAAAAAGCCCATCCATCATGATCATCAAGGTGGTAGCCGTAGAAATGGGGATCGTAGAGGAGGTGGTCGACACAGAGGTGATGGTGGTAATCGAGGTGGTCAAA"
lysine_rich_cdna1 = "CCAGCCAGTGGAACAGCCTGACCGCAGTTGTCCTGGAAAGAAAAGACCAAAAAAAATAAAGCCAAAATATAAAAAAAACAAAGATATAAATGGTCACATGAAATGTAAAAAAACCCACAAAATAACCATCAAACACCCAAGAAAAAAGCCCATCCATCATAATCATCAAGGTGGTAGCCGTAGAAATGGGGATCGTAGAGGAGGTGGTCAACACAGAGGTGATGCTGGTAATCGAGGTGGTCAAATAGATGGTGGTCATCGAATTGGTGGTGGTCAAATTGGTGGTGGTCGTAGAATTAGCTATCACAGTGATGGTAGACGTCGTTATGGTTGA"
lysine_substrings_rich1 <- sapply(1:(nchar(lysine_rich_cdna1) - 26), function(i) substr(lysine_rich_cdna1, i, i + 26))
lysine_substrings_rich2 <- sapply(1:(nchar(lysine_rich_cdna2) - 26), function(i) substr(lysine_rich_cdna2, i, i + 26))
lysine_substrings = unique(c(lysine_substrings_rich1,lysine_substrings_rich2))
counts[anchor%in%lysine_substrings,anchor_type:="lysine"]
counts[anchor%like%"GGTGGTCGTCGAT",anchor_type:="C_terminus_repeat"]


for (counter_dataset in 1:length(list.dataset)){
  print(counter_dataset)
  counts_dataset = counts[dataset==list.dataset[counter_dataset]]
  counts_dataset[,anchor_count:=sum(count),by=list(anchor)]
  counts_dataset[,target_count:=sum(count),by=list(anchor, target)]
  counts_dataset = unique(counts_dataset[,list(target_count,anchor,target, anchor_count, min_anchor_count_per_dataset,anchor_type)])
  
  counts_dataset[,target_frac:=target_count/anchor_count]
  counts_dataset[,target_order:=rank(-target_frac,ties.method = "random"),by=list(anchor)] #ordering targets within each anchor 
  counts_dataset = counts_dataset[target_order<7]
  counts_dataset$target_order = factor(counts_dataset$target_order)
  counts_dataset$anchor=factor(counts_dataset$anchor,levels=levels(final_dt_orig$anchor))
  ggplot(counts_dataset, aes(x=anchor, y=target_count,  color=anchor_type,shape = target_order)) +     geom_point() + theme_bw()  +scale_color_manual(values=brewer.pal(length(unique(counts_dataset$anchor_type)), 'Set1')) +scale_shape_manual(values= c(1,3,8,17,6,9)) 
  ggplot(counts_dataset, aes(x=anchor, y=target_frac,  color=anchor_type,shape = target_order)) +     geom_point() + theme_bw()  +scale_color_manual(values=brewer.pal(length(unique(counts_dataset$anchor_type)), 'Set1')) +scale_shape_manual(values= c(1,3,8,17,6,9)) + ylim(0,0.7)
  
}
