library(data.table)

## DHR123: Archaeocytes (splash sample id 0)
## FLUO: Choanocytes (splash sample id 1)
## NEITHER: cell types other than archaeocytes and choanocytes (splash sample id 2)
granny_FACS_counts = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/Bulk_RNA_runs/Sponge_Jacob_FACS_RNA/satc_granny_anchor/granny_anchor_counts.txt")
granny_FACS_counts[V1==0,cell_type:="Archaeocytes"]
granny_FACS_counts[V1==1,cell_type:="Choanocytes"]
granny_FACS_counts[V1==2,cell_type:="Neither"]
granny_FACS_counts[,V1:=NULL]

granny_FACS_counts_reshape = reshape(granny_FACS_counts[,list(V3,V4,cell_type)], idvar = "V3", timevar = "cell_type", direction = "wide")
granny_FACS_counts_reshape[is.na(granny_FACS_counts_reshape)] = 0

names(granny_FACS_counts)[1:3] = c("anchor","target","count")

final_dt = data.table() # this data table includes the number of true targets estimated for each cell_type
p_seq_err = 0.01
granny_FACS_counts[,anchor_counts_per_cell_type:=sum(count),by=paste(anchor,cell_type)]
granny_FACS_counts[,min_anchor_count_per_cell_type:=min(anchor_counts_per_cell_type)] # this gives the minimum anchor count per cell_type across all cell_types I use this to downsample the anchor count in cell_types that have higher counts
granny_FACS_counts[,target_count_per_cell_type:=sum(count),by=paste(anchor,target,cell_type)]
list.cell_type = unique(granny_FACS_counts$cell_type)
for (counter_cell_type in 1:length(list.cell_type)){
  print(counter_cell_type)
  counts_cell_type = granny_FACS_counts[cell_type==list.cell_type[counter_cell_type]]
  counts_cell_type[,anchor_count:=sum(count),by=list(anchor)]
  counts_cell_type[,target_count:=sum(count),by=list(anchor, target)]
  counts_cell_type = unique(counts_cell_type[,list(target_count,anchor,target, anchor_count, min_anchor_count_per_cell_type)])
  
  counts_cell_type[,target_frac:=target_count/anchor_count]
  counts_cell_type[,target_order:=rank(-target_frac,ties.method = "random"),by=list(anchor)] #ordering targets within each anchor 
  
  list.anchor = unique(counts_cell_type$anchor)
  for (counter_anchor in seq(1,length(list.anchor),1)){
    counts_anchor = counts_cell_type[anchor==list.anchor[counter_anchor]]
    counts_anchor = counts_anchor[order(-target_frac)]
    
    for (counter_iter in 1:20){
      counts_anchor_downsampled = counts_anchor[sample.int(nrow(counts_anchor), unique(counts_anchor$min_anchor_count_per_cell_type), replace = TRUE, prob = counts_anchor$target_frac),]
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
      
      final_dt = rbind(final_dt, data.table(list.anchor[counter_anchor],list.cell_type[counter_cell_type],sum(counts_anchor_downsampled$true_target)))
    }
    
  }
}

names(final_dt) = c("anchor","dataset","num_true_targets")

final_dt[,average_num_target_per_anchor_dataset:=mean(num_true_targets),by=paste(anchor,dataset)]
final_dt = final_dt[!duplicated(paste(anchor,dataset))]
ggplot(final_dt, aes(x=anchor, y=average_num_target_per_anchor_dataset,  color=anchor_type,shape = dataset)) +     geom_point() + theme_bw()  +scale_color_manual(values=brewer.pal(length(unique(final_dt$anchor_type)), 'Set1')) +scale_shape_manual(values= c(1,3,8,17,6,9)) + ylim(0,15)# 6,9



######## now we compute number of true targets without downsampling #####################

final_dt = data.table() # this data table includes the number of true targets estimated for each cell_type
p_seq_err = 0.01
granny_FACS_counts[,anchor_counts_per_cell_type:=sum(count),by=paste(anchor,cell_type)]
granny_FACS_counts[,min_anchor_count_per_cell_type:=min(anchor_counts_per_cell_type)] # this gives the minimum anchor count per cell_type across all cell_types I use this to downsample the anchor count in cell_types that have higher counts
granny_FACS_counts[,target_count_per_cell_type:=sum(count),by=paste(anchor,target,cell_type)]
list.cell_type = unique(granny_FACS_counts$cell_type)
for (counter_cell_type in 1:length(list.cell_type)){
  print(counter_cell_type)
  counts_cell_type = granny_FACS_counts[cell_type==list.cell_type[counter_cell_type]]
  counts_cell_type[,anchor_count:=sum(count),by=list(anchor)]
  counts_cell_type[,target_count:=sum(count),by=list(anchor, target)]
  counts_cell_type = unique(counts_cell_type[,list(target_count,anchor,target, anchor_count, min_anchor_count_per_cell_type)])
  
  counts_cell_type[,target_frac:=target_count/anchor_count]
  counts_cell_type[,target_order:=rank(-target_frac,ties.method = "random"),by=list(anchor)] #ordering targets within each anchor 
  
  counts_cell_type = counts_cell_type[order(-target_frac)]
  counts_anchor_downsampled = counts_cell_type

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
  
  final_dt = rbind(final_dt, data.table(list.anchor[counter_anchor],list.cell_type[counter_cell_type],sum(counts_anchor_downsampled$true_target)))
}

names(final_dt) = c("anchor","dataset","num_true_targets")



## now I am the plot for showing the top targets for each celltype in FACS
granny_FACS_counts_reshape = setorder(granny_FACS_counts_reshape,-V4.Neither)
granny_FACS_counts$target=factor(granny_FACS_counts$target,levels= granny_FACS_counts_reshape$V3)

granny_FACS_counts$cell_type = factor(granny_FACS_counts$cell_type)
ggplot(granny_FACS_counts[target%in%granny_FACS_counts_reshape$V3[1:30]], aes(fill=cell_type, y=count, x=target)) +     geom_bar(position="dodge", stat="identity") + theme_bw() +  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +annotation_logticks()  +theme(axis.text.x=element_blank())
