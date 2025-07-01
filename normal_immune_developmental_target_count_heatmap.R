library(data.table)
library(glmnet)
library(ggplot2)
library(gridExtra)
library(pheatmap)
library(viridis)

counts = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/anchor_counts_sig_immune_granny_mapped_merged_with_normal_dev.tsv")
selected_anchors = unique(counts$anchor)
SPLASH_output = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/Time_Course_runs/Sponge_RNA_immune_response_RQ23179/result.after_correction.scores.tsv")
SPLASH_output = SPLASH_output[anchor%in%selected_anchors]

counts[,Day:=gsub("D","",Day)]
counts$Day = as.numeric(counts$Day)
counts$Day=factor(counts$Day,levels=sort(unique(counts$Day)))

for (counter in 1:nrow(SPLASH_output)){
  tryCatch({
    anchor_interest = SPLASH_output$anchor[counter]
    effectsize_interest = SPLASH_output$effect_size_bin[counter]
    entropy_interest = SPLASH_output$target_entropy[counter]
    
    counts_anchor = counts[anchor==anchor_interest] # target counts for the selected anchor
    counts_anchor[,target_count := sum(count),by=target] # compute total counts for each target
    counts_anchor[,target_frac:=target_count/sum(counts_anchor$count)] # fraction of anchor read for each target
    counts_anchor$target = factor(counts_anchor$target,levels =  counts_anchor[order(-target_frac),list(target,target_frac)][!duplicated(target)]$target) # sort targets based on their counts
    top_targets = unique( sort(counts_anchor$target))
    if (length(top_targets)>4){
      top_targets = top_targets[1:4]
      counts_anchor = counts_anchor[target %in%top_targets]
    }
    counts_anchor[,anchor_count_per_sample:=sum(count),by=sample_name]
    counts_anchor[,target_frac_per_sample:=count/anchor_count_per_sample]
    
    g = ggplot(data = counts_anchor,aes(x = Day, y = target_frac_per_sample, color = target, group=paste(target,Treatment), shape= Treatment)) +
      geom_line(size = 1, alpha = 0.5) +
      geom_point(size = 2) + theme_bw() + ggtitle(anchor_interest)
    
    pdf(file=paste("/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/Time_Course_runs/Sponge_RNA_immune_response_RQ23179/satc_sig_anchor/heatmap_plots/heatmap_plots_for_aligned_anchors_to_peter_cdna/anchor_",anchor_interest,"_entropy_",entropy_interest,"_effectsize_",effectsize_interest,".pdf", sep = ""), width = 9, height = 6)
    print(g)
    dev.off()
    
  },error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}
