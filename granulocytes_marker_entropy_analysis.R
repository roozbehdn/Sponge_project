library(data.table)
library(Biostrings)

high_entropy_anchors_to_granulocyte_markers = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/10X_runs/Sponge_SRP216435_all_samples_together/high_entropy_anchors_to_granulocyte_marker.sam",header=F,fill=T)

high_entropy_anchors = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/10X_runs/Sponge_SRP216435_all_samples_together/result.after_correction.scores.top_target_entropy.tsv")

high_entropy_anchors_to_granulocyte_markers[V2==16,V10:=as.character(reverseComplement(DNAString(V10))),by=V10]

high_entropy_anchors_to_granulocyte_markers = merge(high_entropy_anchors_to_granulocyte_markers,high_entropy_anchors[,list(anchor,target_entropy)],all.x=T,all.y=F,by.x="V10",by.y="anchor")
