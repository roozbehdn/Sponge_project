library(data.table)
library(stringdist)
library(Seurat)

file_path_1002 = "/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/SRP198321/GSM3764781_SRR9051002_latTII1a_raw_gene_bc_matrices_h5.h5"
file_path_1003 = "/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/SRP198321/GSM3764782_SRR9051003_latTII1b_raw_gene_bc_matrices_h5.h5"
file_path_1004 = "/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/SRP198321/GSM3764783_SRR9051004_latTII2_raw_gene_bc_matrices_h5.h5"
file_path_1005 <- "/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/SRP198321/GSM3764784_SRR9051005_larva1_raw_gene_bc_matrices_h5.h5"
file_path_1006 <- "/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/SRP198321/GSM3764785_SRR9051006_larva2_raw_gene_bc_matrices_h5.h5"
file_path_1007 <- "/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/SRP198321/GSM3764786_SRR9051007_larva3_raw_gene_bc_matrices_h5.h5"

data_1002 = Read10X_h5(file_path_1002)
seurat_object_1002 <- CreateSeuratObject(counts = data_1002)
data_1003 = Read10X_h5(file_path_1003)
seurat_object_1003 <- CreateSeuratObject(counts = data_1003)
data_1004 = Read10X_h5(file_path_1004)
seurat_object_1004 <- CreateSeuratObject(counts = data_1004)
data_1005 = Read10X_h5(file_path_1005)
seurat_object_1005 <- CreateSeuratObject(counts = data_1005)
data_1006 = Read10X_h5(file_path_1006)
seurat_object_1006 <- CreateSeuratObject(counts = data_1006)
data_1007 = Read10X_h5(file_path_1007)
seurat_object_1007 <- CreateSeuratObject(counts = data_1007)

counts_1002 = data.frame(seurat_object_1002@meta.data)
counts_1002$cell_barcode=rownames(counts_1002)
counts_1002 = data.table(counts_1002)
counts_1002$sample="SRR9051002"
counts_1003 = data.frame(seurat_object_1003@meta.data)
counts_1003$cell_barcode=rownames(counts_1003)
counts_1003 = data.table(counts_1003)
counts_1003$sample="SRR9051003"
counts_1004 = data.frame(seurat_object_1004@meta.data)
counts_1004$cell_barcode=rownames(counts_1004)
counts_1004 = data.table(counts_1004)
counts_1004$sample="SRR9051004"
counts_1005 = data.frame(seurat_object_1005@meta.data)
counts_1005$cell_barcode=rownames(counts_1005)
counts_1005 = data.table(counts_1005)
counts_1005$sample="SRR9051005"
counts_1006 = data.frame(seurat_object_1006@meta.data)
counts_1006$cell_barcode=rownames(counts_1006)
counts_1006 = data.table(counts_1006)
counts_1006$sample="SRR9051006"
counts_1007 = data.frame(seurat_object_1007@meta.data)
counts_1007$cell_barcode=rownames(counts_1007)
counts_1007 = data.table(counts_1007)
counts_1007$sample="SRR9051007"

counts_all = rbind(counts_1002,counts_1003,counts_1004,counts_1005,counts_1006,counts_1007)
counts_all = counts_all[,list(nCount_RNA,cell_barcode,sample)]
counts_all[,cell_barcode:=strsplit(cell_barcode,split="-")[[1]][1],by=cell_barcode]
counts_all = counts_all[nCount_RNA>5]
counts_all[,cell:=paste(sample,cell_barcode,sep="_"),by=1:nrow(counts_all)]

ciona_1002 = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_SRP198321/SRR9051002/Ayelet_seq_satc/satcOut112.tsv",header = F)
ciona_1003 = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_SRP198321/SRR9051003/Ayelet_seq_satc/satcOut112.tsv",header = F)
ciona_1004 = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_SRP198321/SRR9051004/Ayelet_seq_satc/satcOut112.tsv",header = F)
ciona_1005 = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_SRP198321/SRR9051005/Ayelet_seq_satc/satcOut112.tsv",header = F)
ciona_1006 = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_SRP198321/SRR9051006/Ayelet_seq_satc/satcOut112.tsv",header = F)
ciona_1007 = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_SRP198321/SRR9051007/Ayelet_seq_satc/satcOut112.tsv",header = F)

metadata = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/utility_files/metadata/non_model_organisms/Ciona_intestinalis_SRP198321/ciona_metadata_modified.txt")
metadata[,sample_new:=sample]
metadata[sample=="SRR9051005",sample_new:="SRR9051007"]
metadata[sample=="SRR9051006",sample_new:="SRR9051005"]
metadata[sample=="SRR9051007",sample_new:="SRR9051006"]
metadata[,sample:=NULL]
setnames(metadata,"sample_new","sample")
setnames(metadata,"cell","cell_barcode")
metadata[,cell:=paste(sample,cell_barcode,sep="_"),by=1:nrow(metadata)]
metadata = metadata[sample%in%c("SRR9051002","SRR9051003","SRR9051004","SRR9051005","SRR9051006","SRR9051007")]
metadata[sample%in%c("SRR9051002","SRR9051003","SRR9051004"),dev_stage:="late_tail_bud_II"]
metadata[sample%in%c("SRR9051005","SRR9051006","SRR9051007"),dev_stage:="larva"]
metadata[,cell_type_dev:=paste(cell_type,dev_stage,sep="_"),by=paste(cell_type,dev_stage)]

#######################################################################################################################################
#### below is the metadata file with detailed marker information for mesenchymal celltypes that the authors sent me later #############
mesenchyme = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/utility_files/metadata/non_model_organisms/Ciona_intestinalis_SRP198321/ciona_mesenchyme_annotation_modified.txt")
mesenchyme[,sample_new:=sample]
mesenchyme[sample=="SRR9051005",sample_new:="SRR9051007"]
mesenchyme[sample=="SRR9051006",sample_new:="SRR9051005"]
mesenchyme[sample=="SRR9051007",sample_new:="SRR9051006"]
mesenchyme[,sample:=NULL]
setnames(mesenchyme,"sample_new","sample")
setnames(mesenchyme,"cell","cell_barcode")
mesenchyme[,cell:=paste(sample,cell_barcode,sep="_"),by=1:nrow(mesenchyme)]
#######################################################################################################################################
#######################################################################################################################################

ciona_1002[,target_count_per_sample := sum(V5), by =V4]
ciona_1003[,target_count_per_sample := sum(V5), by =V4]
ciona_1004[,target_count_per_sample := sum(V5), by =V4]
ciona_1005[,target_count_per_sample := sum(V5), by =V4]
ciona_1006[,target_count_per_sample := sum(V5), by =V4]
ciona_1007[,target_count_per_sample := sum(V5), by =V4]

ciona_1002[,sample:="SRR9051002"]
ciona_1003[,sample:="SRR9051003"]
ciona_1004[,sample:="SRR9051004"]
ciona_1005[,sample:="SRR9051005"]
ciona_1006[,sample:="SRR9051006"]
ciona_1007[,sample:="SRR9051007"]
ciona_1002[,dev_stage:="late_tail_bud_II"]
ciona_1003[,dev_stage:="late_tail_bud_II"]
ciona_1004[,dev_stage:="late_tail_bud_II"]
ciona_1005[,dev_stage:="larva"]
ciona_1006[,dev_stage:="larva"]
ciona_1007[,dev_stage:="larva"]

ciona_all = rbind(ciona_1002,ciona_1003,ciona_1004,ciona_1005,ciona_1006,ciona_1007)
ciona_all[,V1:=NULL]
ciona_all[,total_target_count:=sum(V5),by=V4]
setnames(ciona_all,c("V2","V3","V4","V5"),c("cell_barcode","anchor","target","count_per_cell"))

## merging with metadata
ciona_all[,cell:=paste(sample,cell_barcode,sep="_"),by=paste(sample,cell_barcode)]
ciona_all = merge(ciona_all,metadata[,list(cell,cell_type,cell_type_dev)],all.x=T,all.y=F,by.x="cell",by.y="cell")

#### merging with the metadata file for mesenchymal celltypes ################
setnames(mesenchyme, "cell_type", "mesenchymal_cell_type")
ciona_all = merge(ciona_all,mesenchyme[,list(cell,mesenchymal_cell_type)],all.x=T,all.y=F,by.x="cell",by.y="cell")
##############################################################################

################################################################################################
########## now I look at only those cells that are in mesenchymal celltypes ####################
ciona_all_mesenchymal = ciona_all[cell_type%like%"mesen"]

## the mesenchymal clusters that I look at are 32_mesenchyme (larva), 18_mesenchyme (late tail bud), 18_mesenchyme (larva), 32_mesenchyme (late tail bud), 81_mesenchyme (larva), 82_mesenchyme_larva ###
metadata_mesenchyme = metadata[cell_type%like%"mesen"]
metadata_mesenchyme = merge(metadata_mesenchyme, mesenchyme[,list(cell,mesenchymal_cell_type)],all.x=T,all.y=F,by.x="cell",by.y="cell")
metadata_mesenchyme[,num_cell_per_mesenchym_subcluster:=length(unique(cell)),by=paste(cell_type_dev, mesenchymal_cell_type)]
metadata_mesenchyme[,num_cell_per_mesenchyme_cluster:=length(unique(cell)),by=cell_type_dev]
metadata_mesenchyme[,mesenchym_subcluster_fraction:=num_cell_per_mesenchym_subcluster/num_cell_per_mesenchyme_cluster]
metadata_mesenchyme_summary = unique(metadata_mesenchyme[,list(cell_type_dev,mesenchymal_cell_type,num_cell_per_mesenchym_subcluster,num_cell_per_mesenchyme_cluster,mesenchym_subcluster_fraction)])

ciona_all_mesenchymal[,num_anchor_pos_cell_per_mesenchym_subcluster:=length(unique(cell)),by=paste(cell_type_dev, mesenchymal_cell_type)]
ciona_all_mesenchymal[,num_anchor_pos_cell_per_mesenchyme_cluster:=length(unique(cell)),by=cell_type_dev]
ciona_all_mesenchymal[,anchor_pos_mesenchym_subcluster_fraction:=num_anchor_pos_cell_per_mesenchym_subcluster/num_anchor_pos_cell_per_mesenchyme_cluster]
ciona_all_mesenchymal_summary = unique(ciona_all_mesenchymal[,list(cell_type_dev,mesenchymal_cell_type,num_anchor_pos_cell_per_mesenchym_subcluster,num_anchor_pos_cell_per_mesenchyme_cluster,anchor_pos_mesenchym_subcluster_fraction)])

## now I merge the fraction of mesenchymal subclusters from metadata with the fraction of mesychmal subclusters in anchor positive cells to see if there is an enrichment
ciona_all_mesenchymal_summary = merge(ciona_all_mesenchymal_summary,metadata_mesenchyme_summary,all.x=T,all.y = F,by.x=c("cell_type_dev", "mesenchymal_cell_type"),by.y=c("cell_type_dev", "mesenchymal_cell_type"))
ciona_all_mesenchymal_summary[,fraction_expressing_cell_per_cluster:=num_anchor_pos_cell_per_mesenchyme_cluster/num_cell_per_mesenchyme_cluster]

ggplot(ciona_all_mesenchymal_summary[cell_type_dev%like%"larva" & mesenchymal_cell_type%like%"Mech"], aes(x=fraction_expressing_cell_per_cluster, y=mesenchym_subcluster_fraction,size=num_anchor_pos_cell_per_mesenchym_subcluster)) + geom_point() + theme_bw() +geom_text_repel(aes(label = cell_type_dev))

### now I want to look at the difference between Hlx+ and Hlx- clusters
ciona_all_mesenchymal_summary_hlx_pos = ciona_all_mesenchymal_summary[cell_type_dev%like%"larva" & mesenchymal_cell_type=="Hlx+ mesenchyme"]
ciona_all_mesenchymal_summary_hlx_neg = ciona_all_mesenchymal_summary[cell_type_dev%like%"larva" & mesenchymal_cell_type%like%"Hlx-"]


## merging with counts to get the normalized counts per cell
ciona_all = merge(ciona_all, counts_all[,list(nCount_RNA,cell)],all.x=T,all.y=F,by.x="cell",by.y="cell")
ciona_all[,anchor_count_per_cell:=sum(count_per_cell),by=cell]
ciona_all[,normalized_anchor_count_per_cell:=anchor_count_per_cell/nCount_RNA * 10^5]
ciona_all_uniq_cell = ciona_all[!duplicated(cell)]
ciona_all_uniq_cell[,average_normalized_count_per_celltype:=mean(normalized_anchor_count_per_cell),by=cell_type_dev]

## below I want to look at the celltype enrichment for the anchor
a=ciona_all[,length(unique(cell)),by=cell_type_dev]
b=metadata[,length(unique(cell)),by=cell_type_dev]
a=merge(a,b,all.x=T,all.y=F,by.x="cell_type_dev",by.y="cell_type_dev")
names(a)[2:3]=c("num_expressing_cells","num_tot_cells")
a[,frac:=num_expressing_cells/num_tot_cells]
a = merge(a,ciona_all_uniq_cell[!duplicated(cell_type_dev),list(cell_type_dev,average_normalized_count_per_celltype)])
a[,dev_stage:="larva"]
a[grepl("late",cell_type_dev),dev_stage:="late_tail_bud_II"]
a1=a[num_expressing_cells>9 & frac>0.05]

ggplot(a1, aes(x = num_expressing_cells, size = frac, color = dev_stage, y = average_normalized_count_per_celltype)) + geom_point() + theme_bw()+  geom_text_repel(aes(label = cell_type_dev),    # Add text labels
                                                                                                                                                                   size = 3,              # Size of text
                                                                                                                                                                   max.overlaps = 10)  + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +annotation_logticks()


ciona_all[target=="CAACAAAGCCGATGGTTACTATGACAC",model:="1391x2_1392x1"]
ciona_all[target=="CACCAAAGCCGATGGTTACTACGATGG",model:="1391x5_1392x3"]
ciona_all[target=="CACCAAAGCCGATGGTTACTATGACAA",model:="1391"]
ciona_all[target=="CACCAAAGCCGATGGTTACTATGACAC",model:="1391x2_1392x3"]
ciona_all[target=="CACCAAGGCCGATGGTTACTATGACAC",model:="1392"]
ciona_all[target=="GAACAAAGCCGATGGTTACTATGACAC",model:="1391_1392"]

## below i want to compute the minimum hamming distance relative to the targets from the model
model_targets = c("CAACAAAGCCGATGGTTACTATGACAC", "CACCAAAGCCGATGGTTACTACGATGG", "CACCAAAGCCGATGGTTACTATGACAA", "CACCAAAGCCGATGGTTACTATGACAC", "CACCAAGGCCGATGGTTACTATGACAC", "GAACAAAGCCGATGGTTACTATGACAC")
model = c("1391x2_1392x1", "1391x5_1392x3", "1391", "1391x2_1392x3", "1392", "1391_1392")
ciona_all_uniq_target=ciona_all[!duplicated(target)]
ciona_all_uniq_target[,min_ham_dist:=NA]
setorder(ciona_all_uniq_target,-total_target_count)
ciona_all_uniq_target = ciona_all_uniq_target[total_target_count>1]
for (counter in 1:nrow(ciona_all_uniq_target)){
  him_dists_to_ref_targets =stringdist(ciona_all_uniq_target$target[counter],  model_targets)
  ciona_all_uniq_target$min_ham_dist[counter]=min(him_dists_to_ref_targets)
  ciona_all_uniq_target$closest_ref_tar_model[counter]=model[which.min(him_dists_to_ref_targets)]
  ciona_all_uniq_target$closest_ref_tar_sequence[counter]=model_targets[which.min(him_dists_to_ref_targets)]
  print(counter)
}
ciona_all_uniq_target = ciona_all_uniq_target[,list(target,total_target_count,closest_ref_tar_model,closest_ref_tar_sequence, min_ham_dist,model)]

### for the supplementary figure for the ciona cry11bb targets I keep at most 100 targets for each reference target based on their counts #### 
split_ciona_all_uniq_target = split(ciona_all_uniq_target,ciona_all_uniq_target$closest_ref_tar_sequence)

# Keep the top 100 target for each reference target
top_100 <- lapply(split_ciona_all_uniq_target, function(group_df) {
  group_df[order(-group_df$total_target_count), ][1:min(100, nrow(group_df)), ]
})

ciona_all_uniq_target_top100 = do.call(rbind, top_100) 
ciona_all_uniq_target_top100$target = factor(ciona_all_uniq_target_top100$target,levels=ciona_all_uniq_target_top100$target)

library(scales)
ggplot(ciona_all_uniq_target_top100[closest_ref_tar_model=="1391x2_1392x1"], aes(x=target, y=total_target_count,fill="darksalmon")) + geom_bar(stat = "identity") + theme_bw() + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +annotation_logticks()
ggplot(ciona_all_uniq_target_top100[closest_ref_tar_model=="1391x5_1392x3"], aes(x=target, y=total_target_count,fill="darksalmon")) + geom_bar(stat = "identity") + theme_bw() + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +annotation_logticks()
ggplot(ciona_all_uniq_target_top100[closest_ref_tar_model=="1391"], aes(x=target, y=total_target_count,fill="darksalmon")) + geom_bar(stat = "identity") + theme_bw() + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +annotation_logticks()
ggplot(ciona_all_uniq_target_top100[closest_ref_tar_model=="1391x2_1392x3"], aes(x=target, y=total_target_count,fill="darksalmon")) + geom_bar(stat = "identity") + theme_bw() + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +annotation_logticks()
ggplot(ciona_all_uniq_target_top100[closest_ref_tar_model=="1392"], aes(x=target, y=total_target_count,fill="darksalmon")) + geom_bar(stat = "identity") + theme_bw() + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +annotation_logticks()
ggplot(ciona_all_uniq_target_top100[closest_ref_tar_model=="1391_1392"], aes(x=target, y=total_target_count,fill="darksalmon")) + geom_bar(stat = "identity") + theme_bw() + scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +annotation_logticks()

## now I write fasta file for MSA alignment for the target variants of each reference target
write.table(ciona_all_uniq_target_top100[closest_ref_tar_model=="1391x2_1392x1"]$target,"/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/ciona_cry11bb_top_100_target_per_ref_target_1391x2_1392x1.txt",sep="\t",row.names=F,quote=F,col.names = F)
write.table(ciona_all_uniq_target_top100[closest_ref_tar_model=="1391x5_1392x3"]$target,"/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/ciona_cry11bb_top_100_target_per_ref_target_1391x5_1392x3.txt",sep="\t",row.names=F,quote=F,col.names = F)
write.table(ciona_all_uniq_target_top100[closest_ref_tar_model=="1391"]$target,"/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/ciona_cry11bb_top_100_target_per_ref_target_1391.txt",sep="\t",row.names=F,quote=F,col.names = F)
write.table(ciona_all_uniq_target_top100[closest_ref_tar_model=="1391x2_1392x3"]$target,"/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/ciona_cry11bb_top_100_target_per_ref_target_1391x2_1392x3.txt",sep="\t",row.names=F,quote=F,col.names = F)
write.table(ciona_all_uniq_target_top100[closest_ref_tar_model=="1392"]$target,"/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/ciona_cry11bb_top_100_target_per_ref_target_1392.txt",sep="\t",row.names=F,quote=F,col.names = F)
write.table(ciona_all_uniq_target_top100[closest_ref_tar_model=="1391_1392"]$target,"/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/ciona_cry11bb_top_100_target_per_ref_target_1391_1392.txt",sep="\t",row.names=F,quote=F,col.names = F)



#################################################################################################################
########## BELOW I MAKE THE SEQUENCES LOGO FOR THE TARGETS ASSOCIATED TO EACH REFERENCE TARGET ##################
library(ggseqlogo)

create_weighted_pfm <- function(sequences, weights) { ## this functions generates the weights for each target based on its total count
  # Get the length of the sequences
  seq_length <- nchar(sequences[1])
  
  # Define the alphabet (A, T, C, G for DNA)
  alphabet <- c("A", "T", "C", "G")
  
  # Initialize an empty PFM
  pfm <- matrix(0, nrow = length(alphabet), ncol = seq_length, 
                dimnames = list(alphabet, 1:seq_length))
  
  # Fill in the PFM with weighted counts
  for (i in seq_along(sequences)) {
    for (pos in seq_len(seq_length)) {
      base <- substr(sequences[i], pos, pos)
      pfm[base, pos] <- pfm[base, pos] + weights[i]
    }
  }
  
  # Normalize columns to get frequencies
  pfm <- sweep(pfm, 2, colSums(pfm), "/")
  return(pfm)
}

targets = as.character(ciona_all_uniq_target_top100[closest_ref_tar_model=="1391x5_1392x3"]$target)
target_counts = ciona_all_uniq_target_top100[closest_ref_tar_model=="1391x5_1392x3"]$total_target_count

# Create weighted PFM
weighted_pfm <- create_weighted_pfm(targets, target_counts)

# Generate sequence logo
ggseqlogo(weighted_pfm, method = "custom") +
  ggtitle("Weighted Sequence Logo")

##############################################################################################
##############################################################################################
