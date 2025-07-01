library(data.table)
library(pheatmap)
library(viridis)
library(RColorBrewer)
library(ggplot2)

f = fread("/oak/stanford/groups/horence/julias/splash_postprocessing/granny_analysis/cells_with_granny_transcripts.tsv",header=F,fill=T)
metadata = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/10X_runs/Sponge_SRP216435_anchor_count_3/metadata_modified.tsv")


umitools_directory = "/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Sponge_10X_SRP216435/"
whitelist_files = list.files(path=umitools_directory,pattern = "_whitelist_Granulocyte_minimum_50umis.txt",recursive = TRUE)
cell_barcodes = data.table()
for (counter in 1:length(whitelist_files)){
  whitelist_file_counter = paste(umitools_directory,whitelist_files[counter],sep="")
  fastq_file_name = strsplit(whitelist_files[counter],split="_")[[1]][1]
  whitelist_counter = fread(whitelist_file_counter,select=c("V1","V3"))
  names(whitelist_counter) = c("barcode","num_UMI")
  whitelist_counter[,fastq_file:= fastq_file_name]
  cell_barcodes = rbind(whitelist_counter,cell_barcodes) 
}

cell_barcodes[fastq_file%in%c("SRR9841059","SRR9841060","SRR9841061","SRR9841062","SRR9841063","SRR9841064","SRR9841065","SRR9841066"),sample:="S1"]
cell_barcodes[fastq_file%in%c("SRR9841067","SRR9841068","SRR9841069","SRR9841070","SRR9841071","SRR9841072","SRR9841073","SRR9841074"),sample:="S2"]
cell_barcodes[fastq_file%in%c("SRR9841075","SRR9841076","SRR9841077","SRR9841078","SRR9841079","SRR9841080","SRR9841081","SRR9841082"),sample:="S3"]
cell_barcodes[fastq_file%in%c("SRR9841083","SRR9841084","SRR9841085","SRR9841086"),sample:="S4"]
cell_barcodes[,num_UMI_per_cell:=sum(num_UMI),by=paste(barcode,sample)]
cell_barcodes = unique(cell_barcodes[,list(barcode,sample,num_UMI_per_cell)])


f[,UMI:=strsplit(V1,split="_")[[1]][7],by=V1]
f[,cell_barcode:=strsplit(V1,split="_")[[1]][6],by=V1]
f[,fastq_file_name:=strsplit(V1,split=":")[[1]][2],by=V1]
f[,fastq_file_name:=strsplit(fastq_file_name,split=".",fixed = T)[[1]][1],by=fastq_file_name]
f[fastq_file_name%in%c("SRR9841059","SRR9841060","SRR9841061","SRR9841062","SRR9841063","SRR9841064","SRR9841065","SRR9841066"),sample:="S1"]
f[fastq_file_name%in%c("SRR9841067","SRR9841068","SRR9841069","SRR9841070","SRR9841071","SRR9841072","SRR9841073","SRR9841074"),sample:="S2"]
f[fastq_file_name%in%c("SRR9841075","SRR9841076","SRR9841077","SRR9841078","SRR9841079","SRR9841080","SRR9841081","SRR9841082"),sample:="S3"]
f[fastq_file_name%in%c("SRR9841083","SRR9841084","SRR9841085","SRR9841086"),sample:="S4"]
f1 = f[!duplicated(paste(V3,sample,cell_barcode,UMI))]
f1 = f1[,list(sample,cell_barcode,UMI,fastq_file_name,V3,V12)]
f1[,num_granny_transcripts_per_cell:=.N,by=paste(sample,cell_barcode)]

f1[V3=="Group5_Jpb124_cDNA",Group5_count_per_cell:=.N,by=paste(sample,cell_barcode)]
f1[V3=="Group4_Jpb134_cDNA",Group4_count_per_cell:=.N,by=paste(sample,cell_barcode)]
f1[V3=="Group3_Jpb72_cDNA",Group3_count_per_cell:=.N,by=paste(sample,cell_barcode)]
f1[V3=="Group2_Jpb100_cDNA",Group2_count_per_cell:=.N,by=paste(sample,cell_barcode)]
f1[V3=="Group1_Jpb102_cDNA",Group1_count_per_cell:=.N,by=paste(sample,cell_barcode)]
f1 = merge(f1,cell_barcodes,all.x=T,all.y=F,by.x=c("sample","cell_barcode"),by.y=c("sample","barcode")) ## now I add the total number of UMI per cell

f1[,Group4_count_per_cell:=max(Group4_count_per_cell,na.rm = T),by=paste(sample,cell_barcode)]
f1[,Group5_count_per_cell:=max(Group5_count_per_cell,na.rm = T),by=paste(sample,cell_barcode)]
f1[,Group3_count_per_cell:=max(Group3_count_per_cell,na.rm = T),by=paste(sample,cell_barcode)]
f1[,Group2_count_per_cell:=max(Group2_count_per_cell,na.rm = T),by=paste(sample,cell_barcode)]
f1[,Group1_count_per_cell:=max(Group1_count_per_cell,na.rm = T),by=paste(sample,cell_barcode)]

f1[is.na(Group5_count_per_cell),Group5_count_per_cell:=0]
f1[is.na(Group4_count_per_cell),Group4_count_per_cell:=0]
f1[is.na(Group3_count_per_cell),Group3_count_per_cell:=0]
f1[is.na(Group2_count_per_cell),Group2_count_per_cell:=0]
f1[is.na(Group1_count_per_cell),Group1_count_per_cell:=0]

setnames(f1,"num_granny_transcripts_per_cell","Granny_count_per_cell")

#### now merging with metadata file ########
metadata[,sample_cell:=paste(sample_name,cell)]
f1[,sample_cell:=paste(sample,cell_barcode)]
f1 = merge(f1,metadata[,list(sample_cell,cell_type)],all.x=T,all.y=F,by.x="sample_cell",by.y="sample_cell")


## now I want to make a dot plot that shows the cell_type (color), top transcript (shape), read count (y-axis), and fraction of top transcript (size of dot)
f1_uniq = f1[!duplicated(sample_cell)]
f1_uniq[,num_cell_per_celltype:=length(unique(sample_cell)),by=cell_type]
f1_uniq=f1_uniq[num_cell_per_celltype>30]
f1_uniq = f1_uniq[Granny_count_per_cell>4]
f1_uniq_only_transcript_counts = f1_uniq[,list(Group1_count_per_cell,Group2_count_per_cell,Group3_count_per_cell,Group4_count_per_cell,Group5_count_per_cell)]
f1_uniq_only_transcript_counts[, max:= do.call(pmax, .SD)]
f1_uniq_only_transcript_counts[Group1_count_per_cell==max,top_transcript:="Group1"]
f1_uniq_only_transcript_counts[Group2_count_per_cell==max,top_transcript:="Group2"]
f1_uniq_only_transcript_counts[Group3_count_per_cell==max,top_transcript:="Group3"]
f1_uniq_only_transcript_counts[Group4_count_per_cell==max,top_transcript:="Group4"]
f1_uniq_only_transcript_counts[Group5_count_per_cell==max,top_transcript:="Group5"]
f1_uniq_only_transcript_counts[,tot_count:=Group1_count_per_cell+Group2_count_per_cell+Group3_count_per_cell+Group4_count_per_cell+Group5_count_per_cell,by=1:nrow(f1_uniq_only_transcript_counts)]
f1_uniq_only_transcript_counts[,frac_top_transcript:=max/tot_count]
f1_uniq = cbind(f1_uniq,f1_uniq_only_transcript_counts[,list(top_transcript,max,frac_top_transcript)])
f1_uniq[is.na(cell_type),cell_type:="NA"]
f1_uniq = setorder(f1_uniq,-num_cell_per_celltype)
f1_uniq[,normalized_granny_expression:=Granny_count_per_cell/num_UMI_per_cell]
f1_uniq$cell_type = factor(f1_uniq$cell_type,levels = unique(f1_uniq$cell_type))
f1_uniq = setorder(f1_uniq,cell_type,-normalized_granny_expression)
f1_uniq$sample_cell = factor(f1_uniq$sample_cell,levels = f1_uniq$sample_cell)
ggplot(f1_uniq[num_cell_per_celltype>2 & cell_type!="NA" & cell_type!="7"], aes(x=sample_cell, y=Granny_count_per_cell,  color=top_transcript,shape = cell_type,size = frac_top_transcript)) +     geom_point() + theme_bw()  +scale_color_manual(values=brewer.pal(length(unique(f1_uniq$top_transcript)), 'Paired'))
ggplot(f1_uniq[num_cell_per_celltype>6 & cell_type!="NA" & cell_type!="7"], aes(x=sample_cell, y=normalized_granny_expression,  color=top_transcript,shape = cell_type,size = frac_top_transcript)) +     geom_point() + theme_bw()  +scale_color_manual(values=brewer.pal(length(unique(f1_uniq$top_transcript)), 'Paired'))
