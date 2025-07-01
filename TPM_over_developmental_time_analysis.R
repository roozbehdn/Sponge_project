library(data.table)
library(ggplot2)

f = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/kallisto_runs/aggregated_abundances.tsv")


## below I want to make the expression plot based on bowtie alignment to the cDNAs Peter has assembled from D12 data
f= fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/alignment_to_peter_GranRep_cDNAs_IllumD12cd/aggregated_normal_immune.sam",header=F,fill=T)
f[,sample_name:=strsplit(V1,split=":")[[1]][1],by=V1]
f[sample_name%like%"D5",Day:=5]
f[sample_name%like%"D8",Day:=8]
f[sample_name%like%"D12",Day:=12]
f[sample_name%like%"developmental",treatment:="normal"]
f[sample_name%like%"LPS",treatment:="LPS"]
f[sample_name%like%"CD",treatment:="CD"]

f[,allelic_count_per_sample:=.N,by=paste(V3,sample_name)]
f_uniq = unique(f[,list(V3,Day,treatment,allelic_count_per_sample)])

## now below I add the sequencing_depth to compute the normalized expression
num_reads_per_sample  = c(150445421, 166068975, 159671742, 191772289, 213922402, 178972205, 172084785)
sample = c(normal_D5, normal_D8, normal_D12, D8_CD, D8_LPS, D12_CD, D12_LPS)

f_uniq[Day==5 & treatment=="normal",seq_depth:=150445421]
f_uniq[Day==8 & treatment=="normal",seq_depth:=166068975]
f_uniq[Day==12 & treatment=="normal",seq_depth:=159671742]
f_uniq[Day==8 & treatment=="CD",seq_depth:=191772289]
f_uniq[Day==8 & treatment=="LPS",seq_depth:=213922402]
f_uniq[Day==12 & treatment=="CD",seq_depth:=178972205]
f_uniq[Day==12 & treatment=="LPS",seq_depth:=172084785]
f_uniq[,gene:=strsplit(V3,split="_")[[1]][1],by=V3]
f_uniq[,gene_count_per_sample:=sum(allelic_count_per_sample),by=paste(Day,treatment,gene)]

f_uniq[,allelic_ratio:=allelic_count_per_sample/gene_count_per_sample]
f_uniq[,normalized_count:=gene_count_per_sample/seq_depth*10^5]
f_uniq_gene = unique(f_uniq[!duplicated(paste(gene,Day,treatment))])
f_uniq_gene$Day = factor(f_uniq_gene$Day,levels=c(5,8,12))
f_uniq_gene$treatment = factor(f_uniq_gene$treatment,levels=c("CD","LPS","normal"))
f_uniq_gene$gene = factor(f_uniq_gene$gene,levels=c("GranRep1","GranRep2","GranRep3","GranRep5"))

ggplot(f_uniq_gene, aes(x = Day, y = normalized_count, fill = gene)) +
  facet_wrap(~treatment) + theme_bw() + geom_point()

#### now I want to compute the fold change between control and treatment
f_uniq_gene_control = f_uniq_gene[treatment == "normal"]  
f_uniq_gene = merge(f_uniq_gene,unique(f_uniq_gene_control[,list(gene,Day,normalized_count)]),all.x=T,all.y=F,by.x=c("Day","gene"),by.y=c("Day","gene"))
setnames(f_uniq_gene,c("normalized_count.x","normalized_count.y"),c("normalized_count","normalized_count_control"))
f_uniq_gene[,fold_change:=normalized_count/normalized_count_control]


#### now I want to compute the total granny expression ##############
f_uniq_gene[,total_Granny_expression:=sum(normalized_count),by=paste(Day,treatment)]
f_uniq_gene_total_expression = unique(f_uniq_gene[,list(total_Granny_expression,Day,treatment)])
f_uniq_gene_total_expression_control = f_uniq_gene_total_expression[treatment == "normal"] 
f_uniq_gene_total_expression  = merge(f_uniq_gene_total_expression, f_uniq_gene_total_expression_control[,list(Day,total_Granny_expression)],all.x=T,all.y=F,by.x=c("Day"),by.y=c("Day"))
setnames(f_uniq_gene_total_expression,c("total_Granny_expression.x","total_Granny_expression.y"),c("total_Granny_expression","total_Granny_expression_control"))
f_uniq_gene_total_expression[,fold_change:=total_Granny_expression/total_Granny_expression_control]

ggplot(f_uniq_gene_total_expression, aes(x = Day, y = total_Granny_expression)) +  geom_point()+ facet_wrap(~treatment) + theme_bw()


#####################################################################################################################
########### Now I want to compute the TPMs for granulcoyte markers in RNA immune/developmental data #################
f= fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/alignment_to_granulocyte_markers/aggregated_normal_immune.sam",header=F,fill=T)
f[,sample_name:=strsplit(V1,split=":")[[1]][1],by=V1]
f[sample_name%like%"D5",Day:=5]
f[sample_name%like%"D8",Day:=8]
f[sample_name%like%"D12",Day:=12]
f[sample_name%like%"developmental",treatment:="normal"]
f[sample_name%like%"LPS",treatment:="LPS"]
f[sample_name%like%"CD",treatment:="CD"]
f[,gene_name1:=strsplit(V3,split="_")[[1]][1],by=V3]
f[,gene_name2:=strsplit(V3,split="_")[[1]][2],by=V3]
f$gene_name=paste(f$gene_name1,f$gene_name2,sep="_") ## gene name for c104489_g2_i1.mrna5 is c104489_g2

f[,gene_count_per_sample:=.N,by=paste(gene_name,sample_name)]
f_uniq = unique(f[,list(gene_name,Day,treatment,gene_count_per_sample)])

f_uniq[Day==5 & treatment=="normal",seq_depth:=150445421]
f_uniq[Day==8 & treatment=="normal",seq_depth:=166068975]
f_uniq[Day==12 & treatment=="normal",seq_depth:=159671742]
f_uniq[Day==8 & treatment=="CD",seq_depth:=191772289]
f_uniq[Day==8 & treatment=="LPS",seq_depth:=213922402]
f_uniq[Day==12 & treatment=="CD",seq_depth:=178972205]
f_uniq[Day==12 & treatment=="LPS",seq_depth:=172084785]

f_uniq[,normalized_count:=gene_count_per_sample/seq_depth*10^5]


sponge_marker_genes = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/utility_files/metadata/non_model_organisms/Sponge/Sponge_marker_gene_clusters.csv",sep=",")
names(sponge_marker_genes)[3] = "gene_name"
sponge_marker_genes[,gene_name_1:=strsplit(gene_name,split=" ")[[1]][1],by=gene_name]
sponge_marker_genes = sponge_marker_genes[!is.na(`Cell Type`)]
sponge_marker_genes[,num_celltype_per_marker:=length(unique(`Cell Type`)),by=gene_name_1]
sponge_granulocyte_exclusive_markers = sponge_marker_genes[num_celltype_per_marker==1 & `Cell Type` %like%"Granulocytes"]
sponge_granulocyte_exclusive_markers[,gene_name_1:=gsub("-","_",gene_name_1),by=gene_name_1]
f_uniq = f_uniq[gene_name%in%sponge_granulocyte_exclusive_markers$gene_name_1]

f_uniq_12 = f_uniq[Day==12]
f_uniq_12_wide <- reshape(f_uniq_12[,list(normalized_count,treatment,gene_name)], 
                          idvar = "gene_name", 
                          timevar = "treatment", 
                          direction = "wide")


f_uniq_8 = f_uniq[Day==8]
f_uniq_8_wide <- reshape(f_uniq_8[,list(normalized_count,treatment,gene_name)], 
                          idvar = "gene_name", 
                          timevar = "treatment", 
                          direction = "wide")

f_uniq_8_wide[,fold_change_CD:=normalized_count.CD/normalized_count.normal]
f_uniq_8_wide[,fold_change_LPS:=normalized_count.LPS/normalized_count.normal]
f_uniq_8_wide = f_uniq_8_wide[!is.na(normalized_count.normal)]
f_uniq_8_wide = f_uniq_8_wide[normalized_count.normal>0.05]
f_uniq_12_wide[,fold_change_LPS:=normalized_count.LPS/normalized_count.normal]
f_uniq_12_wide[,fold_change_CD:=normalized_count.CD/normalized_count.normal]
f_uniq_12_wide = f_uniq_12_wide[normalized_count.normal>0.05]

f_uniq_12_wide = f_uniq_12_wide[normalized_count.normal < 1]
f_uniq_8_wide = f_uniq_8_wide[gene_name%in%f_uniq_12_wide$gene_name]
f_uniq_12_wide = f_uniq_12_wide[gene_name%in%f_uniq_8_wide$gene_name]


data_table_to_plot = cbind(f_uniq_8_wide[,list(fold_change_LPS,fold_change_CD)],f_uniq_12_wide[,list(fold_change_LPS,fold_change_CD)])
names(data_table_to_plot ) = c("LPS_8","CD_8","LPS_12","CD_12")
data_table_to_plot[,gene_name:=1:12] 
dt_long <- melt(data_table_to_plot, 
                id.vars = "gene_name",   # Column to keep fixed
                variable.name = "condition",  # New column for melted variable names
                value.name = "Value")  # New column for values
ggplot(dt_long, aes(x=condition, y=Value)) +  geom_boxplot() + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + geom_jitter()