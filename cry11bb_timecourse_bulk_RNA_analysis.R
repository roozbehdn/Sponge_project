library(data.table)

sorted_dev_stages = c("stage 1 (fertilized egg)", "stage 2 (2 cell stage)", "stage 4 (8 cell stage)", "stage 5 (16 cell stage)", "stage 6 (32 cell stage)", "stage 8 (64 cell stage)", "stage 10 (initial gastrula)", "stage 12 (mid gastrula)", "stage 14 (early Neurula)", "stage 16 (late Neurula)", "stage 19 (early Tailbud)", "stage 22 (mid Tailbud)", "stage 24 (late Tailbud)", "stage 27 (early swimming larva)", "stage 29 (late swimming larva)", "stage 35 (early rotation)", "stage 37 (late rotation)", "stage 38 (early juvenile I)", "stage 40 (mid juvenile)", "stage Late Juvenile", "stage Adult")

satc_count = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_DRP003810/Ayelet_seq_satc/Ayelet_anchor_count.tsv",header=F)

metadata = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/bulk_RNA/DRP003810/Metadata_DRP003810.csv", sep=",", header=T, select = c("Run", "AvgSpotLen", "Bases", "dev_stage", "Organism"))
metadata = metadata[Organism=="Ciona intestinalis"]
metadata[,num_reads:=Bases/AvgSpotLen,by=Run]

sample_name_to_id_conversion = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_DRP003810/sample_name_to_id.mapping.txt",header=F)


names(satc_count) = c("sample_id","anchor","target","count")
names(sample_name_to_id_conversion) = c("Run","sample_id")
satc_count = merge(satc_count,sample_name_to_id_conversion,all.x=T,all.y=F,by.x="sample_id",by.y="sample_id")

## below I manually add the counts for those samples that have 0 or fewer than anchor count threshold and therefore are not in splash satc file
new = data.table(Run=c("DRR032683", "DRR032687", "DRR032694", "DRR032695", "DRR032705", "DRR032707","DRR032717", "DRR032718", "DRR032719", "DRR032721", "DRR032722", "DRR032724", "DRR032725"),count=c(1,2,0,0,0,1,0,0,0,3,0,0,0))
satc_count = rbind(satc_count,new,fill=T)
satc_count = merge(satc_count,metadata[,list(Run,num_reads,dev_stage)],all.x=T,all.y=F,by.x="Run",by.y="Run")
satc_count[,anchor_count_per_sample:=sum(count),by=Run]
satc_count[,normalized_count:=anchor_count_per_sample/num_reads*(10^5),by=Run] 

count_uniq_sample = satc_count[!duplicated(Run),list(Run,dev_stage,normalized_count,anchor_count_per_sample)]
count_uniq_sample$dev_stage = factor(count_uniq_sample$dev_stage,levels=sorted_dev_stages)

## below I make the plot for the normalized expression per developmental stage
library(dplyr)
df.summary <- count_uniq_sample %>%
  group_by(dev_stage) %>%
  summarise(
    sd = sd(normalized_count, na.rm = TRUE),
    normalized_count = mean(normalized_count)
  )

ggplot(count_uniq_sample, aes(dev_stage, normalized_count)) +  geom_jitter(position = position_jitter(0.2), color = "darksalmon") + theme_bw() + geom_pointrange(aes(ymin = normalized_count-sd, ymax = normalized_count+sd),data = df.summary) +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))



############## NOW I DO THIS FOR THE OTHER DATASET SRP339256 ####################
sorted_dev_stages = c("prehatched larvae", "just after hatched larvae", "4 hours post hatched larvae", "6 hours post hatched larvae", "12 hours post hatched larvae", "larvae at tail regression completed stage")

satc_count = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_SRP339256/Ayelet_seq_satc/Ayelet_anchor_count.tsv",header=F)

metadata = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/bulk_RNA/SRP339256/Metadata_SRP339256.csv", sep=",", header=T, select = c("Run", "AvgSpotLen", "Bases", "dev_stage", "Organism"))
metadata[,num_reads:=Bases/AvgSpotLen,by=Run]

sample_name_to_id_conversion = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_SRP339256/sample_name_to_id.mapping.txt",header=F)


names(satc_count) = c("sample_id","anchor","target","count")
names(sample_name_to_id_conversion) = c("Run","sample_id")
satc_count = merge(satc_count,sample_name_to_id_conversion,all.x=T,all.y=F,by.x="sample_id",by.y="sample_id")
satc_count[,Run:=strsplit(Run,split="_")[[1]][1],by=Run] # I do this to collapse across R1 and R2 for each sample

satc_count = merge(satc_count,metadata[,list(Run,num_reads,dev_stage)],all.x=T,all.y=F,by.x="Run",by.y="Run")
satc_count[,anchor_count_per_sample:=sum(count),by=Run]
satc_count[,normalized_count:=anchor_count_per_sample/num_reads*(10^5),by=Run] 

count_uniq_sample = satc_count[!duplicated(Run),list(Run,dev_stage,normalized_count,anchor_count_per_sample)]
count_uniq_sample$dev_stage = factor(count_uniq_sample$dev_stage,levels=sorted_dev_stages)

## below I make the plot for the normalized expression per developmental stage
library(dplyr)
df.summary <- count_uniq_sample %>%
  group_by(dev_stage) %>%
  summarise(
    sd = sd(normalized_count, na.rm = TRUE),
    normalized_count = mean(normalized_count)
  )

ggplot(count_uniq_sample, aes(dev_stage, normalized_count)) +  geom_jitter(position = position_jitter(0.2), color = "darksalmon") + theme_bw() + geom_pointrange(aes(ymin = normalized_count-sd, ymax = normalized_count+sd),data = df.summary) +  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
