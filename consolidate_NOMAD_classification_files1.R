library(data.table)

directory = "/oak/stanford/groups/horence/Roozbeh/Sponge_project/runs/10X_runs/Mouse_lemur_10x"

list_files = list.files(path=directory,pattern = "compactors_generic_repeat.tsv",recursive = TRUE)

nomad = data.table()
for (counter in 1:length(list_files)){
  tissue_counter = strsplit(list_files[counter],split= "/")[[1]][2] 
  donor_counter = strsplit(list_files[counter],split= "/")[[1]][1] 
  nomad_counter = fread(paste(directory,list_files[counter],sep="/"))
  nomad_counter[,tissue:=tissue_counter]
  nomad_counter[,donor:=donor_counter]
  nomad = rbind(nomad,nomad_counter,fill = TRUE)
  print(counter)
}

#nomad = rbind(nomad_TSP10_Blood,nomad_TSP10_Fat,nomad_TSP10_Skin,nomad_TSP11,nomad_TSP12,nomad_TSP14_Blood,nomad_TSP14_BoneMarrow,nomad_TSP14_Fat,nomad_TSP14_Lung,nomad_TSP14_LymphNode,nomad_TSP14_Muscle,nomad_TSP14_Spleen,nomad_TSP14_Vasculature,nomad_TSP7_Blood,nomad_TSP7_LymphNode,nomad_TSP7_SalivaryGland,nomad_TSP7_Spleen,nomad_TSP7_Tongue)
#nomad_splicing_corrected = nomad[anchor_event=="splicing" & pval_rand_init_alt_max_corrected < 0.05 & effect_size_bin>0.1]
#nomad_splicing_uncorrected = nomad[anchor_event=="splicing" & pval_rand_init_alt_max < 0.05 & effect_size_bin>0.1]

#listInput_upsetplot = list( TSP7_Blood = nomad_splicing[donor=="TSP7" & tissue == "Blood"]$anchor_gene,TSP10_Blood = nomad_splicing[donor=="TSP10" & tissue == "Blood"]$anchor_gene,TSP14_Blood = nomad_splicing[donor=="TSP14" & tissue == "Blood"]$anchor_gene)
#upset(fromList(listInput_upsetplot),sets = c("TSP7_Blood","TSP10_Blood","TSP14_Blood"))

#listInput_upsetplot = list( TSP10_Fat = nomad_splicing[donor=="TSP10" & tissue == "Fat"]$anchor_gene,TSP14_Fat = nomad_splicing[donor=="TSP14" & tissue == "Fat"]$anchor_gene)
#upset(fromList(listInput_upsetplot),sets = c("TSP10_Fat","TSP14_Fat"))

write.table(nomad,"/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/Lemur_10x_compactors_generic_repeat.tsv",sep="\t",row.names = FALSE,quote = FALSE)
