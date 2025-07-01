library(Seurat)
library(rhdf5)
library(data.table)

## as there is no cell annotation for sample SRR9051002 in metadata for alte tail bud I only use samples SRR9051003 and SRR9051004

file_path_1005 <- "/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/SRP198321/GSM3764784_SRR9051005_larva1_raw_gene_bc_matrices_h5.h5"
file_path_1006 <- "/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/SRP198321/GSM3764785_SRR9051006_larva2_raw_gene_bc_matrices_h5.h5"
file_path_1007 <- "/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/SRP198321/GSM3764786_SRR9051007_larva3_raw_gene_bc_matrices_h5.h5"

#### below I want to get the counts for the ayelet anchor ##########################
ciona_1002 = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_SRP198321/SRR9051002/Ayelet_seq_satc/satcOut112.tsv",header = F)
ciona_1003 = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_SRP198321/SRR9051003/Ayelet_seq_satc/satcOut112.tsv",header = F)
ciona_1004 = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_SRP198321/SRR9051004/Ayelet_seq_satc/satcOut112.tsv",header = F)
ciona_1005 = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_SRP198321/SRR9051005/Ayelet_seq_satc/satcOut112.tsv",header = F)
ciona_1006 = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_SRP198321/SRR9051006/Ayelet_seq_satc/satcOut112.tsv",header = F)
ciona_1007 = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_SRP198321/SRR9051007/Ayelet_seq_satc/satcOut112.tsv",header = F)

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
ciona_all[,cell_barcode:=paste(cell_barcode,"-1",sep=""),by=cell_barcode]



metadata = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/utility_files/metadata/non_model_organisms/Ciona_intestinalis_SRP198321/ciona_metadata_modified.txt")
metadata[,sample_new:=sample]
metadata[sample=="SRR9051005",sample_new:="SRR9051007"]
metadata[sample=="SRR9051006",sample_new:="SRR9051005"]
metadata[sample=="SRR9051007",sample_new:="SRR9051006"]
metadata[,sample:=NULL]
setnames(metadata,"sample_new","sample")
setnames(metadata,"cell","cell_barcode")
metadata[,cell_barcode:=paste(cell_barcode,"-1",sep=""),by=cell_barcode]
metadata[,cell:=paste(sample,cell_barcode,sep="_"),by=1:nrow(metadata)]
metadata = metadata[sample%in%c("SRR9051002","SRR9051003","SRR9051004","SRR9051005","SRR9051006","SRR9051007")]
metadata[sample%in%c("SRR9051002","SRR9051003","SRR9051004"),dev_stage:="late_tail_bud_II"]
metadata[sample%in%c("SRR9051005","SRR9051006","SRR9051007"),dev_stage:="larva"]
metadata[,cell_type_dev:=paste(cell_type,dev_stage,sep="_"),by=paste(cell_type,dev_stage)]

## merging with metadata
ciona_all[,cell:=paste(sample,cell_barcode,sep="_"),by=paste(sample,cell_barcode)]
ciona_all = merge(ciona_all,metadata[,list(cell,cell_type)],all.x=T,all.y=F,by.x="cell",by.y="cell")
ciona_all[,cell_type_dev:=paste(cell_type,dev_stage,sep="_"),by=paste(cell_type,dev_stage)]
#######################################################################################

# Read datasets from the HDF5 file
data_1005 = Read10X_h5(file_path_1005)
seurat_object_1005 <- CreateSeuratObject(counts = data_1005)
data_1006 = Read10X_h5(file_path_1006)
seurat_object_1006 <- CreateSeuratObject(counts = data_1006)
data_1007 = Read10X_h5(file_path_1007)
seurat_object_1007 <- CreateSeuratObject(counts = data_1007)

## to get the number of cells expressing a gene
sum(GetAssayData(object = seurat_object, slot = "counts")["KH2012:KH.L18.59",]>0)

###############################################################################
##### first I want to look at the unannotated cells in sample SRR9051005 ######
###############################################################################
anchor_negative_barcodes = metadata[sample=="SRR9051005"]$cell_barcode ## this is the list of all annotated cell barcodes that I want to use as control set
#all_barcodes = metadata[sample=="SRR9051007" & cell_type=="32_mesenchyme"]$cell_barcode ## this is the list of all cell barcodes that are from cluster 18_mesenchyme
anchor_negative_barcodes = anchor_negative_barcodes[which(!anchor_negative_barcodes%in%unique(ciona_all[sample=="SRR9051005"]$cell_barcode))] #removing those annotated barcodes that are anchor+ 

anchor_positive_barcodes = unique(ciona_all[sample=="SRR9051005" & is.na(cell_type)]$cell_barcode) # unannotated cells containing the anchor

all_barcodes = c(anchor_negative_barcodes,anchor_positive_barcodes)
barcodes_from_seurat_object=rownames(seurat_object@meta.data)
barcodes_from_seurat_object = data.table(barcodes_from_seurat_object) # I get these barcodes to make sure that all barcodes from all_barcodes and anchor_positive_barcodes are in the seurat object
names(barcodes_from_seurat_object) ="V1"
all_barcodes= all_barcodes[which(all_barcodes%in%barcodes_from_seurat_object$V1)]

## as there might still be duplicated barcodes because of ge well id, I remove all those from the list of all barcodes
duplicated_barcodes = all_barcodes[which(duplicated(all_barcodes))]
all_barcodes = all_barcodes[which(!all_barcodes%in%duplicated_barcodes)]

seurat_object_filtered = subset(seurat_object_1005, cells = all_barcodes)

###### Subset the sparse matrix to cells from 18_mesenchyme #############
#### now I manually define clusters based on the presence or absence of the anchor ##################
manual_clusters = rep("cluster_0", length(all_barcodes))                        # cluster_0 means cells w/o anchor
manual_clusters[ which(all_barcodes%in%anchor_positive_barcodes)] = "cluster_1" # cluster_1 means cells with anchor
names(manual_clusters) = all_barcodes
table(manual_clusters)
manual_clusters <- manual_clusters[match(rownames(seurat_object_filtered@meta.data), names(manual_clusters))] # I sort the manual cluster vector so that the same columns in clusters vector and sparse matrix correspond to the same cell barcode
seurat_object_filtered$manual_clusters <- manual_clusters
seurat_object_filtered <- SetIdent(seurat_object_filtered, value = manual_clusters)
##################################################

# function to log normalize, find variable features, and scale the data
normVarScaleData<-function(seuratOb,numFeat=5000){
  seuratOb <- NormalizeData(seuratOb)
  seuratOb <- FindVariableFeatures(seuratOb, selection.method = "vst", nfeatures=numFeat)
  seuratOb <- ScaleData(seuratOb, features = rownames(seuratOb))
  return(seuratOb)
}
seurat_object_filtered <- normVarScaleData(seurat_object_filtered) #normalization step

# Check metadata
head(seurat_object_filtered@meta.data)
table(seurat_object_filtered$manual_clusters)  # Verify manual clusters


###### Differential expression between cluster_0 and cluster_1 #####
markers <- FindMarkers(
  object = seurat_object_filtered,
  ident.1 = "cluster_1",
  ident.2 = "cluster_0",
  group.by = "manual_clusters",  # Specify the custom cluster column
  logfc.threshold = 0.25,
  min.pct = 0.25,
  test.use = "wilcox"
)

markers$gene <- rownames(markers)
markers$Significant <- ifelse(
  markers$p_val_adj < 0.05 & abs(markers$avg_log2FC) > 0.25,
  "Yes",
  "No"
)
markers = data.table(markers)
markers = markers[Significant=="Yes"]
markers = setorder(markers, avg_log2FC)

gene_description = fread("/oak/stanford/groups/horence/SPLASH_reference_fastas/eukaryote/non_model_organsisms/Ciona_intestinalis_sea_squirt/KH_gene_id_description.txt",sep="\t")
markers[,gene1:=gsub("KH2012:","",gene),by=gene]
markers = merge(markers,gene_description,all.x=T,all.y=F,by.x="gene1",by.y="gene")
markers [,gene1:=NULL]
write.table(markers,"/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/unannotated_larva_SRR9051006_cry11bb_marker_genes.tsv",sep="\t",row.names=F,quote=F)

VlnPlot(seurat_object, features = c("KH2012:KH.L171.6"), slot = "data", log = FALSE)