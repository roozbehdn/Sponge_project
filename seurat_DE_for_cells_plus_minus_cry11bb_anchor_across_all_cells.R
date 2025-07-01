library(Seurat)
library(rhdf5)
library(data.table)

process_h5_file <- function(file_path) {
  
  # Read datasets from the HDF5 file
  barcodes <- h5read(file_path, "/ci_kh_ghostdb_mchsv40_cfpsv40_ensembl_mt/barcodes")
  barcodes <- gsub("-1", "", barcodes)
  gene_names <- h5read(file_path, "/ci_kh_ghostdb_mchsv40_cfpsv40_ensembl_mt/gene_names")
  data <- h5read(file_path, "/ci_kh_ghostdb_mchsv40_cfpsv40_ensembl_mt/data")
  indices <- h5read(file_path, "/ci_kh_ghostdb_mchsv40_cfpsv40_ensembl_mt/indices")
  indptr <- h5read(file_path, "/ci_kh_ghostdb_mchsv40_cfpsv40_ensembl_mt/indptr")
  shape <- h5read(file_path, "/ci_kh_ghostdb_mchsv40_cfpsv40_ensembl_mt/shape")
  
  # Convert to a sparse matrix
  sparse_matrix <- sparseMatrix(
    i = indices + 1,        # Convert 0-based to 1-based indexing
    p = indptr,
    x = data,
    dims = shape
  )
  
  # Add row and column names
  rownames(sparse_matrix) <- gene_names
  colnames(sparse_matrix) <- barcodes
  return(sparse_matrix)
}


file_path_1002 <- "/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/SRP198321/GSM3764781_SRR9051002_latTII1a_raw_gene_bc_matrices_h5.h5"
file_path_1003 <- "/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/SRP198321/GSM3764782_SRR9051003_latTII1b_raw_gene_bc_matrices_h5.h5"
file_path_1004 <- "/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/SRP198321/GSM3764783_SRR9051004_latTII2_raw_gene_bc_matrices_h5.h5"
file_path_1005 <- "/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/SRP198321/GSM3764784_SRR9051005_larva1_raw_gene_bc_matrices_h5.h5"
file_path_1006 <- "/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/SRP198321/GSM3764785_SRR9051006_larva2_raw_gene_bc_matrices_h5.h5"
file_path_1007 <- "/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/SRP198321/GSM3764786_SRR9051007_larva3_raw_gene_bc_matrices_h5.h5"

file_paths = c(file_path_1002, file_path_1003, file_path_1004, file_path_1005, file_path_1006, file_path_1007)

##############################################################################
#### below I want to get the counts for the ayelet anchor ####################
##############################################################################
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

#######################################
######### merging with metadata #######
#######################################
metadata = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/utility_files/metadata/non_model_organisms/Ciona_intestinalis_SRP198321/ciona_metadata_modified.txt")
setnames(metadata,"cell","cell_barcode")
metadata[,cell:=paste(sample,cell_barcode,sep="_"),by=1:nrow(metadata)]
metadata = metadata[sample%in%c("SRR9051002","SRR9051003","SRR9051004","SRR9051005","SRR9051006","SRR9051007")]
metadata[sample%in%c("SRR9051002","SRR9051003","SRR9051004"),dev_stage:="late_tail_bud_II"]
metadata[sample%in%c("SRR9051005","SRR9051006","SRR9051007"),dev_stage:="larva"]
metadata[,cell_type_dev:=paste(cell_type,dev_stage,sep="_"),by=paste(cell_type,dev_stage)]

ciona_all[,cell:=paste(sample,cell_barcode,sep="_"),by=paste(sample,cell_barcode)]
ciona_all = merge(ciona_all,metadata[,list(cell,cell_type)],all.x=T,all.y=F,by.x="cell",by.y="cell")
ciona_all[,cell_type_dev:=paste(cell_type,dev_stage,sep="_"),by=paste(cell_type,dev_stage)]
######################################
######################################
######################################



#################################################
############# Create a Seurat object ############
#################################################
sparse_matrices <- lapply(file_paths, process_h5_file)
combined_matrix <- do.call(cbind, sparse_matrices)

colnames(sparse_matrix)=1:length(barcodes)
rownames(sparse_matrix)=1:length(gene_names)
seurat_object <- CreateSeuratObject(counts = sparse_matrix)


#################################################################################################################
##### first I want to look at the cells from cluster 18_mesenchyme +/- cry11bb anchor in sample SRR9051003 ######
#################################################################################################################
all_barcodes = metadata[sample=="SRR9051003" & cell_type=="18_mesenchyme"]$cell_barcode ## this is the list of all cell barcodes that are from cluster 18_mesenchyme
anchor_positive_barcodes = unique(ciona_all[sample=="SRR9051003" & cell_type=="18_mesenchyme"]$cell_barcode)

###### Subset the sparse matrix to cells from 18_mesenchyme #############
# Check which columns match
columns_to_keep <- colnames(sparse_matrix) %in% all_barcodes

filtered_matrix <- sparse_matrix[, columns_to_keep]
#########################################################################

#### now I manually define clusters based on the presence or absence of the anchor ##################
manual_clusters = rep("cluster_0", length(all_barcodes))                        # cluster_0 means cells w/o anchor 
manual_clusters[ which(all_barcodes%in%anchor_positive_barcodes)] = "cluster_1" # cluster_1 means cells with anchor

names(manual_clusters) = all_barcodes

manual_clusters <- manual_clusters[match(colnames(filtered_matrix), names(manual_clusters))] # I sort the manula cluster vector so that the same columns in clusters vector and sparse matrix correspond to the same cell barcode

all_barcodes_sorted = names(manual_clusters)

## now that I know both cluster vector and filtered matrix are sorted based on the same set of cell barcodes i can just renames them based on numbers 
names(manual_clusters) = 1:length(all_barcodes)
colnames(filtered_matrix) = 1:length(all_barcodes)
rownames(filtered_matrix)=1:length(gene_names)

seurat_object <- CreateSeuratObject(counts = filtered_matrix)

# function to log normalize, find variable features, and scale the data
normVarScaleData<-function(seuratOb,numFeat=2000){
  seuratOb <- NormalizeData(seuratOb)
  seuratOb <- FindVariableFeatures(seuratOb, selection.method = "vst", nfeatures=numFeat)
  seuratOb <- ScaleData(seuratOb, features = rownames(seuratOb))
  return(seuratOb)
}
seurat_object <- NormalizeData(seurat_object) #normalization step

# Add to metadata
seurat_object$manual_clusters <- manual_clusters

# Check metadata
head(seurat_object@meta.data)
table(seurat_object$manual_clusters)  # Verify manual clusters

#Assign Manual Clusters to Idents
seurat_object <- SetIdent(seurat_object, value = manual_clusters)


seurat_object <- FindVariableFeatures(seurat_object, assay = "RNA")
seurat_object <- ScaleData(seurat_object, assay = "RNA")

###### Differential expression between cluster_0 and cluster_1 #####
markers <- FindMarkers(
  object = seurat_object,
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

markers = setorder(markers, -avg_log2FC)

markers[,gene_name:=gene_names[as.numeric(gene)]]

rownames(seurat_object[["RNA"]]) = gene_names

## plot the expression of marker genes in anchor+/- cells 
VlnPlot(seurat_object, features = markers$gene_name[1:10], slot = "data", log = FALSE)
