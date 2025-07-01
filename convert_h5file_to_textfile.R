library(Seurat)
library(rhdf5)

file_path <- "/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/SRP198321/GSM3764786_larva3_raw_gene_bc_matrices_h5.h5"


barcodes <- h5read(file_path, "/ci_kh_ghostdb_mchsv40_cfpsv40_ensembl_mt/barcodes")
gene_names <- h5read(file_path, "/ci_kh_ghostdb_mchsv40_cfpsv40_ensembl_mt/gene_names")
data <- h5read(file_path, "/ci_kh_ghostdb_mchsv40_cfpsv40_ensembl_mt/data")
indices <- h5read(file_path, "/ci_kh_ghostdb_mchsv40_cfpsv40_ensembl_mt/indices")
indptr <- h5read(file_path, "/ci_kh_ghostdb_mchsv40_cfpsv40_ensembl_mt/indptr")
shape <- h5read(file_path, "/ci_kh_ghostdb_mchsv40_cfpsv40_ensembl_mt/shape")
sparse_matrix <- sparseMatrix(
i = indices + 1,        # Convert 0-based to 1-based indexing
p = indptr,
x = data,
dims = shape
)

# Add row and column names
rownames(sparse_matrix) <- gene_names
colnames(sparse_matrix) <- barcodes

# Convert triplet to data.table
dt <- data.table(
row = triplet$i,      # Row indices
col = triplet$j,      # Column indices
value = triplet$x     # Non-zero values
)
triplet <- summary(sparse_matrix)
# Convert triplet to data.table
dt <- data.table(
row = triplet$i,      # Row indices
col = triplet$j,      # Column indices
value = triplet$x     # Non-zero values
)
dt[, rowname := rownames(sparse_matrix)[row]]
dt[, colname := colnames(sparse_matrix)[col]]

names(dt)[3:5]=c("count","gene_name","cell")
setcolorder(dt,c("gene_name","cell","count","row","col"))
dt[,cell:=gsub("-1","",cell),by=cell]
write.table(dt,"/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/SRP198321/SRR9051007_larva3_rep3_gene_counts.tsv",sep="\t",row.names=F,quote=F)
