file_path_1003 = "/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/SRP198321/GSM3764782_SRR9051003_latTII1b_raw_gene_bc_matrices_h5.h5"
file_path_1004 = "/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/SRP198321/GSM3764783_SRR9051004_latTII2_raw_gene_bc_matrices_h5.h5"

data_1003 = Read10X_h5(file_path_1003)
seurat_object_1003 <- CreateSeuratObject(counts = data_1003)

data_1004 = Read10X_h5(file_path_1004)
seurat_object_1004 <- CreateSeuratObject(counts = data_1004)

anchor_positive_barcodes = unique(ciona_all[sample=="SRR9051004"]$cell_barcode)
anchor_positive_barcodes = paste(anchor_positive_barcodes,"-1",sep="")
seurat_object_anchor_positive_barcodes_1003 = subset(seurat_object_1003, cells = anchor_positive_barcodes)
seurat_object_anchor_positive_barcodes_1004 = subset(seurat_object_1004, cells = anchor_positive_barcodes)

sum(GetAssayData(object = seurat_object_anchor_positive_barcodes_1003, slot = "counts")["KH2012:KH.L18.59",]>0)


v_1003=as.data.frame(seurat_object_1003@meta.data)
v_1003$gene_name=rownames(v_1003)
v_1003=data.table(v_1003)

v_1004=as.data.frame(seurat_object_1004@meta.data)
v_1004$gene_name=rownames(v_1004)
v_1004=data.table(v_1004)

metadata[,cell_barcode_1:=paste(cell_barcode,"-1",sep=""),by=cell_barcode]

metadata = merge(metadata,v_1003,all.x=T,all.y=F,by.x="cell_barcode_1",by.y="gene_name")
metadata = merge(metadata,v_1004,all.x=T,all.y=F,by.x="cell_barcode_1",by.y="gene_name")
setnames(metadata,c("nCount_RNA.x","nCount_RNA.y"),c("nCount_RNA_1003","nCount_RNA_1004"))
setnames(metadata,c("nFeature_RNA.x","nFeature_RNA.y"),c("nFeature_RNA_1003","nFeature_RNA_1004"))


metadata[,mean(nCount_RNA_1003),by=sample]
metadata[,mean(nCount_RNA_1004),by=sample]