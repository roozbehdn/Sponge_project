library(data.table)

alignment_file = "/oak/stanford/groups/horence/Roozbeh/Sponge_project/data/Ciona_intestinalis/pacbio_wgs/Ciona_SRP198321_cry11bb_compactors_to_gene_modelsAligned.out.sam"
compactors_file ="/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/Ciona_SRP198321_cry11bb_compactors.fasta"
compactors_sample_info_file = "/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/Ciona_SRP198321_cry11bb_compactors.tsv"

alignment = fread(alignment_file,header=FALSE,skip="NH:",select = c("V1","V2","V3","V4","V6","V16"))
setnames(alignment,c("V2","V3","V4","V6","V16"),c("STAR_flag","STAR_chr","STAR_coord","STAR_CIGAR","num_mismatches"))
alignment$V1=paste(">",alignment$V1,sep="")

compactors_fasta = fread(compactors_file,header=F)
compactors_fasta = data.table(compactors_fasta[grepl("seq",V1)]$V1,compactors_fasta[!grepl("seq",V1)]$V1)

compactors_fasta = merge(compactors_fasta, alignment[!duplicated(V1)],all.x=T,all.y=F,by.x="V1",by.y="V1")

compactors_sample_info = fread(compactors_sample_info_file)

compactors_sample_info = merge(compactors_sample_info,compactors_fasta,all.x=T,all.y=F,by.x="compactor",by.y="V2")
compactors_sample_info[!is.na(STAR_CIGAR) & !grepl("S",STAR_CIGAR),num_STAR_soft_clipped:=0]
compactors_sample_info[,num_STAR_soft_clipped:=NULL]
compactors_sample_info[!is.na(STAR_CIGAR),num_STAR_soft_clipped:=max(explodeCigarOpLengths(STAR_CIGAR, ops=c("S"))[[1]]),by=STAR_CIGAR]
compactors_sample_info[is.na(num_STAR_soft_clipped),num_STAR_soft_clipped:=1000]
compactors_sample_info[,number_mismatches:=as.numeric(strsplit(num_mismatches,split=":",fixed=T)[[1]][3]),by=num_mismatches]

compactors_sample_info$num_STAR_soft_clipped=factor(compactors_sample_info$num_STAR_soft_clipped,levels=unique(sort(compactors_sample_info$num_STAR_soft_clipped)))
