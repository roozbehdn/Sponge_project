# Here I want to align unaligned reads for lemur to human and then show that they may align to MHC
library(data.table)

riffle <- function(a, b) {   # this function interleaves the elements of two vectors into a vector
  seqmlab <- seq(length=length(a))
  c(rbind(a[seqmlab], b[seqmlab]), a[-seqmlab], b[-seqmlab])
}

anchors=fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/Lemur_10x_classification.tsv")
anchors=anchors[is.aligned_STAR==0]

STAR_index_path = "/oak/stanford/groups/horence/Roozbeh/T2T_human_genome/STAR_index_files"
Bowtie_index_path = "/oak/stanford/groups/horence/Roozbeh/T2T_human_genome/Bowtie_index_files/chm13v2_bt2"
known_splice_sites_file = "/oak/stanford/groups/horence/Roozbeh/T2T_human_genome/T2T_chm13.draft_v2.0_known_splice_sites.txt"
known_exon_boundaries_file = "/oak/stanford/groups/horence/Roozbeh/T2T_human_genome/T2T_CAT_CHM13v2_exon_coordinates.bed"
gene_coordinates_file = "/oak/stanford/groups/horence/Roozbeh/T2T_human_genome/T2T_genes_cleaned.bed"
paralogs_file = "/oak/stanford/groups/horence/Roozbeh/T2T_human_genome/human_paralogs.tsv"


## now align unaligned extendors to human
anchors[,extendor_index:=paste(donor,tissue,extendor_index,sep="_")]
anchors$extendor_index = paste(">",anchors$extendor_index,sep="")
extendors_fasta = riffle(anchors$extendor_index,anchors$extendor)
extendors_fasta  = data.table(extendors_fasta) # the fasta file containing all extendor sequences

write.table(extendors_fasta,"/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/Lemur_10x_unaligned_extendors.fa", quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t")
directory = "/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/"
system(paste("/oak/stanford/groups/horence/Roozbeh/software/STAR-2.7.5a/bin/Linux_x86_64/STAR --runThreadN 4 --genomeDir  ",STAR_index_path, " --readFilesIn ", directory,"Lemur_10x_unaligned_extendors.fa", " --outFileNamePrefix ", directory,"STAR_alignment/Lemur_10x_unaligned_extendors_to_human"," --twopassMode Basic --alignIntronMax 1000000 --limitOutSJcollapsed 3000000 --chimJunctionOverhangMin 10 --chimSegmentReadGapMax 0 --chimOutJunctionFormat 1 --chimSegmentMin 12 --chimScoreJunctionNonGTAG -4 --chimNonchimScoreDropMin 10 --outSAMtype SAM --chimOutType SeparateSAMold --outSAMunmapped None --clip3pAdapterSeq AAAAAAAAA --outSAMattributes NH HI AS nM NM ",sep = ""))


alignment_info_extendors = fread(paste(directory,"STAR_alignment/Lemur_10x_unaligned_extendors_to_humanAligned.out.sam",sep=""),header=FALSE,skip="NH:",select = c("V1","V2","V3","V4","V6","V16")) ## now grabbing alignment information for extendor sequences after running STAR
alignment_info_extendors$V1=paste(">",alignment_info_extendors$V1,sep="")
alignment_info_extendors[,num_human_alignments:=.N,by=V1]
alignment_info_extendors[,STAR_human_num_mismatches:=as.numeric(strsplit(V16,split=":")[[1]][3]),by=V16]

anchors = merge(anchors,alignment_info_extendors[!duplicated(V1),list(V1,V2,V3,V4,V6,num_alignments,STAR_num_mismatches)],all.x=TRUE,all.y=FALSE,by.x="extendor_index",by.y="V1")
setnames(anchors,c("V2","V3","V4","V6","num_alignments"),c("STAR_human_flag","STAR_human_chr","STAR_human_coord","STAR_human_CIGAR","STAR_human_num_alignments"))


### below I add the flags for the STAR alignment status of each extendor
anchors[,extendor_index:=gsub(">","",extendor_index),by=extendor_index]
anchors[,is.human_aligned_STAR:=0]   # whether extendor has been mapped by STAR
anchors[,is.human_STAR_SJ:=0]        # whether STAR reports a gapped alignment (splice junction for the extendor)

anchors[STAR_human_CIGAR%like%"N", is.human_STAR_SJ:=1]    
anchors[!is.na(STAR_human_chr), is.human_aligned_STAR:=1]  # the extendors with non-NA STAR alignment are flagged as mapped by STAR

system(paste("/home/groups/horence/applications/samtools-0.1.19/samtools view -S -b ",directory,"STAR_alignment/extendorsAligned.out.sam > ",directory,"STAR_alignment/extendorsAligned.out.bam",sep=""))
system(paste("/home/groups/horence/applications/bedtools2/bin/bedtools bamtobed -split -i ",directory,"STAR_alignment/extendorsAligned.out.bam | sort -k1,1 -k2,2n > ", directory,"STAR_alignment/called_exons.bed",sep=""))
system(paste("/home/groups/horence/applications/bedtools2/bin/bedtools intersect -a ",directory,"STAR_alignment/called_exons.bed -b ", gene_coordinates_file, " -wb -loj | cut -f 4,10   | /home/groups/horence/applications/bedtools2/bin/bedtools groupby -g 1 -c 2 -o distinct  > ", directory,"STAR_alignment/extendor_genes.txt",sep=""))

extendor_genes = fread(paste(directory,"STAR_alignment/extendor_genes.txt",sep=""),sep="\t",header=FALSE)
names(extendor_genes) = c("extendor_index","extendor_human_gene")

anchors = merge(anchors,extendor_genes[!duplicated(extendor_index)],all.x=TRUE,all.y=FALSE,by.x="extendor_index",by.y="extendor_index")
anchors[extendor_human_gene==".",extendor_human_gene:=NA]
