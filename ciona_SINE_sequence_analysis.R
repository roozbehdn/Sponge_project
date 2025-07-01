library(data.table)



riffle <- function(a, b) {   # this function interleaves the elements of two vectors into a vector
  seqmlab <- seq(length=length(a))
  c(rbind(a[seqmlab], b[seqmlab]), a[-seqmlab], b[-seqmlab])
}



ciona_samples = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_SRP198321/samples_Ciona_intestinalis_SRP198321_only_names.txt",header=F)

SINE_seq_count  = data.table()
for (counter in 1:nrow(ciona_samples)){
  sample_counter = ciona_samples$V1[counter]
  counts_sample = fread(paste("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_SRP198321/",sample_counter,"/SINE_seq_satc/satcOut120.tsv",sep=""),header=F)
  names(counts_sample) = c("sample_id","barcode","anchor","target","count")
  counts_sample$sample =   sample_counter
  SINE_seq_count = rbind(SINE_seq_count,counts_sample) # this data table has the counts for SINE sequence across all ciona samples
}
SINE_seq_count[,total_target_count:=sum(count),by=target]
SINE_seq_count$extendor = paste(SINE_seq_count$anchor,SINE_seq_count$target,sep="")
SINE_seq_count[,target_frac:=total_target_count/sum(SINE_seq_count[!duplicated(target)]$total_target_count)]

## making a fasta file of targets for SINE anchor for blast alignment ###########
SINE_seq_count_uniq = SINE_seq_count[!duplicated(target)]
setorder(SINE_seq_count_uniq,-target_frac)
SINE_seq_count_uniq[,num_CATT_repeat:=str_count(extendor,"CATT"),by=target]
targets = SINE_seq_count_uniq[ !grepl("A{5,}", target)]$target
targets_fasta = riffle(paste(">",1:length(targets),sep=""),targets)
targets_fasta = data.table(targets_fasta)
write.table(targets_fasta,"/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/ciona_SINE_targets_after_removing_those_with5A.fa",sep="\t",row.names=F,quote=F)
###################################################


## now below after blasting the targets with less than 5 polyA, i want to analyze its blast results and see how many of them mapped to ciona
## I used this command for blast: /share/software/user/open/ncbi-blast+/2.11.0/bin/blastn -outfmt "6 qseqid stitle sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore sgi sacc slen staxids stitle" -query test.fa  -db nt -out blast_out_test2.txt -task blastn -dust no -taxids 7719 -max_target_seqs 5
blast_output = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/ciona_SINE_targets_after_removing_those_with5A_blast_output.txt")
names(blast_output)=c("qseqid", "stitle", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sgi", "sacc", "slen", "staxids")

blast_output_called = blast_output[evalue<0.05]
targets_dt = data.table(1:length(targets),targets)
names(targets_dt)[2] = "target"
blast_output_called = merge(blast_output_called,targets_dt,all.x=T,all.y=F,by.x="qseqid",by.y="V1") 
blast_output_called = merge(blast_output_called,SINE_seq_count_uniq[,list(target,total_target_count,target_frac,num_samples_per_target,num_CATT_repeat)],all.x=T,all.y=F,by.x="target",by.y="target")

blast_output_called_mapped_to_predicted_genes = blast_output_called[ grepl("PREDICTED", stitle) | grepl("LOC", stitle)]
blast_output_called_mapped_to_predicted_genes[,gene:=strsplit(stitle,split=",")[[1]][1],by=stitle]

## to get the number of unique targets mapped to a gene
length(unique(blast_output_called_mapped_to_predicted_genes$qseqid))
## to get the number of targets mapped by blast
length(unique(blast_output_called$qseqid))
## to get the number of unique genes
length(unique(blast_output_called_mapped_to_predicted_genes$gene))


## below I want to make the boxplot to show the distribution of target fractions for each value
blast_output_called$evalue_categorical = factor(blast_output_called$evalue,levels = sort(unique(blast_output_called$evalue)))
setorder(blast_output_called,evalue)
 ggplot(blast_output_called[!duplicated(target)], aes(evalue_categorical, target_frac)) + geom_boxplot() + theme_bw() +  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +annotation_logticks()
 
 
 
## now below I want to merge with metadata for celltypes
metadata=fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/utility_files/metadata/non_model_organisms/Ciona_intestinalis_SRP198321/ciona_metadata_modified.txt")
metadata$cell_id=paste(metadata$cell,metadata$sample,sep="")
SINE_seq_count$cell=paste(SINE_seq_count$barcode,SINE_seq_count$sample,sep="")
SINE_seq_count = merge(SINE_seq_count,metadata[,list(cell_id,cell_type)],all.x=T,all.y=F,by.x="cell",by.y="cell_id")
