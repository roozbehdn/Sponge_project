library(data.table)



riffle <- function(a, b) {   # this function interleaves the elements of two vectors into a vector
  seqmlab <- seq(length=length(a))
  c(rbind(a[seqmlab], b[seqmlab]), a[-seqmlab], b[-seqmlab])
}



ciona_samples = fread("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_SRP198321/samples_Ciona_intestinalis_SRP198321_only_names.txt",header=F)

leader_seq_count  = data.table()
for (counter in 1:nrow(ciona_samples)){
  sample_counter = ciona_samples$V1[counter]
  counts_sample = fread(paste("/oak/stanford/groups/horence/Roozbeh/NOMAD_10x/runs/marine_organisms/Ciona_intestinalis/Ciona_intestinalis_SRP198321/",sample_counter,"/leader_seq_satc/satcOut111.tsv",sep=""),header=F)
  names(counts_sample) = c("sample_id","barcode","anchor","target","count")
  counts_sample$sample =   sample_counter
  leader_seq_count = rbind(leader_seq_count,counts_sample) # this data table has the counts for leader sequence across all ciona samples
}
leader_seq_count[,total_target_count:=sum(count),by=target]
leader_seq_count$extendor = paste(leader_seq_count$anchor,leader_seq_count$target,sep="")
leader_seq_count_uniq_target = leader_seq_count[!duplicated(target)]

## removing those targets with a polyA tail of loner than 4 (have at least 5As)
compute_polya_length <- function(sequence) {
  # Reverse the sequence to count 'A's from the end
  reversed_seq <- rev(strsplit(sequence, NULL)[[1]])
  
  # Count the length of consecutive 'A's from the end
  tail_length <- 0
  for (base in reversed_seq) {
    if (base == "A") {
      tail_length <- tail_length + 1
    } else {
      break
    }
  }
  return(tail_length)
}
leader_seq_count_uniq_target[, polya_length := sapply(target, compute_polya_length)]

leader_seq_count_uniq_target = leader_seq_count_uniq_target[polya_length<5]
extendors_fasta = riffle(paste(">",1:nrow(leader_seq_count_uniq_target),sep=""),leader_seq_count_uniq_target$target)
extendors_fasta = data.table(extendors_fasta)

write.table(extendors_fasta,"/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/ciona_leader_targets_after_removing_those_with5A.fa",sep="\t",row.names=F,quote=F)


blast_output = fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/ciona_leader_targets_after_removing_those_with5A_blast_output.txt")
names(blast_output)=c("qseqid", "stitle", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "sgi", "sacc", "slen", "staxids")

blast_output_called = blast_output[evalue<0.05]

targets_dt = data.table(1:length(leader_seq_count_uniq_target$target),leader_seq_count_uniq_target$target)
names(targets_dt)[2] = "target"
blast_output_called = merge(blast_output_called,targets_dt,all.x=T,all.y=F,by.x="qseqid",by.y="V1") 
blast_output_called = merge(blast_output_called,leader_seq_count_uniq[,list(target,total_target_count)],all.x=T,all.y=F,by.x="target",by.y="target")