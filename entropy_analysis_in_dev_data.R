library(DescTools)
library(data.table)
library(ez)

f=fread("/oak/stanford/groups/horence/Roozbeh/Sponge_project/processed_files/sig_RNA_developmental_anchor_counts_granny_mapping.tsv")
f1=f[dataset=="RNA_normal_development"]
f1[,cell_barcode:=NULL]
f1[,cell_type:=NULL]

f1[,entropy:=Entropy(count),by=paste(anchor,Day)] # here we compute the entropy for each anchor as the entropy of the target counts for each anchor and day 
f1[,num_day_per_anchor:=length(unique(Day)),by=anchor]
f2=unique(f1[,list(anchor,Day,entropy)])
anova_results <- ezANOVA(
  data = f2,
  dv = entropy,
  wid = anchor,
  within = Day,
  type = 3  # Type 3 sums of squares
)
ggplot(f2, aes(x = Day, y = entropy, group = anchor, color = anchor)) +
  geom_line() +
  geom_point() +
  labs(title = "Repeated Measures ANOVA: Group Means Over Time", x = "Time", y = "Value") +
  theme_minimal()

setkey(f2,Day,anchor)
diff=f2[Day==12]$entropy-f2[Day==5]$entropy
diff_dt=data.table(f2[Day==12]$anchor,diff)
