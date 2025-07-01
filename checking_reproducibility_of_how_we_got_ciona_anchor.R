library(data.table)

f=fread("/oak/stanford/groups/horence/george/protein_domain_project/workflow/ent_geq3_amazon_plus_splash_stats_repeats_nonoverlapping_pfam_hits.tsv")

## Julia's original command as mentioned in her email to Ayelet on how she found the ciona anchor:
f[order(-target_entropy)][periodic_repeat_size>20][generic_repeat_3mer_entropy>1][generic_repeat_4mer_entropy>1][generic_repeat_5mer_entropy>1][effect_size_bin>.4][target_entropy>3][sequencing_technology %like% "10x"][1]

## I can essentially remove some criteria and still get the same anchor
f[order(-target_entropy)][periodic_repeat_size>20][generic_repeat_5mer_entropy>1][effect_size_bin>.4][target_entropy>3][sequencing_technology %like% "10x"][1]

## but I can also assume that even without [generic_repeat_5mer_entropy>1] I can still get the ciona anchor as the top hit because the top anchor was a simple GT repeat which we can assume can be filteredout by artifact removal
f[order(-target_entropy)][periodic_repeat_size>20][effect_size_bin>.4][target_entropy>3][sequencing_technology %like% "10x"][1]

f[order(-target_entropy)][periodic_repeat_size>20][generic_repeat_3mer_entropy>1][generic_repeat_4mer_entropy>1][generic_repeat_5mer_entropy>1][effect_size_bin>.4][target_entropy>3][sequencing_technology %like% "10x"][!duplicated(anchor)][,list(anchor,target_entropy,Experiment,Organism)]
