#!/bin/sh
#############################
# File Name : make_separate_fasta_files.sh
#
# Purpose : [???]
#
# Creation Date : 19-02-2024
#
# Last Modified : Mon 19 Feb 2024 08:58:28 PM PST
#
# Created By : Roozbeh Dehghannasiri
#
##############################

mkdir new
cd new
csplit -s -z /path/to/INPUT.FA '/>/' '{*}'
for i in xx* ; do \
  n=$(sed 's/>// ; s/ .*// ; 1q' "$i") ; \
  mv "$i" "$n.fa" ; \
done
