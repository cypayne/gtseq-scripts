#!/bin/bash
# usage: bash ./bulk-gtseq-genos.sh path_to_gtseq_scripts locusinfo.txt fasta_files 

GTSEQ_PATH=$1
LOCUSINFO=$2
FQ="${@:3}"
for file in $FQ; do
  perl $GTSEQ_PATH/GTseq_Genotyper_v2.pl $LOCUSINFO $file > ${file}.genos
done

