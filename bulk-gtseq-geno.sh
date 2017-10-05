#!/bin/bash
# usage: bash ./bulk-gtseq-genos.sh locusinfo.txt fasta_files 

LOCUSINFO=$1
FQ="${@:2}"
for file in $FQ; do
  perl /scratch/PI/spalumbi/Conch/Conch_SNP_Test/GTseq_Scripts03072017/14535R_R1/pipeline/GTseq_Genotyper_v2.pl $LOCUSINFO $file > ${file}.genos
done

