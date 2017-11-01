#!/bin/bash
# usage: bash ./bulk-gtseq-priminter.sh Primer_Interaction_input_file hash_files 

INFILE=$1
HASH="${@:2}"
for file in $HASH; do
  BASE=$( basename $file .hash )
  perl /scratch/PI/spalumbi/Conch/Conch_SNP_Test/GTseq_Scripts03072017/pipeline/GTseq_Primer-Interaction-Test_v2.pl $INFILE $file > ${BASE}.primint.txt
done
