#!/bin/bash
# usage: bash ./bulk-gtseq-MOD-seqtest.sh assayinfo.txt hash_files 

ASSAYINFO=$1
FQ="${@:2}"
for file in $FQ; do
  BASE=$( basename $file .hash )
  perl /scratch/PI/spalumbi/Conch/Conch_SNP_Test/GTseq_Scripts03072017/pipeline/GTseq_seqtest_v2_modified.pl $ASSAYINFO $file > ${BASE}.mod_seqtest.csv
done
