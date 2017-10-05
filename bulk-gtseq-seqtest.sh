#!/bin/bash
# usage: bash ./bulk-gtseq-seqtest.sh assayinfo.txt hash_files 

ASSAYINFO=$1
FQ="${@:2}"
for file in $FQ; do
  BASE=$( basename $file .hash )
  perl /scratch/PI/spalumbi/Conch/Conch_SNP_Test/GTseq_Scripts03072017/14535R_R1/pipeline/GTseq_SeqTest_v2.pl $ASSAYINFO $file > ${BASE}.seqtest.csv
done
