#!/bin/bash
# usage: bash ./bulk-gtseq-seqtest.sh path_to_gtseq_scripts assayinfo.txt hash_files 

GTSEQ_PATH=$1
ASSAYINFO=$2
FQ="${@:3}"
for file in $FQ; do
  BASE=$( basename $file .hash )
  perl $GTSEQ_PATH/GTseq_SeqTest_v2.pl $ASSAYINFO $file > ${BASE}.seqtest.csv
done
