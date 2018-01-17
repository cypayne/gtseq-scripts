#!/bin/bash
# usage: bash ./bulk-gtseq-MOD-seqtest.sh path_to_gtseq_scripts assayinfo.txt hash_files 

GTSEQ_PATH=$1
ASSAYINFO=$2
FQ="${@:3}"
for file in $FQ; do
  BASE=$( basename $file .hash )
  perl $GTSEQ_PATH/GTseq_seqtest_v2_modified.pl $ASSAYINFO $file > ${BASE}.mod_seqtest.csv
done
