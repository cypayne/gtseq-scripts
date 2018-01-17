#!/bin/bash
# usage: bash ./bulk-gtseq-priminter.sh path_to_gtseq_scripts Primer_Interaction_input_file hash_files 

GTSEQ_PATH=$1
INFILE=$2
HASH="${@:3}"
for file in $HASH; do
  BASE=$( basename $file .hash )
  perl $GTSEQ_PATH/GTseq_Primer-Interaction-Test_v2.pl $INFILE $file > ${BASE}.primint.txt
done
