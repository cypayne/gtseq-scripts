#!/bin/bash
# usage: bash ./bulk-gtseq-hash.sh path_to_gtseq_scripts fastq_files_to_hash

GTSEQ_PATH=$1
FQ="${@:2}"
for file in $FQ; do
  perl $GTSEQ_PATH/HashSeqs.pl $file > ${file}.hash 
done
