#!/bin/bash
# usage: bash ./bulk-gtseq-hash.sh fastq_files_to_hash

for file in $@; do
  perl /scratch/PI/spalumbi/Conch/Conch_SNP_Test/GTseq_Scripts03072017/14535R_R1/pipeline/HashSeqs.pl $file > ${file}.hash 
done
