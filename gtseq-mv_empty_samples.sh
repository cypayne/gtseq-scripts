#!/bin/bash

# First run one of the following to get a 'remove_samples.txt' file:
#   gtseq-rm_samples_genolib.py 
#   gtseq-genolib_cleanup.py                               
# usage: bash ./gtseq-mv_empty_samples.sh remove_samples.txt
# CYP 10/11/2017

# make empty_samples dir if there isn't one
if [ ! -d "./empty_samples" ]; then
  mkdir "./empty_samples"
fi

for file in $(<$1); do
  mv "$file" "./empty_samples/"
  mv "${file}.hash" "./empty_samples/"
done
