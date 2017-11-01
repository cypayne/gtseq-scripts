#!/usr/bin/env python

''' Removes any sample that has missing data ('00')
    for all loci

    usage: ./gtseq-rm_samples_genolib.py
    
    CYP 10/11/2017
'''

import sys
import csv

NUM_SKIP = 5
samples_to_remove = []
samples_to_keep = []

with open(sys.argv[1], 'rb') as csvfile:
  rdr = csv.reader(csvfile)
  samples = list(rdr)

# want to remove all samples from
for sample in samples:
  if (all(x == '00' for x in sample[NUM_SKIP:-1])) == True:
    samples_to_remove.append(sample[0])
  else: samples_to_keep.append(sample[0])

with open('removed_samples.txt', 'w') as out:
  for sample in samples_to_remove:
    out.write(sample + '\n')

print(samples_to_remove)

