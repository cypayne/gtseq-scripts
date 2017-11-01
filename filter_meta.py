#!/usr/bin/env python

''' Creates new meta file that does not
    contain specified sample meta data

    usage: ./filter_meta.py metafile.txt samples_to_remove.txt

    CYP 10/11/2017
'''
import sys

infile = sys.argv[1]
outfile = (infile.split('.'))[0] + '_filtered.txt'
out = open(outfile, 'w')

with open(sys.argv[2], 'r') as rm_file:
  rm_list = [(line.split('.')[0]).strip() for line in rm_file]
print(rm_list)

with open(infile, 'r') as metafile:
  for line in metafile:
    sample = line.split('\t')[0]
    if not sample in rm_list: out.write(line)
out.close()
