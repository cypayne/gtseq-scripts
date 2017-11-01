#!/usr/bin/env python

'''
    Generates a meta file for multiplex data,
    to be used as input for GT-Seq-SNP.R
    usage: ./generate_metafile.py multiplex_info_file
    output: multiplex_info_file_meta.txt
    CYP 10/09/2017
'''

import re 
import sys

locs = ['Aruba', 'Belize', 'Bahamas', 'Florida', 'Jamaica', 'Statia']
def get_site(x, locs):
  return {
      'AR': locs[0],
      'BZ': locs[1],
      'EX': locs[2],
      'FK': locs[3],
      'PB': locs[4],
      'SE': locs[5] 
  }.get(x, 'Other')

infile = sys.argv[1]
outfile = (infile.split('.'))[0] + '_meta.txt'
out = open(outfile, 'w')

first = True
with open(infile, 'r') as inf:
  for line in inf:
    if first: 
      first = False
      continue
    items = line.split('\t')
    id_name = items[1]
    samp_name = items[2]
    
    location = re.split('(\d+)',samp_name)
    site = get_site(location[0], locs)

    # if the first round of site determination
    # isn't specific enough, check for the name 
    # of the actual site in the sample name
    if site == 'Other':
      for s in locs:
        temp = re.match(s, samp_name)
        if temp: site = temp.group(0)

    out.write(id_name + '\t' + site + '\n')

out.close()
