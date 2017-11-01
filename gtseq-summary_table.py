#!/usr/bin/env python

''' Takes GTseq seqtest.csv output files and collects
    read count information for each locus per sample.
    
    input: *.seqtest.csv - have format:
              LOCUS,FWD_PRIM,PROBE,OT_RATIO,FWD-PROBE-REV,GHOST_LOCI,NUM_LENGTHS,MOST_FREQ_LENGTH,ALL_LENGTHS:FREQ

    usage: ./gtseq-otread_counts.py *.seqtest.csv
    CYP 10/26/2017
'''
import sys
import csv

info_dict = {} # keep dictionary of locus info
length_dict = {}

first = True # special processing for first file

# read in all modified gtseq seqtest files
for f in sys.argv[1:]:
  header = True
  with open(f, 'r') as seqf:
    sample_name = f.split('.')[0]
    for line in seqf:
      if header: header = False; continue
      items = line.strip().split(',')
      locus = items[0]
      if first: 
        # info_dict should have following fmt:
        # locus: [FWD_PRIM, PROBE, FWD-PROBE-REV, GHOST_LOCI]
        info_dict[locus] = [int(items[1]),int(items[2]),int(items[4]),int(items[5])]
        # length_dict should have following fmt:
        # locus: {length:freq}
        length_dict[locus] = {}
      else:
        # if this is not the first file, then add all relevant 
        # info to existing locus info
          info_dict[locus][0] += int(items[1])
          info_dict[locus][1] += int(items[2])
          info_dict[locus][2] += int(items[4])
          info_dict[locus][3] += int(items[5])

      # add all possible fragment lengths to length_dict
      if int(items[6]) > 0:
        len_freqs = items[8][1:-2].split('/')
        if len(len_freqs) > 0:
          for pair in len_freqs:
            length = pair.split(':')[0]
            freq = pair.split(':')[1]
            if length in length_dict[locus]:
              # if the length is already in the length_dict for this
              # locus, add the frequency count to it
              length_dict[locus][length] += int(freq)
            else:
              # otherwise just add new length
              length_dict[locus][length] = int(freq) 


  first = False 
      
  
with open('./locus_summary_info.csv', 'wt') as out:
  writer = csv.writer(out, delimiter=',')
  to_write = ['LOCUS','FWD_PRIM','PROBE','OT_RATIO','FWD-PROBE-REV','GHOST_LOCI','NUM_LENGTHS','MOST_FREQ_LENGTH','ALL_LENGTHS:FREQ']
  writer.writerow(to_write)
  # write out one locus name per row
  for locus in info_dict: 
    most_freq_length = 0
    if info_dict[locus][0] == 0: 
      otread_ratio = 0
      len_freq_str = '[]'
    else: 
      otread_ratio = float(info_dict[locus][1]) / float(info_dict[locus][0])
      top = True
      len_freq_str = '['
      for l in sorted(length_dict[locus], key=length_dict[locus].get, reverse=True):
        if top: most_freq_length = l; top = False
        len_freq_str += ':'.join([l,str(length_dict[locus][l])]) 
        len_freq_str += '/'
      len_freq_str += ']'
    num_lengths = len(length_dict[locus])
    to_write = [locus,info_dict[locus][0],info_dict[locus][1],otread_ratio,info_dict[locus][2],info_dict[locus][3],num_lengths,most_freq_length,len_freq_str]

    writer.writerow(to_write)
