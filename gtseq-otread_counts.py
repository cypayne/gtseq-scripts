#!/usr/bin/env python

''' Takes GTseq seqtest.csv output files and collects
    read count information for each locus per sample.
   
    usage: ./gtseq-otread_counts.py *.seqtest.csv
    CYP 10/17/2017
'''

import sys
import csv

# create dictionary of samples, list of sample info per locus
sample_dict = {}
loc_rcounts_dict = {}
samp_rcounts_dict = {}
locus_name_list = ['Sample']
first = True

# read in all gtseq seqtest files
for f in sys.argv[1:]:
  header = True
  samp_otreads = 0
  samp_total_reads = 0
  
  with open(f, 'r') as seqf:
    sample_name = f.split('.')[0]
    locus_ratio_list = []
    for line in seqf:
      if header: header = False; continue

      # collect locus_name, otreads, otratio values
      items = line.strip().split(',')
      locus_name = items[0]
      otreads = int(items[3])
      otratio = items[4]

      # if there are no reads for a locus (i.e.
      # if fwd_prim read count is 0.1), set
      # total_reads to 0
      if float(items[1]) < 1.0:
        total_reads = 0
      else: total_reads = int(items[1])

      if first:
        # collect locus names for outfile header
        locus_name_list.append(locus_name)
        # collect initial otreads and total_reads values per locus
        loc_rcounts_dict[locus_name] = [otreads, total_reads]
      else:
        # add otreads and total_reads counts for each locus
        # in each sample to the running total
        loc_rcounts_dict[locus_name][0] += otreads
        loc_rcounts_dict[locus_name][1] += total_reads 

      # keep track of each locus ot_ratio per sample
      locus_ratio_list.append(otratio)

      # collect running totals for sample total_read and otread
      # counts per sample as well
      samp_total_reads = int(total_reads) + samp_total_reads
      samp_otreads = int(items[3]) + samp_otreads

    # add list of read counts to dictionary with sample as key
    samp_rcounts_dict[sample_name] = [samp_otreads, samp_total_reads]

    # add list of locus ratios to a dictionary, with sample as key
    sample_dict[sample_name] = locus_ratio_list 
    first = False

# write to csv file  
with open('./otread_counts.csv', 'wt') as out:
  writer = csv.writer(out, delimiter=',')

  # write out header of locus names
  new_header = locus_name_list+['#_real_reads','total_#_reads']
  writer.writerow(new_header)

  # write ot_ratios for every locus, per sample
  for sample in sample_dict: 
    to_write = [sample] + sample_dict[sample] + samp_rcounts_dict[sample]
    writer.writerow(to_write)

  # write #_real_reads (i.e. # reads with both primer and probe hits) per locus
  otreads_list = ['#_real_reads']
  for locus in locus_name_list[1:]: otreads_list.append(loc_rcounts_dict[locus][0])
  writer.writerow(otreads_list)

  # write total_#_reads per locus (i.e. all reads with a primer hit for this locus)
  tot_reads_list = ['Total reads']
  for locus in locus_name_list[1:]: tot_reads_list.append(loc_rcounts_dict[locus][1])
  writer.writerow(tot_reads_list)

with open('./locus_read_info.csv', 'wt') as out:
  writer = csv.writer(out, delimiter=',')
  writer.writerow(['Locus','OT Reads','Failed Reads','Total Reads','OT Read Ratio'])
  # write out one locus name per row
  for locus in locus_name_list[1:]:
    failed_reads = loc_rcounts_dict[locus][1] - loc_rcounts_dict[locus][0]
    if loc_rcounts_dict[locus][0] == 0: otread_ratio = 0
    else: otread_ratio = float(loc_rcounts_dict[locus][0]) / float(loc_rcounts_dict[locus][1])
    to_write = [locus, loc_rcounts_dict[locus][0], failed_reads, loc_rcounts_dict[locus][1], otread_ratio]
    writer.writerow(to_write)
