#!/usr/bin/env python

# usage: ./gtseq-make_prim-inter_infile.py primer_file

import sys
from Bio.Seq import Seq

primer_file = sys.argv[1]
out = open('prim-int-input.txt', 'w')

header = True
with open(primer_file) as prim:
  for line in prim:
    # skip header
    if header: 
      header = False
      continue
    # collect the contig#, SNP position, and primer seq
    # from the primer file (format: Contig#-POS,primer)
    items = line.split(',')
    snp_name = items[0]
    fwd_primer = items[1]
    rev_primer = items[2]

    #XXX shave off first 20bp (adapter) of each primer
    fwd_primer = fwd_primer.replace('CGACAGGTTCAGAGTTCTACAGTCCGACGATC','')
    rev_primer = rev_primer.replace('GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT','')

#    rev_primer_rc = Seq(rev_primer).reverse_complement()

#    out.write('\t'.join([snp_name, fwd_primer, str(rev_primer_rc)]) + '\n')
    out.write('\t'.join([snp_name, fwd_primer, rev_primer]) + '\n')

out.close()
