#!/usr/bin/env python

''' Script that creates input files (AssayInfo & LocusInfo)
    required by GTSeq pipeline
          
    input:  contig file (fastq) 
            vcf file (vcf)
            primer file (csv):
              Contains the contig #, SNP position, and 
              primer seq (FWD primer only for paired-end 
              seqs), with the format:
                    Contig#-POS,FWD-primer
              where Contig#-POS is the SNP name

    output: assayinfo.txt - first input file
            locusinfo.txt - second input file
            missing_snps.txt - file with snps that were
                               not added to prior 2 output files

    NOTE: make sure snp/contig names match between files 

    usage: ./gtseq-make_infiles.py contigs.fa snp-info.vcf primers.csv

    CYP 09/26/2017
'''
import re
import sys
from Bio import SeqIO
from Bio.Seq import Seq

if len(sys.argv) != 4:
  print("ERROR: Please specify a contigs (fasta), vcf, and " + \
          "primer file (csv).")
  sys.exit()

# global variables
fasta_file = open(sys.argv[1], 'r')
vcf_file = sys.argv[2]
primer_file = sys.argv[3]
vcf_dict = {}
missing_list = []
ref_allele, alt_allele = '', ''
ref_probe, alt_probe = '', ''

# change these variables as needed
FLANK = 3 # num of bps that should flank snp in probe
DEBUG = False # specify True to display DEBUG msgs 

# Process FASTA file
# create dictionary of contigs from fasta file
fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))	
fasta_file.close()

# Process VCF file
# create dictionary of all snps and their ref and 
# alt alleles found in VCF file 
with open(vcf_file, 'r') as vcf:
  for line in vcf:
    if not line.startswith('#'):
      items = line.split('\t') 
      # the full snp name (as it appears in 
      # primer file) is the dict key
      name = items[0]+'-'+items[1]
      # items[3] = ref allele
      # items[4] = alt allele
      vcf_dict[name] = [items[3],items[4]]

# file with contents of vcf_dict 
if DEBUG:
  with open('vcf_dict.txt', 'w') as out:
    for key in vcf_dict:
      out.write(key+': ')
      for value in vcf_dict[key]:
        out.write(value+' ')
      out.write('\n')

assayinfo = open("assayinfo.txt", 'w')
locusinfo = open("locusinfo.csv", 'w')
header = True
with open(primer_file) as prim:
  for line in prim:
    # skip header
    if header: 
      header = False
      continue
    # collect the contig#, SNP position, and primer seq
    # from the primer file (format: Contig#-POS,primer)
    snp_name = line.split(',')[0]
    name_pos = snp_name.split('-') 
    contig_name = name_pos[0] 
    snp_pos = int(name_pos[1])
    if DEBUG: print(snp_pos)
    primer = line.split(',')[1]
#    rev_primer = line.split(',')[2] #XXX for including rev in assayinfo

    #XXX shave off first 20bp (adapter) of each primer
    primer = primer.replace('CGACAGGTTCAGAGTTCTACAGTCCGACGATC','')
#    rev_primer = rev_primer.replace('GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCT','')
#    rev_primer_rc = Seq(rev_primer).reverse_complement()

    # set start and end positions for probe sequence
    probe_start = snp_pos-FLANK-1
    probe_end = snp_pos+FLANK-1
    
    if contig_name in fasta_dict:
        # collect probe sequence from fasta contigs dict
        fasta = fasta_dict[contig_name].seq
        probeseq = fasta[probe_start:probe_end+1] 

    if not snp_name in vcf_dict:
      # if snp name cannot be found in vcf snp dictionary,
      # notify user, add to list of missing snps
      print('%s could not be found in vcf file' % snp_name)
      missing_list.append(snp_name)
    else:
      # otherwise find snp name in vcf dictionary to 
      # collect reference and alternate alleles
      ref_allele = vcf_dict[snp_name][0]
      alt_allele = vcf_dict[snp_name][1]

      # create ref,alt probe sequences using the ref,alt alleles 
      # and the probeseq fragment collect from the fasta file
      ref_probe = str(probeseq[:FLANK]+ref_allele+probeseq[FLANK+1:]) 
      alt_probe = str(probeseq[:FLANK]+alt_allele+probeseq[FLANK+1:]) 

      # snp_dict[snp_name] = [primer,ref_allele,alt_allele,ref_probe,alt_probe]
      # write all relevant values in specific format to each outfile
      assayinfo.write(snp_name+'\t'+primer+'\t'+ref_probe+'\t'+alt_probe+'\n')
#      assayinfo.write(snp_name+'\t'+primer+'\t'+ref_probe+'\t'+alt_probe+'\t'+str(rev_primer_rc)+'\n')
      locusinfo.write(snp_name+','+ref_allele+','+alt_allele+','+ref_probe+','+alt_probe+','+primer+'\n')

# write all snps that weren't added to out files
# to a file called missing_snps.txt
with open('missing_snps.txt', 'w') as out:
  for snp in missing_list:
    out.write(snp+'\n')

assayinfo.close()
locusinfo.close()
