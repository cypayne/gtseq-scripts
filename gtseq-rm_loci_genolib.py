#!/usr/bin/env python

''' For GTseq: gtseq-rm_loci_genolib.py
    Description: Takes GTseq Library_Genotypes.csv outfile, finds
                 all loci with no data (i.e. '00' -missing data- for
                 all samples) and writes Library_Genotypes.csv
                 contents without these loci to a new file called
                 Library_Genotypes_filtered.csv

    input: Library_Genotypes.csv -
                 a GTseq output file, from which we will find
                 loci with missing data for all samples (and
                 which therefore should be removed)

    output: Library_Genotypes_filtered.csv -
                 Library_Genotypes.csv contents without empty
                 loci
            removed_loci.txt -
                 a list of empty loci that were removed from 
                 Library_Genotypes.csv

    usage: ./gtseq-rm_loci_genolib.py Library_Genotypes.csv 

    CYP 10/04/2017
'''

import sys
from csv import DictReader
from csv import writer

# specify # of columns with no loci data (to skip)
NUM_SKIP = 5
header = True
loci_to_remove = []
loci_to_keep = []

# open the file in universal line ending mode 
with open(sys.argv[1], 'rU') as infile:
  # read the file as a dictionary for each row ({ header : value })
  reader = DictReader(infile)
  fields = reader.fieldnames
  skip_cols = fields[0:NUM_SKIP]
  data = {}
  # create a dictionary of lists for each column/loci name 
  # { locus : [snps] }
  for row in reader:
    for header, value in row.items():
      try:
        data[header].append(value)
      except KeyError:
        data[header] = [value]

  for locus in fields[NUM_SKIP:]:
    if (all(x == '00' for x in data[locus][NUM_SKIP:])) == True:
      # if every SNP value at this locus is missing ('00')
      # then add to list of loci to remove, & remove from dict
      loci_to_remove.append(locus)
      del data[locus]
    else:
      loci_to_keep.append(locus)
  
# write all non-empty loci to LGenos_filtered.csv
new_header = skip_cols + loci_to_keep

# write filtered Library_Genotypes.csv contents to 
# Library_Genotypes_filtered.csv
infile_list = (sys.argv[1]).split('.')
filtered_file = infile_list[0] + '_filtered.' + infile_list[-1]
with open(filtered_file, 'wb') as outfile:
  writer = writer(outfile, delimiter=',')
  writer.writerow(new_header)
  writer.writerows(zip(*[data[key] for key in new_header]))

# dump names of removed loci in removed_loci.txt
with open('removed_loci.txt', 'w') as out:
  for locus in loci_to_remove:
    out.write(locus + '\n')

# Print some info to STDOUT
print('\t************************************************************')
print('\tLGenos_filtered.csv: genotypes library without empty loci. \n \
    \tremoved_loci.txt: contains removed loci names.\n')
print('\t# of loci to start with in Library_Genotypes.csv file: ' + str(len(fields[NUM_SKIP:])))
print('\t# (empty) loci removed: ' + str(len(loci_to_remove)))
print('\t# loci kept: ' + str(len(loci_to_keep)))
print('\t************************************************************')
