#!/usr/bin/env python

''' For GTseq: gtseq-gtseq-genolib_cleanup.py
    Description: Takes GTseq Library_Genotypes.csv outfile, finds
                 all loci with no data (i.e. '00' -missing data- for
                 all samples) and writes Library_Genotypes.csv
                 contents without these loci to a new file called
                 Library_Genotypes_filtered.csv

    input: Library_Genotypes.csv -
                 a GTseq output file, from which we will find
                 loci with missing data for all samples (and
                 which therefore should be removed)
           %_threshold - 
                 integer value, corresponds to the percentage of
                 samples a locus must have data for in order for
                 a locus to be kept
                 ex. if %_threshold is 90, then you want to keep
                     any locus that has data for 90% or more of the
                     samples (i.e. only 10% missing data, or less, 
                     is allowed for a locus to be kept)

    output: Library_Genotypes_filtered.csv -
                 Library_Genotypes.csv contents without empty
                 loci
            removed_loci.txt -
                 a list of loci that were removed from 
                 Library_Genotypes.csv, and that should be
                 removed from input files
            removed_samples.txt -
                a list of samples that were removed from
                Library_Genotypes.csv, and that should be 
                removed from analysis pool

    usage: ./gtseq-rm_loci_genolib.py Library_Genotypes.csv %_threshold_missing_data 

    CYP 10/04/2017
'''


import sys
import csv 

# specify # of columns with no loci data (to skip)
NUM_SKIP = 5 

genolib_file = sys.argv[1]

''' Returns list of samples with 100% missing data '''
def get_empty_samples():
  samples_to_remove = []
  samples_to_keep = []

  with open(genolib_file, 'rb') as csvfile:
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

  # note: samples list includes header (so subtract 1 from length)
  print('\tinitial # of samples: ' + str(len(samples)-1))
  print'\t# samples kept: ' + str(len(samples_to_keep))
  print('\t# samples removed: ' + str(len(samples_to_remove)))

  return samples_to_remove


''' Removes loci that have less non-missing data than the 
    specified % threshold, and simultaneously removes 
    empty samples (returned in get_empty_samples)
    Parameter: samples_to_remove - list of samples to remove
'''
# pass in list of samples to remove
def remove_empty_loci(samples_to_remove):

  # percentage of data that is missing to qualify
  # getting rid of a loci
  threshold = (100.0 - int(sys.argv[2]))

  header = True
  loci_to_remove = []
  loci_to_keep = []

  row_count = 0.0
  with open(genolib_file, 'rU') as infile:
    # read the file as a dictionary for each row ({ header : value })
    reader = csv.DictReader(infile)
    fields = reader.fieldnames
    # skip columns that don't contain genotypes
    skip_cols = fields[0:NUM_SKIP]
    data = {}

    # create a dictionary of lists for each column/loci name 
    # { locus : [snps] }
    for row in reader:
      row_count += 1
      # skip rows with samples that should be removed
      if row['Sample'] in samples_to_remove: continue
      for header, value in row.items():
        try:
          data[header].append(value)
        except KeyError:
          data[header] = [value]

    for locus in fields[NUM_SKIP:]:
      remove = False
      missing_data_ctr = 0.0
      for x in data[locus][NUM_SKIP:]:
        # if every SNP value at this locus is missing ('00')
        # then add to list of loci to remove, & remove from dict
        if x == '00': missing_data_ctr += 1
        missing_ratio = (missing_data_ctr/row_count)*100.0 
        if missing_ratio >= threshold: 
          remove = True
          break
      if remove:
        loci_to_remove.append(locus)
        del data[locus]
      else:
        loci_to_keep.append(locus)
    
  # write all non-empty loci to LGenos_filtered.csv
  new_header = skip_cols + loci_to_keep

  # write filtered Library_Genotypes.csv contents to 
  # Library_Genotypes_filtered.csv
  infile_list = genolib_file.split('.')
  filtered_file = infile_list[0] + '_filtered.' + infile_list[-1]
  with open(filtered_file, 'wb') as outfile:
    writer = csv.writer(outfile, delimiter=',')
    writer.writerow(new_header)
    writer.writerows(zip(*[data[key] for key in new_header]))

  # dump names of removed loci in removed_loci.txt
  with open('removed_loci.txt', 'w') as out:
    for locus in loci_to_remove:
      out.write(locus + '\n')

  # Print some info to STDOUT
  print('\tinitial # of loci: ' + str(len(fields[NUM_SKIP:])))
  print('\t# loci kept: ' + str(len(loci_to_keep)))
  print('\t# (empty) loci removed: ' + str(len(loci_to_remove)))


### MAIN ###

print('\t************************************************************')
print('\tLGenos_filtered.csv: genotypes library without empty loci. \n \
    \tremoved_loci.txt: contains removed loci names.\n \
    \tremoved_samples.txt: contains removed sample names.\n')

samples_to_remove = get_empty_samples()
remove_empty_loci(samples_to_remove)

print('\t************************************************************')
