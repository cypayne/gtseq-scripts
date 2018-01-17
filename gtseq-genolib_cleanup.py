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
                 integer value, corresponds to the maximum percentage of 
                 missing genotypes a locus can have in order to be kept
                 ex. if %_threshold is 10, then you want to keep
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

    usage: ./gtseq-rm_loci_genolib.py Library_Genotypes.csv %_threshold_loci 
            %_threshold_sample 

    CYP 10/04/2017
'''

import sys
import csv 

# specify # of columns with no loci data (to skip)
NUM_SKIP = 5 

genolib_file = sys.argv[1]

''' Returns list of samples with 100% missing data 
    NEW - also include samples with 90% or more missing data
def get_empty_samples():
  samples_to_remove = []
  samples_to_keep = []
  threshold = 50 

  with open(genolib_file, 'rb') as csvfile:
    rdr = csv.reader(csvfile)
    samples = list(rdr)

  first = True 
  # want to remove all samples from
  for sample in samples:
    if first: first=False; continue
#    if (all(x == '00' for x in sample[NUM_SKIP:-1])) == True:
    # XXX new addition, remove any samp with 90%+ missing data
    missing_ratio = 0.0
    missing_data_ctr = 0.0
    remove = False
    for x in sample[NUM_SKIP:-1]:
      if x == '00': missing_data_ctr += 1
      missing_ratio = (missing_data_ctr*100.0)/len(sample[NUM_SKIP:-1]) 
      if missing_ratio >= threshold: 
        remove = True
        break
    if remove:
      samples_to_remove.append(sample[0])
    else:
      samples_to_keep.append(sample[0])

  with open('removed_samples.txt', 'w') as out:
    for sample in samples_to_remove:
      out.write(sample + '\n')

  # note: samples list includes header (so subtract 1 from length)
  print('\tinitial # of samples: ' + str(len(samples)-1))
  print'\t# samples kept: ' + str(len(samples_to_keep))
  print('\t# samples removed: ' + str(len(samples_to_remove)))

  return samples_to_remove
'''

''' Input:
              empty_locus_index_list: list of indices corresponding
                                      to empty loci
              nonempty_locus_count:   integer value of number of loci
                                      that contain at least one genotype
    Returns:
              rm_sample_indices: list of sample indices to be ignored
                                 when evaluating loci
              samples_to_remove: list of sample names to be removed
'''
def samps_to_remove(empty_locus_index_list, nonempty_locus_count):  
#  print(empty_locus_index_list)
#  print(nonempty_locus_count)
  rm_sample_indices = []
  samples_to_remove = []
  samples_to_keep = []

  threshold = float(sys.argv[3]) 

  # collect sample data from Lib_Geno input file
  with open(genolib_file, 'rb') as csvfile:
    rdr = csv.reader(csvfile)
    samples = list(rdr)

  first = True 
  for samp_index,sample in enumerate(samples):
    # ignore first line (header)
    if first: first=False; continue

    # remove any samp with threshold%+ missing data
    missing_ratio = 0.0
    missing_data_ctr = 0.0
    remove = False
    for i,x in enumerate(sample[NUM_SKIP:-1]):
      # if the current index matches one of the 
      # removed locus indices, skip to simulate removal
      if i in empty_locus_index_list: continue 

      # otherwise include locus genotype
      if x == '00': 
        missing_data_ctr += 1
        # re-calculate percent of missing locus data in this sample
        missing_ratio = (missing_data_ctr*100.0)/nonempty_locus_count
        if missing_ratio >= threshold: remove = True; break

    if remove:
      samples_to_remove.append(sample[0])
      rm_sample_indices.append(samp_index)
    else:
      samples_to_keep.append(sample[0])

  # write removed sample names to output file
  with open('removed_samples.txt', 'w') as out:
    for sample in samples_to_remove:
      out.write(sample + '\n')

  # note: samples list includes header (so subtract 1 from length)
  print('\tinitial # of samples: ' + str(len(samples)-1))
  print'\t# samples kept: ' + str(len(samples_to_keep))
  print('\t# samples removed: ' + str(len(samples_to_remove)))

  return rm_sample_indices, samples_to_remove


''' Removes loci that have less non-missing data than the 
    specified % threshold, and simultaneously removes 
    empty samples (returned in get_empty_samples)
'''
def remove_empty_loci():

  # percentage of data that is missing to qualify
  # getting rid of a locus
  threshold = float(sys.argv[2])

  rm_loci_indices = []
  loci_to_remove = []
  loci_to_keep = []

  row_count = 0.0
  with open(genolib_file, 'rU') as infile:
    # read the file as a dictionary for each row ({ header : value })
    reader = csv.DictReader(infile)
    fields = reader.fieldnames
    # skip columns that don't contain locus genotypes
    skip_cols = fields[0:NUM_SKIP]
    data = {}

    # create a dictionary of lists for each column/loci name 
    # { locus : [snps] }
    for row in reader:
      row_count += 1
      for header, value in row.items():
        try:
          data[header].append(value)
        except KeyError:
          data[header] = [value]

    # first record and remove every locus that is completely empty
    for index,locus in enumerate(fields[NUM_SKIP:-1]):
      if (all(x == '00' for x in data[locus])) == True:
        loci_to_remove.append(locus)
        rm_loci_indices.append(fields[NUM_SKIP:-1].index(locus))
        del data[locus]
    nonempty_locus_count = (len(data)-NUM_SKIP-1)

    # get list of samples, list of indices to remove 
    rm_sample_indices, samples_to_remove = samps_to_remove(rm_loci_indices, nonempty_locus_count)

    # remove loci (non-empty) with more than threshold of missing data
    for locus in fields[NUM_SKIP:-1]:
      if locus in loci_to_remove: continue
      remove = False
      missing_data_ctr = 0.0
      for index,x in enumerate(data[locus]):
        # simulate sample removal by skipping sample indices
        if index in rm_sample_indices: continue 
        # if every SNP value at this locus is missing ('00')
        # then add to list of loci to remove, & remove from dict
        if x == '00': 
          missing_data_ctr += 1
          missing_ratio = (missing_data_ctr*100.0)/row_count 
          if missing_ratio >= threshold: remove = True; break
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

    # transpose locus dictionary to rewrite proper format in
    # output file
    samp_list = (zip(*[data[key] for key in new_header]))
    for samp in samp_list:
      # if a sample name appears in rm list, don't write
      if samp[0] in samples_to_remove: continue
      else: writer.writerow(samp)

  # dump names of removed loci in removed_loci.txt
  with open('removed_loci.txt', 'w') as out:
    for locus in loci_to_remove:
      out.write(locus + '\n')

  # Print some info to STDOUT
  print('\tinitial # of loci: ' + str(len(fields[NUM_SKIP:-1])))
  print('\t# loci kept: ' + str(len(loci_to_keep)))
  print('\t# (empty) loci removed: ' + str(len(loci_to_remove)))


### MAIN ###

print('\t************************************************************')
print('\tLGenos_filtered.csv: genotypes library without empty loci. \n \
    \tremoved_loci.txt: contains removed loci names.\n \
    \tremoved_samples.txt: contains removed sample names.\n')

remove_empty_loci()

print('\t************************************************************')
