# gtseq-scripts

These scripts were written to process a specific dataset with the GTSeq genotyping pipeline, however
they can be easily modified to be used with other datasets. Please contact the repo owner if modification
help is needed.

## How to use

1) use gtseq-make_infiles.py to generate the required AssayInfo and LocusInfo input files
```
./gtseq-make_infiles.py fastq_file vcf_file primer_file
```
