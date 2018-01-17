# Modified GTseq Pipeline Protocol

In conjunction with Nate Campbell's [GTseq genotyping pipeline](https://github.com/GTseq/GTseq-Pipeline#gtseq-pipeline)
this protocol can be used to efficiently collect genotype profiles of individual samples from multiplex PCR NGS sequences. 
This pipeline is primarily useful if no good reference genome exists for your study organism. While these scripts
were specifically written to process one dataset, they can be easily modified to be used with other datasets. 
Please contact the repo owner if modification help is needed.


## Before starting this protocol: 
+ Download the GTseq pipeline from GitHub. (Palumbi lab: The pipeline scripts are  available on Sherlock @ `$PI_SCRATCH/gtseq_pipeline`) 
+ Please read through the [GTseq pipeline instructions](https://github.com/GTseq/GTseq-Pipeline/blob/master/GTseq_Pipeline.txt)
to better understand the steps in this protocol.  


## How to use

1) Use `gtseq-make_infiles.py` to automatically generate the GTseq pipeline's required AssayInfo and LocusInfo input files
```
./gtseq-make_infiles.py fastq_file vcf_file primer_file
```
+ Output: `AssayInfo.txt`, `LocusInfo.csv`, `missing_snps.txt`
+ `AssayInfo.txt` is a tab-delimited file with the format: `locus_name  forward_primer_sequence allele_1_probe_sequence allele_2_probe_sequence` 
+ `LocusInfo.csv` is a comma-delimited file with the format: `locus_name,allele_1,allele_2,allele_1_probe_sequence,allele_2_probe_sequence`
+ This file produces a `missing_snps.txt` file, which contains the names of any loci that did not make it into the input files 
(most likely due to a name mismatch between fastq/vcf/primer infiles).  

2) [optional] If `missing_snps.txt` contains loci, you can use the following to help you manually input those SNPs into your AssayInfo and LocusInfo files:
+ Use `extract_contig_by_name.py` to extract a full sequence from your raw fasta/fastq file by its header name
```
./extract_contig_by_name.py raw_data.fastq example_seq_name
```
+ Then remove the non-basepair characters from sequence file (i.e. the seq header and newline characters)
```
tail -n+2 example_seq_name.txt | tr -d '\n' > example_seq_name_raw.txt
```
+ Use `extract_seqs.py` to extract the probe sequence from your raw sequence file. 
```
./extract_seqs.py example_contig_name_raw.txt index-FLANK index+FLANK+1 
```
+ FLANK is the number of basepairs you want to flank your target SNP.
For example, if your target SNP is at index 10 and you want your probe to be 7 basepairs long (3 bp flanking each side of your 
SNP), then you would use the following: `./extract_seqs.py example_contig_name_raw.txt 7 14` 

3) Run GTseq HashSeqs on all your sample .fastq files using:
```
bash bulk-gtseq-hash.sh *.fastq
```
+ Output: one .hash file per .fastq file
+ This script runs each sample through the `GTseq_HashSeqs.pl` script, which collects and counts unique reads in the sequence file.

4) Run GTseq SeqTest on all your .hash files using:
```
bash bulk-gtseq-seqtest.sh AssayInfo.txt *.hash
```
+ Output: one .seqtest.csv file per .hash file
+ This script runs each sample through the `GTseq_SeqTest_v2.pl` script, which counts the occurrence of the forward primer and probe sequences, and the number of times they occur in the same sequence. 

5) Run GTseq Genotyper on all your sample files (.fastq) using:
```
bash bulk-gtseq-geno.sh LocusInfo.txt *.fastq
```
+ Output: one .geno file per .fastq file
+ This script runs each sample through the `GTseq_Genotyper_v2.pl` script, which collects and organizes genotype data and statistics from each fastq .hash file.

6) Run GTseq GenoCompile in the directory containing all of the .geno files generated in step 5:
```
perl GTseq_GenoCompile_v2.pl > Library_Genotypes_allele.csv
```
+ Output: `Library_Genotypes_allele.csv`, which is the desired file containing genotypes
+ This script from the GTseq pipeline compiles genotype data from all .geno files into a single data file.
