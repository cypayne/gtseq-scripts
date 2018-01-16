# gtseq-scripts

In conjunction with Nate Campbell's [GTseq genotyping pipeline](https://github.com/GTseq/GTseq-Pipeline#gtseq-pipeline)
these scripts can be used to efficiently collect individual genotypes from multiplex PCR NGS sequences. This pipeline
is particularly useful if a good reference genome for you study organism does not exist. While these scripts
were specifically written to process one dataset, they can be easily modified to be used with other datasets. 
Please contact the repo owner if modification help is needed.


## How to use

1) Use gtseq-make_infiles.py to automatically generate the GTseq pipeline's required AssayInfo and LocusInfo input files
```
./gtseq-make_infiles.py fastq_file vcf_file primer_file
```
This file produces a ``missing_snps.txt``` file, which contains the names of any loci that did not make it into the input files 
(most likely due to a name mismatch between fastq/vcf/primer infiles)

2) If missing_snps.txt contains loci, you can use the following to manually input those SNPs into your AssayInfo and LocusInfo files:
Use extract_contig_by_name.py to extract a full sequence from your raw fasta/fastq file by its header name
```
./extract_contig_by_name.py raw_data.fastq example_seq_name
```
Then remove non-basepair characters from sequence file (i.e. seq header and newline characters)
```
tail -n+2 example_contig_name.txt | tr -d '\n' > example_seq_name_raw.txt
```
Use extract_seqs.py to extract the probe sequence from your raw sequence file. 
```
./extract_seqs.py example_contig_name_raw.txt index-FLANK index+FLANK+1 
```
For example, if your target SNP is at index 10 and you want your probe to be 7 basepairs long (3 bp flanking each side of your 
SNP), then you would use the following: ./extract_seqs.py example_contig_name_raw.txt 7 14 

3) Run GTseq hashing on all your sample .fastq files using:
```
bash bulk-gtseq-hash.sh *.fastq
```

4) Run GTseq SeqTest on all your .hash files using:
```
bash bulk-gtseq-seqtest.sh assayinfo.txt *.hash
```

5) Run GTseq Genotyper on all your sample .fastq files using:
```
bash bulk-gtseq-geno.sh locusinfo.txt *.fastq
```
