#!/usr/bin/perl
#Test .hash file to see how many times loci fwd primer plus probe seqs occur

# pass in tab-delimited assayinfo file with RC rev primer appended to end
# usage: perl ./GTseq_seqtest_v2_modified.pl assyinfo_with_rev-primer.txt
#
#
# modifications added by CYP 10/26/2017

use strict; use warnings;

die "usage: provide <tab delimited txt file of locus, fwd, probe1, probe2 seqs, rev_RC> and <*.hash>\n" unless @ARGV == 2;

my %fwd_seq = ();
my %fwd_seqAssay = ();
my %probe1 = ();
my %probe2 = ();
my %probe1RC = ();
my %probe2RC = ();
my %rev_seq = (); #XXX
#TODO my %rev_seqAssay = (); #XXX
my %fwd_count = ();
my %probe_count = ();
my %both_count = ();
my %OT_Ratio = ();

# XXX
my %both_fwd_rev_count = (); 
my %lengths = ();
my %ghost_locus_count = ();
my %additional_lengths = ();
my %with_ghosts = ();

#read in assay and allele information and push to arrays...
open(SEQ, "<$ARGV[0]") or die "error reading $ARGV[0]\n";

while (<SEQ>) {
	chomp;
	my @info = split(/\t/, $_);
	my $tring = substr $info[1], 0, 14;
  # storing locus name here
	$fwd_seq{$info[0]} = $tring;
	$fwd_seqAssay{$tring} = $info[0];
	$probe1{$info[0]} = $info[2];
	$probe2{$info[0]} = $info[3];
	my $p1RC = reverse $info[2];
	my $p2RC = reverse $info[3];
	$p1RC =~ tr/ACGT][/TGCA[]/;
	$p2RC =~ tr/ACGT][/TGCA[]/;
	$probe1RC{$info[0]} = $p1RC;
	$probe2RC{$info[0]} = $p2RC;

  #XXX new addition
  my $tring_rev = substr $info[4], 0, 14;
  $rev_seq{$info[0]} = $tring_rev;
#TODO  $rev_seqAssay{$tring_rev} = $info[4];
}
close SEQ;

#initialize stats at zero...
foreach my $assays (sort keys %fwd_seq){
$fwd_count{$assays} = 0;
$probe_count{$assays} = 0;
$both_count{$assays} = 0;
$OT_Ratio{$assays} = 0;
# XXX
$both_fwd_rev_count{$assays} = 0;
$ghost_locus_count{$assays} = 0;
$lengths{$assays} = 0;
$additional_lengths{$assays} = {};
$with_ghosts{$assays} = [];
}
open(HASH, "<$ARGV[1]") or die "error reading $ARGV[1]\n";
while (<HASH>) {
  my $hash = $_;
  my $R1_seq = <HASH>;
  my $seq_start = substr $R1_seq, 0, 14;
  chomp ($hash);
  my @info = split(/;/, $hash);
  my $count = $info[2];
  # XXX
  my $probe_found = '0'; 
  my $fragment = '';
  my $len = 0;

  if (exists $fwd_seqAssay{$seq_start}) {
    # Increment a locus' fwd primer count by the amount 
    # specified in hash seq header if hash seq begins
    # with that locus 
    # $assay_hit = locus name
    my $assay_hit = $fwd_seqAssay{$seq_start};
    $count = $fwd_count{$assay_hit} + $count;
    $fwd_count{$assay_hit} = $count;
    $count = $info[2];

    if ($R1_seq =~ m/$probe1{$assay_hit}|$probe2{$assay_hit}|$probe1RC{$assay_hit}|$probe2RC{$assay_hit}/){
      # if alternate or reference probe or RC of alt or ref probe is found in seq 
      # after fwd primer, increment probe count by amount specified in hash seq header
      $count = $probe_count{$assay_hit} + $count;
      $probe_count{$assay_hit} = $count;
      $probe_found = '1'; #XXX
    }
    $count = $info[2];
    if (($R1_seq =~ m/^$fwd_seq{$assay_hit}/) && ($R1_seq =~ m/$probe1{$assay_hit}|$probe2{$assay_hit}|$probe1RC{$assay_hit}|$probe2RC{$assay_hit}/)) {
      # if both fwd primer and probe occur in same seq, increment
      # both count by amount specified in hash seq header
      $count = $both_count{$assay_hit} + $count;
      $both_count{$assay_hit} = $count;
    }

    # XXX
    $count = $info[2];
    if ($R1_seq =~ m/$rev_seq{$assay_hit}/){
      if ($probe_found) {
        # if rev primer found in seq, increment both_fwd_rev_count
        # by amount specified in hash seq header
        $count = $both_fwd_rev_count{$assay_hit} + $count;
        $both_fwd_rev_count{$assay_hit} = $count;
      }
      else {
        $count = $ghost_locus_count{$assay_hit} + $count;
        $ghost_locus_count{$assay_hit} = $count;
        # Note: to print this info, uncomment "ghosts" block below
        push(@{$with_ghosts{$assay_hit}},$R1_seq);
      } 
      # get length of fragment between fwd and rev, hash for locus lengths
      ($fragment) = $R1_seq =~ /$fwd_seq{$assay_hit}(.*?)$rev_seq{$assay_hit}/;
      $len = length $fragment;
#      print "$assay_hit\t";
#      print "$len\n";
#      print "$fragment;\n"; 
#      print "$R1_seq;~\n";
#      if (($lengths{$assay_hit} ne "0") && ($lengths{$assay_hit} ne $len) && !(grep {$_ eq $len} @{$additional_lengths{$assay_hit}})) { 
      if (exists $additional_lengths{$assay_hit}{$len}) { 
        # if this length was already recorded for this locus, increment its count
#       push(@{$additional_lengths{$assay_hit}},$len);
        $additional_lengths{$assay_hit}{$len}++;
#       print "$assay_hit\t$len\n";
      }
      else {
        # otherwise make new instance of length in lengths hash table
#        $lengths{$assay_hit} = $len;
        $additional_lengths{$assay_hit}{$len} = 1; 
      }
    }
  }
}
close HASH;

### GHOSTS ###

#print "************ghosts******************";
#foreach my $locus (sort keys %fwd_seq) {
#  if (@{$with_ghosts{$locus}} ne []) { print "$locus\t@{$with_ghosts{$locus}};\n";}
#}

### LENGTHS ###
#for my $locus (sort keys %additional_lengths) {
#  my $count = keys %{$additional_lengths{$locus}};
#  print "$locus: ";
#  for my $len (sort keys %{ $additional_lengths{$locus} }) {
#    print "$len=$additional_lengths{$locus}{$len} ";
#  }
#  print "~ count: $count";
#  print "\n";
#}

print "LOCUS,FWD_PRIM,PROBE,OT_RATIO,FWD-PROBE-REV,GHOST_LOCI,NUM_LENGTHS,MOST_FREQ_LENGTH,ALL_LENGTHS:FREQ\n";

foreach my $locus (sort keys %fwd_seq) {
  if ($fwd_count{$locus} == 0) {$OT_Ratio{$locus} = 0; }
  else {
    # calculate on target ratio (probe count / total read count per locus)
    $OT_Ratio{$locus} = $probe_count{$locus}/$fwd_count{$locus};
  }
  my $num_lengths = keys %{$additional_lengths{$locus}};
    
  print "$locus,$fwd_count{$locus},$probe_count{$locus},$OT_Ratio{$locus},$both_fwd_rev_count{$locus},$ghost_locus_count{$locus},$num_lengths,";

  # report all fragment lengths in order
  my $first = 1;
  if ($num_lengths == 0) { print "0,0,["; }
  else {
    foreach my $len (reverse sort { $additional_lengths{$locus}{$a} <=> $additional_lengths{$locus}{$b} or $a cmp $b} keys %{$additional_lengths{$locus}}) {
      my $count = $additional_lengths{$locus}{$len};
      if ($first) {
        $first = 0;
        print "$len,[";
      }
      print "$len:$count/";
    }
  }
  print "]\n";
}


### ORIGINAL ###

#print "LOCUS,FWD_PRIMER,PROBE,BOTH,OT_RATIO\n";
#foreach my $assays2 (sort keys %fwd_seq) {
#	if ($fwd_count{$assays2} == 0) {$fwd_count{$assays2} = 0.1}
#	$OT_Ratio{$assays2} = $probe_count{$assays2}/$fwd_count{$assays2};
#	print "$assays2,$fwd_count{$assays2},$probe_count{$assays2},$both_count{$assays2},$OT_Ratio{$assays2}\n";
#}
