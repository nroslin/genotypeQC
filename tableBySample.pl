#!/usr/bin/perl -w

use strict;

# Created 10 July 2023.
# Last modified:  20 Dec 2023

# Make a final table of the QC stats generated per sample.

my $prefix = $ARGV[0];
my $batch = $ARGV[1];  #not consistent, so need to read it in

my $inclin = "../../recoded/${prefix}.providedSex.txt";
	#we will see if this hard-coded file name works


my $insex1 = "${prefix}_sexCheck.txt";  #raw stats
my $insex2 = "${prefix}_inferredSex.txt";   #inference
my $inmiss = "${prefix}_missing_step2.imiss";   #missingness
my $inhet = "${prefix}_autoHet.txt";   #heterozygosity
my $inanc = "${prefix}.1kgpca.closestAncestry.txt";   #pop inference
my $insr = "${prefix}_SRancestry.txt";   #self-reported ancestry
my $inrem = "${prefix}.remove.txt";   #IDs failing QC
my $infam = "../../recoded/${prefix}.fam";    #genotyped samples
my $inplate = "${prefix}_idsWithPlate.txt";  #plate labels for samples
my $output = "${prefix}_tableBySample.txt";

my ( %datahash, %sexfail, %srsex, %sranc );
my ( %anchash, %failhash, %reasonhash );

#original clinical sex
open CLIN, "$inclin" or die "Cannot open file $inclin:  $!";
while (<CLIN>) {
  chomp;
  my ( $id0, $clinsex ) = split;
  $datahash{$id0} = $clinsex;
  $srsex{$id0} = $clinsex;   #SR sex, to be compared to inferred sex
}
close CLIN;


### sex inference ###

#get chr X het and chr Y call rate
open SEX1, "$insex1" or die "Cannot open file $insex1:  $!";
while (<SEX1>) {
  chomp;
  next if /^FID/;
  my ( $id1, $xhet, $ycr ) = (split)[1,5,8];
  if ( exists $datahash{$id1} ) {
    $datahash{$id1} = join "\t", $datahash{$id1}, $xhet, $ycr;
  }
  else { $datahash{$id1} = join "\t", "NA", $xhet, $ycr; }
}
close SEX1;

#get inferred sex and chromosome counts
open SEX2, "$insex2" or die "Cannot open file $insex2:  $!";
while (<SEX2>) {
  chomp;
  next if /^FID/;
  my ( $id2, $pedsex2, $infsex, $nx, $ny ) = (split)[1,2,3,4,5];
  my $sexchr;
  if ( $nx eq "NA" || $ny eq "NA" ) { $sexchr = "unknown"; }
  elsif ( $nx == 1 && $ny == 1 ) { $sexchr = "XY"; }
  elsif ( $nx == 2 && $ny == 0 ) { $sexchr = "XX"; }
  elsif ( $nx == 2 && $ny == 1 ) { $sexchr = "XXY"; }   #add other options?
  else { $sexchr = "unknown"; }
  $datahash{$id2} = join "\t", $datahash{$id2}, $infsex, $sexchr;

  #IDs where inferred sex ne pedsex and pedsex not missing/unknown
  #(the data in $pedsex2 comes from the genotyping lab and is not accurate)
  my $selfsex;
  if ( exists $srsex{$id2} ) {
	$selfsex = $srsex{$id2};   #sex from self-reported clinical data
  }
  elsif ( $id2 !~ /^TAG/ ) {   #non-TAG, assume is NA12878 (female)
	$selfsex = 2;
  }
  else { print "No SR sex for $id2\n"; $selfsex = 0; }

  #compare self-reported sex to inferred sex
  if ( $infsex eq "NA" ) { $sexfail{$id2} = 1; }
  elsif ( $selfsex eq "NA" ) { next; }  #do nothing if SR sex missing
  elsif ( $selfsex == 1 && $infsex == 2 ) { $sexfail{$id2} = 1; }
  elsif ( $selfsex == 2 && $infsex == 1 ) { $sexfail{$id2} = 1; }
	#fail if inferred sex is ambiguous
	#OR if SR sex is 1/2 (don't fail if non-standard response)
	#AND inferred sex is also known
}
close SEX2;


### call rate ###
open MISS, "$inmiss" or die "Cannot open file $inmiss:  $!";
while (<MISS>) {
  chomp;
  next if /FID/;
  my ( $id3, $fmiss ) = (split)[1,5];   #missing rate
  my $cr = sprintf "%.3f", 1 - $fmiss;
  $datahash{$id3} = join "\t", $datahash{$id3}, $cr;
}
close MISS;

### heterozygosity ###
open HET, "$inhet" or die "Cannot open file $inhet:  $!";
while (<HET>) {
  chomp;
  next if /Filter/;
  my ( $id4, $ngt, $nhet ) = (split)[2,3,7];
  my $fhet = sprintf "%.3f", $nhet / $ngt;
  $datahash{$id4} = join "\t", $datahash{$id4}, $fhet;
}
close HET;

### ancestry ###

#self-reported ancestry
#this info comes from clinical file, so will probably have more IDs than were
#genotyped
open SRA, "$insr" or die "Cannot open file $insr:  $!";
while (<SRA>) {
  chomp;
  my ( $id6, $sr ) = split;
  if ( exists $datahash{$id6} ) {
	$datahash{$id6} = join "\t", $datahash{$id6}, $sr;
	$sranc{$id6} = $sr;   #reported ancestry, to test against inferred
  }
#	if ( $sr eq "NA" ) {
#	  print "Warning:  ID $id6 has missing self-reported ancestry\n";
#	  $anchash{$id6} = 0;   #don't call this a mismatch
#	}
#	elsif ( $infpop{$id6} ne $sr ) { $anchash{$id6} = 1; }   #mismatch
#	else { $anchash{$id6} = 0; }
#  }
}
close SRA;

#inferred ancestry
open ANC, "$inanc" or die "Cannot open file $inanc:  $!";
while (<ANC>) {
  chomp;
  next if /^FID/;
  my ( $id5, $pop ) = (split)[1,2];
  $datahash{$id5} = join "\t", $datahash{$id5}, $pop;

  #test reported vs. inferred
  if ( $id5 !~ /^TAG/ ) { $sranc{$id5} = "EUR"; }   #these should be controls
  if ( $sranc{$id5} eq "NA" ) {
	print "Warning:  ID $id5 has missing self-reported ancestry\n";
	$anchash{$id5} = 0;   #don't call this a mismatch
  }
  elsif ( $sranc{$id5} ne $pop ) { $anchash{$id5} = 1; }  #mismatch
  else { $anchash{$id5} = 0; }
}
close ANC;

#### plate and batch labels
open PL, "$inplate" or die "Cannot open file $inplate:  $!";
while (<PL>) {
  chomp;
  my ( $id9, $plate ) = split;
  if ( exists $datahash{$id9} ) {
	$datahash{$id9} = join "\t", $datahash{$id9}, $plate, $batch;
  }
  else { print "Warning:  No QC stats for genotyped ID $id9\n"; }
}
close PL;

### list of IDs failing QC ###
open REM, "$inrem" or die "Cannot open file $inrem:  $!";
while (<REM>) {
  chomp;
  next if /^FID/;
  my ( $id7, $reason ) = (split)[1,2];
  $failhash{$id7} += 1;   #failed QC
  if ( exists $reasonhash{$id7} ) {
	$reasonhash{$id7} = join ",", $reasonhash{$id7}, $reason;
  }
  else { $reasonhash{$id7} = $reason; }
}
close REM;

#if failed QC and also has sex or ancestry mismatch, add sex/ancestry to list
#of reasons
foreach my $id8 ( keys %failhash ) {
  if ( exists $sexfail{$id8} ) {
	$reasonhash{$id8} = join ",", $reasonhash{$id8}, "SEX";
  }
  if ( exists $anchash{$id8} && $anchash{$id8} == 1 ) {
	$reasonhash{$id8} = join ",", $reasonhash{$id8}, "ANCESTRY";
  }   #samples not in %anchash are controls
}
  


#write it out, only for TAG samples
open OUT, ">$output" or die "Cannot write to file $output:  $!";
print OUT "family_id\tparticipant_id\tReportedSex\tXhet\tYcallRate\tGeneticSex\tSexChr\tCallRate\tAutoHet\tReportedAncestry\tClosestAncestry\tGtPlate\tGtBatch\tSexMismatch\tAncestryMismatch\tQCfail\tReasons\n";

open FAM, "$infam" or die "Cannot open file $infam:  $!";
while (<FAM>) {
  chomp;
  my ( $fid, $iid ) = (split)[0,1];
  next unless $iid =~ /^TAG/;
  if ( exists $datahash{$iid} ) {
	$sexfail{$iid} = 0 unless exists $sexfail{$iid};
	my $fail = 0;
	if ( exists $failhash{$iid} ) { $fail = 1; }
	my $why = "NA";
	if ( exists $reasonhash{$iid} ) { $why = $reasonhash{$iid}; }
	print OUT "$fid\t$iid\t$datahash{$iid}\t$sexfail{$iid}\t$anchash{$iid}\t$fail\t$why\n";
  }
  else { print "No QC results for $iid\n"; }
}
close FAM;
close OUT;
