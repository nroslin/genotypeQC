#!/usr/bin/perl -w

use strict;

# Created 10 July 2023.
# Last modified:  30 Oct 2023

# Make a final table of the QC stats generated per sample.

my $prefix = $ARGV[0];

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
my $output = "${prefix}_tableBySample.txt";

my ( %datahash, %sexhash, %infpop, %anchash, %failhash, %reasonhash );

#original clinical sex
open CLIN, "$inclin" or die "Cannot open file $inclin:  $!";
while (<CLIN>) {
  chomp;
  my ( $id0, $clinsex ) = split;
  $datahash{$id0} = $clinsex;
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
  #if ( $pedsex2 != $infsex && $pedsex2 != 0 ) { $sexhash{$id2} += 1; }
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
open ANC, "$inanc" or die "Cannot open file $inanc:  $!";
while (<ANC>) {
  chomp;
  next if /^FID/;
  my ( $id5, $pop ) = (split)[1,2];
  $datahash{$id5} = join "\t", $datahash{$id5}, $pop;
  $infpop{$id5} = $pop;   #inferred pop, to test against SR ancestry
}
close ANC;

#self-reported ancestry
#this info comes from clinical file, so will probably have more IDs than were
#genotyped
open SRA, "$insr" or die "Cannot open file $insr:  $!";
while (<SRA>) {
  chomp;
  my ( $id6, $sr ) = split;
  if ( exists $infpop{$id6} ) {
	$datahash{$id6} = join "\t", $datahash{$id6}, $sr;
	if ( $sr eq "NA" ) {
	  print "Warning:  ID $id6 has missing self-reported ancestry\n";
	  $anchash{$id6} = 0;   #don't call this a mismatch
	}
	elsif ( $infpop{$id6} ne $sr ) { $anchash{$id6} = 1; }   #mismatch
	else { $anchash{$id6} = 0; }
  }
}
close SRA;

### list of IDs failing QC ###
#put any sex mismatch in a separate hash from other QC failures
open REM, "$inrem" or die "Cannot open file $inrem:  $!";
while (<REM>) {
  chomp;
  next if /^FID/;
  my ( $id7, $reason ) = (split)[1,2];
  if ( $reason =~ /Sex/ ) { $sexhash{$id7} = 1; }
  else {   #failed for reasons other than sex
	$failhash{$id7} += 1;   #failed QC
	if ( exists $reasonhash{$id7} ) {
	  $reasonhash{$id7} = join ",", $reasonhash{$id7}, $reason;
	}
	else { $reasonhash{$id7} = $reason; }
  }
}
close REM;
  


#write it out, only for TAG samples
#### need to include column for sexhash problems
open OUT, ">$output" or die "Cannot write to file $output:  $!";
print OUT "FID\tIID\tClinSex\tXhet\tYcallRate\tGeneticSex\tSexChr\tCallRate\tAutoHet\tClosestAncestry\tSRancestry\tSexMismatch\tAncestryMismatch\tQCfail\tReasons\n";

open FAM, "$infam" or die "Cannot open file $infam:  $!";
while (<FAM>) {
  chomp;
  my ( $fid, $iid ) = (split)[0,1];
  next unless $iid =~ /^TAG/;
  if ( exists $datahash{$iid} ) {
	$sexhash{$iid} = 0 unless exists $sexhash{$iid};
	my $fail = 0;
	if ( exists $failhash{$iid} ) { $fail = 1; }
	my $why = "NA";
	if ( exists $reasonhash{$iid} ) { $why = $reasonhash{$iid}; }
	print OUT "$fid\t$iid\t$datahash{$iid}\t$sexhash{$iid}\t$anchash{$iid}\t$fail\t$why\n";
  }
  else { print "No QC results for $iid\n"; }
}
close FAM;
close OUT;