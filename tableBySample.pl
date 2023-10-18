#!/usr/bin/perl -w

use strict;

# Created 10 July 2023.
# Last modified:  02 Aug 2023

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
my $intwin = "${prefix}_mztwins.txt";   #any twins/duplicates
my $inrem = "${prefix}.remove.txt";   #IDs failing QC
my $infam = "../../recoded/${prefix}.fam";    #genotyped samples
my $output = "${prefix}_tableBySample.txt";

my ( %datahash, %sexhash, %srpop, %anchash, %twinhash, %twinidhash );
my ( %failhash, %reasonhash );

#original clinical sex
open CLIN, "$inclin" or die "Cannot open file $inclin:  $!";
while (<CLIN>) {
  chomp;
  my ( $id0, $clinsex ) = split;
  $datahash{$id0} = $clinsex;
}
close CLIN;
#control samples (not in clinical database) should be NA12878, female EUR



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
  else { 
	$datahash{$id1} = join "\t", 2, $xhet, $ycr; 
	print "Warning:  ID $id1 has genotypes but is not in the clinical file\n";
  }
	#assume anyone genotyped but not in the clinical database are the
	#controls, who should be female (NA12878)
	#but also give warning
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
#self-reported ancestry
#this info comes from clinical file, so will probably have more IDs than were
#genotyped
open SRA, "$insr" or die "Cannot open file $insr:  $!";
while (<SRA>) {
  chomp;
  my ( $id5, $sr ) = split;
  $srpop{$id5} = $sr;   #to test against inferred
  if ( exists $datahash{$id5} ) {
	$datahash{$id5} = join "\t", $datahash{$id5}, $sr;
  }
}
close SRA;

open ANC, "$inanc" or die "Cannot open file $inanc:  $!";
#inferred ancestry
while (<ANC>) {
  chomp;
  next if /^FID/;
  my ( $id6, $pop ) = (split)[1,2];
  if ( exists $srpop{$id6} ) {   #in clinical file
	$datahash{$id6} = join "\t", $datahash{$id6}, $pop;
	if ( $srpop{$id6} ne $pop ) { $anchash{$id6} = 1; }  #mismatch
	else { $anchash{$id6} = 0; }
  }
  else {   #not in clinical file, hopefully a control ("SR" EUR)
	$datahash{$id6} = join "\t", $datahash{$id6}, "EUR", $pop;
	if ( $pop ne "EUR" ) { $anchash{$id6} = 1; }
	else { $anchash{$id6} = 0; }
  }
	
}
close ANC;


### twins/duplicates ###
#this is a work in progress
open TWIN, "$intwin" or die "Cannot open file $intwin:  $!";
while (<TWIN>) {
  chomp;
  my ( $name, $tid ) = split;   #$name is FID_IID, $tid is unique twin ID
  my ( $fid7, $id7 ) = split /_/, $name;
#  if ( exists $twinhash{$tid} ) {    #get list of TAG ids for each tid
#	$twinhash{$tid} = join ",", $twinhash{$tid}, $id7;
#  }
#  else { $twinhash{$tid} = $id7; }  #not sure how helpful this is, since tid
					#only exists here

  if ( exists $twinhash{$tid} ) { #find matching twin pairs
	my $othertwin = $twinhash{$tid}; 
	$twinidhash{$id7} = $othertwin;
	$twinidhash{$othertwin} = $id7;
  }
  else { $twinhash{$tid} = $id7; }
  #this won't work properly for larger clusters (eg, triplets)
}
close TWIN;


### list of IDs failing QC ###
#put any sex mismatch in a separate hash from other QC failures
open REM, "$inrem" or die "Cannot open file $inrem:  $!";
while (<REM>) {
  chomp;
  next if /^FID/;
  my ( $id8, $reason ) = (split)[1,2];
  if ( $reason =~ /Sex/ ) { $sexhash{$id8} = 1; }
  else {   #failed for reasons other than sex
	$failhash{$id8} += 1;   #failed QC
	if ( exists $reasonhash{$id8} ) {
	  $reasonhash{$id8} = join ",", $reasonhash{$id8}, $reason;
	}
	else { $reasonhash{$id8} = $reason; }
  }
}
close REM;
  

#write it out, only for TAG samples; no, include controls but flag
#### need to include column for sexhash problems
open OUT, ">$output" or die "Cannot write to file $output:  $!";
print OUT "FID\tIID\tClinSex\tXhet\tYcallRate\tGeneticSex\tSexChr\tCallRate\tAutoHet\tSRancestry\tClosestAncestry\tSexMismatch\tAncestryMismatch\tQCfail\tReasons\n";

open FAM, "$infam" or die "Cannot open file $infam:  $!";
while (<FAM>) {
  chomp;
  my ( $fid, $iid ) = (split)[0,1];
#  next unless $iid =~ /^TAG/;
  if ( exists $datahash{$iid} ) {
	$sexhash{$iid} = 0 unless exists $sexhash{$iid};
	my $fail = 0;
	if ( exists $failhash{$iid} ) { $fail = 1; }
	#add non-TAG samples to list of those who fail
	#unless ( $iid =~ /^TAG/ ) { $fail = 1; }
	my $why = "NA";
	if ( exists $reasonhash{$iid} ) { 
	  $why = $reasonhash{$iid}; 
	  if ( $why =~ /TWIN/ ) {
		$why = join ",", $why, $twinidhash{$iid};
	  }
	}
	print OUT "$fid\t$iid\t$datahash{$iid}\t$sexhash{$iid}\t$anchash{$iid}\t$fail\t$why\n";
  }
  else { print "No QC results for $iid\n"; }
}
close FAM;
close OUT;
