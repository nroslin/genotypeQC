#!/usr/bin/perl -w

use strict;

# Created 20 January 2022.
# Last modified:  25 Jan 2023

# Read in output of bcftools for chr X and Y data.  Super annoying to try in R.

my $inchrx = "$ARGV[0]_c23.bcfout";   #chrX, output of bcftools
my $inchry = "$ARGV[0]_c24.imiss";   #chrY missingness from PLINK
my $incrx = "$ARGV[0]_c23.imiss";   #chrX missingness from PLINK
my $inhety = "$ARGV[0]_c24.bcfout";   #chrY het
my $inped = "$ARGV[1]/$ARGV[2].fam";   #ped file with gender info
my $output = "$ARGV[2]_sexCheck.txt";
my ( %datahash, $ynsnp );


# chr X:  need to calculate heterozygosity
#id here is FID_IID
open CX, "$inchrx" or die "Cannot open file $inchrx:  $!";
while (<CX>) {
  chomp;
  next if /^#/;
  next if /^CMD/;
  next if /^DEF/;
  next if /^SITE0/;
  my ( $xid, $xnsnp, $xnhet ) = (split)[1,2,6];
  my $xhet;
  if ( $xnsnp == 0 ) {
	$xhet = 0;  #or should it be NA?  no calls on this chr
	print "No valid SNPs on chrX for $xid\n";
  }
  else { $xhet = sprintf "%.4f", $xnhet / $xnsnp; }
  $datahash{$xid} = join "\t", $xnsnp, $xnhet, $xhet;
}
close CX;


#chr Y, get call rate
#also need to attach FID so will match bcftools output
open CY, "$inchry" or die "Cannnot open file $inchry:  $!";
while (<CY>) {
  chomp;
  next if /FID/;
  my ( $fid, $yid, $ynmiss, $ynsnp, $fmiss ) = (split)[0,1,3,4,5];
  my $sample = join "_", $fid, $yid;
  my $yosnp = $ynsnp - $ynmiss;   #number of non-missing genotypes
  my $ycr = sprintf "%.4f", 1 - $fmiss;   #call rate
  $datahash{$sample} = join "\t", $datahash{$sample}, $ynsnp, $yosnp, $ycr;
}
close CY;

#chr X, get call rate
#not sure yet how helpful this will be
if ( -e $incrx ) {
  open CRX, "$incrx" or die "Cannot open file $incrx:  $!";
  while (<CRX>) {
	chomp;
	next if /FID/;
	my ( $cfid, $ciid, $cnmiss, $cnsnp, $cfmiss ) = (split)[0,1,3,4,5];
	my $xsample = join "_", $cfid, $ciid;
	#my $xosnp = $cnsnp - $cnmiss;   #number of non-missing genotypes
	my $xcr = sprintf "%.4f", 1 - $cfmiss;   #call rate
	$datahash{$xsample} = join "\t", $datahash{$xsample}, $xcr;
  }
}

#chrY het
#id here is FID_IID
open YHET, "$inhety" or die "Cannot open file $inhety:  $!";
while (<YHET>) {
  chomp;
  next if /^#/;
  next if /^CMD/;
  next if /^DEF/;
  next if /^SITE0/;
  my ( $yhetid, $yhnsnp, $ynhet ) = (split)[1,2,6];
  my $yhet;
  if ( $yhnsnp == 0 ) {
	$yhet = 0;  #should I make it NA instead?
	print "No valid SNPs on chrY for $yhetid\n";
  }
  else { $yhet = sprintf "%.4f", $ynhet / $yhnsnp; }
  $datahash{$yhetid} = join "\t", $datahash{$yhetid}, $yhet;
}
close YHET;


#print everything out
### assume that if have both or none of chrX call rate and chrY het
open OUT, ">$output" or die "Cannot write to file $output:  $!";
print OUT "FID\tIID\tPedSex\tChrXnSNP\tChrXnHet\tChrXhetRate\tChrYnSNP\tChrYnCalled\tChrYcr\tChrXcr\tChrYhet\n";

open PED, "$inped" or die "Cannot open file $inped:  $!";
while (<PED>) {
  chomp;
  my ( $fid, $iid, $sex ) = (split)[0,1,4];
  my $label = join "_", $fid, $iid;
### add NA if X chr call rate file not present
  my @array = split $datahash{$label};
  if ( scalar @array == 6 ) {   #if no chrX call rate or chrY het
	print OUT "$fid\t$iid\t$sex\t$datahash{$label}\tNA\tNA\n";
  }
  else { print OUT "$fid\t$iid\t$sex\t$datahash{$label}\n"; }
}
close PED;
close OUT;
