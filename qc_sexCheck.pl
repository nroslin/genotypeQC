#!/usr/bin/perl -w

use strict;

# Created 20 January 2022.
# Last modified:  12 Jan 2023

# Read in output of bcftools for chr X and Y data.  Super annoying to try in R.

my $inchrx = "$ARGV[0]_c23.bcfout";   #chrX, output of bcftools
my $inchry = "$ARGV[0]_c24.imiss";   #chrY missingness from PLINK
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
  my $xhet = sprintf "%.4f", $xnhet / $xnsnp;
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

#print everything out
open OUT, ">$output" or die "Cannot write to file $output:  $!";
print OUT "FID\tIID\tPedSex\tChrXnSNP\tChrXnHet\tChrXhetRate\tChrYnSNP\tChrYnCalled\tChrYcr\n";

open PED, "$inped" or die "Cannot open file $inped:  $!";
while (<PED>) {
  chomp;
  my ( $fid, $iid, $sex ) = (split)[0,1,4];
  my $label = join "_", $fid, $iid;
  print OUT "$fid\t$iid\t$sex\t$datahash{$label}\n";
}
close PED;
close OUT;
