#!/usr/bin/perl -w

use strict;

# Created 20 January 2022.
# Last modified:  21 Jan 2022

# Read in output of bcftools for chr X and Y data.  Super annoying to try in R.

my $inchrx = "$ARGV[0]_c23.bcfout";   #chrX, output of bcftools
my $inchry = "$ARGV[0]_c24.bcfout";   #chrY, output from bcftoos
my $inped = "$ARGV[1]/$ARGV[2].fam";   #ped file with gender info
my $inlog = "$ARGV[0]_c24.log";   #PLINK log file when created chrY vcf file
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

#chr Y:  need call rate, so start with total number of SNPs
open LOG, "$inlog" or die "Cannot open file $inlog:  $!";
while (<LOG>) {
  chomp;
  next unless /pass/;   #<blank> variants and <blank> people pass filters...
  $ynsnp = (split)[0];
}
close LOG;

#chr Y, get call rate
#ID here is FID_IID
open CY, "$inchry" or die "Cannnot open file $inchry:  $!";
while (<CY>) {
  chomp;
  next if /^#/;
  next if /^CMD/;
  next if /^DEF/;
  next if /^SITE0/;
  my ( $yid, $yosnp ) = (split)[1,2];  #observed count
  my $ycr = sprintf "%.4f", $yosnp / $ynsnp;
  $datahash{$yid} = join "\t", $datahash{$yid}, $ynsnp, $yosnp, $ycr;
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
