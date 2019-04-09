#! perl

use strict;
use warnings;
use Cwd;
use File::Basename;

##impose filters on Pisces VCF format
##this outputs a 'raw.vcf' and 'pass.snv.vcf', 'pass.indel.vcf'

##inputs:
##-ID = Tumour ID (VCF must be single sample)
##-DP = Required depth
##-MD = Required min depth of ALT
##-VCF = VCF input file (in current wd or full path)

##strategy: Tumour col is [-1]
##make Tumour subs, input and done
##print to FILTER new PASS_BMFILT

my ($no,$id,$dp,$md,$vcf)="";
if(scalar(@ARGV)!=4){
  print "Require 4 inputs:\n-ID=Tumour ID (VCF must be single sample)\n-DP=Required depth (>=)\n-MD=Required min depth of ALT (>=)\n-VCF=VCF input file (in current dir or full path)\n";
  exit;
}

for(my $i=0;$i<@ARGV;$i++){
  if($ARGV[$i]=~m/ID/){($no,$id)=split(/\=/,$ARGV[$i]);}
  if($ARGV[$i]=~m/DP/){($no,$dp)=split(/\=/,$ARGV[$i]);}
  if($ARGV[$i]=~m/MD/){($no,$md)=split(/\=/,$ARGV[$i]);}
  if($ARGV[$i]=~m/VCF/){($no,$vcf)=split(/\=/,$ARGV[$i]);}
}

if($vcf!~m/\//){
  my $dir=getcwd();
  my $vcf1=$dir . "/" . $vcf;
  $vcf=$vcf1;
}

##open VCF, run filter
my $tcol=-1;
my $header="";
my $snv="";
my $indel="";
my $raw="";

open(VCF,$vcf);
while(<VCF>){
  chomp;
  my @sp=split(/\t/);

  ##header
  if($_=~m/^#/){
    $header.=$_ . "\n";
    next;
  }

  ##variants
  else{
    $raw.=$_ . "\n";
    if($sp[6] eq "PASS"){
      my $tflag=&tumourFilter(@sp);
      my $iflag=&indelFilter(@sp);
      my $siflag=&indelOrSNV(@sp);
      if($siflag==0){
        if($tflag == 1){
          $snv.=$_ . "\n";
          next;
        }
      }
      else{
        if(($tflag == 1) && ($iflag == 1)){
          $indel.=$_ . "\n";
          next;
        }
      }
    }
  }
}

my $outDir=dirname($vcf);
my $outName=$outDir . "/" . $id . ".pisces";
my $snvName=$outName . ".snv.pass.vcf";
open(OUT,">$snvName");
print OUT $header;
print OUT $snv;
close OUT;

my $indelName=$outName . ".indel.pass.vcf";
open(OUT,">$indelName");
print OUT $header;
print OUT $indel;
close OUT;

my $rawName=$outName . ".raw.vcf";
open(OUT,">$rawName");
print OUT $header;
print OUT $raw;
close OUT;

##subroutines
sub tumourFilter {
  my @info=split(/\:/,$_[$tcol]);
  my @refalt=split(/\,/,$info[1]);
  my $tot=$refalt[0]+$refalt[1];
  if(($tot > $dp) && ($refalt[1] > $md)){
    return(1);
  }
  else{
    return(0);
  }
}

sub indelFilter {
  my @info=split(/\:/,$_[$tcol]);
  my @refalt=split(/\,/,$info[1]);
  my $tot=$refalt[0]+$refalt[1];
  if(($tot > $dp) && ($refalt[1] < $md)){
    return(1);
  }
  else{
    return(0);
  }
}

sub indelOrSNV {
  my @ref=split(//,$_[3]);
  my @alt=split(//,$_[4]);
  if((scalar(@ref) > 1) || (scalar(@alt) > 1)){
      return(1);
  }
  else{
      return(0);
  }
}
