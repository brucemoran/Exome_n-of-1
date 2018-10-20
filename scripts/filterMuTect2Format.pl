#! perl

use strict;
use warnings;
use Cwd;
use File::Basename;

##impose filters on MuTect2 VCF format
##this outputs a 'raw.vcf' and 'pass.snv.vcf', 'pass.indel.vcf'

##inputs:
##-ID = Tumour ID (VCF must be somatic, 2 sample, Tumour ID -> other is Germline)
##-DP = Required depth Tumour and Germline
##-MD = Required min depth of ALT in Tumour
##-VCF = VCF input file (in current wd or full path)
##NB hardcoded for 0 ALT in Germline

##strategy: find Tumour col, either [-1] or [-2]
##make Germline, Tumour subs, input and done
##print to FILTER new PASS_BMFILT

my ($no,$id,$dp,$md,$vcf)="";
if(scalar(@ARGV)!=4){
  print "Require 4 inputs:\n-ID=Tumour ID (VCF must be somatic, 2 sample, Tumour ID -> other is Germline)\n-DP=Required depth Tumour and Germline (>=)\n-MD=Required min depth of ALT in Tumour (>=)\n-VCF=VCF input file (in current dir or full path)\n";
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
my $gcol=-2;
my $header="";
my $snv="";
my $indel="";
my $raw="";

open(VCF,$vcf);
while(<VCF>){
  chomp;
  my @sp=split(/\t/);

  ##header
  if($_=~m/^##/){
    $header.=$_ . "\n";
    next;
  }

  if($_=~m/^#CHROM/){
    if($sp[-1] eq $id){
      $header.=$_ . "\n";
      next;
    }
    else{
      $tcol=-2;$gcol=-1;
      $header.=$_ . "\n";
      next;
    }
  }

  ##variants
  else{
    $raw.=$_ . "\n";
    if(($sp[6] eq "PASS") || ($sp[6] eq ".")){
      my $tflag=&tumourFilter(@sp);
      my $gflag=&germlineFilter(@sp);
      my $iflag=&indelgermlineFilter(@sp);
      my $siflag=&indelOrSNV(@sp);
      if($siflag==0){
        if(($tflag == 1) && ($gflag == 1)){
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
my $outName=$outDir . "/" . $id . ".mutect2";
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

sub germlineFilter {
  my @info=split(/\:/,$_[$gcol]);
  my @refalt=split(/\,/,$info[1]);
  my $tot=$refalt[0]+$refalt[1];
  if(($tot > $dp) && ($refalt[1] == 0)){
    return(1);
  }
  else{
    return(0);
  }
}

sub indelgermlineFilter {
  my @info=split(/\:/,$_[$gcol]);
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
