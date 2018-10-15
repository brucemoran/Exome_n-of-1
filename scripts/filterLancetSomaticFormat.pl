#! /bin/perl5

use strict;
use warnings;
use Cwd;
use File::Basename;

##Lancet VCF format does not contain FORMAT AF
##inputs:
##-ID = Tumour ID (VCF must be somatic, 2 sample, Tumour ID -> other is Germline)
##-DP = Required depth Tumour and Germline
##-MD = Required min depth of ALT in Tumour
##-VCF = VCF input file (in current wd or full path)
##NB hardcoded for 0 ALT in Germline

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
my $filtflag=0;
my $header="";
my $snv="";
my $indel="";
my $raw="";

open(VCF,$vcf);
while(<VCF>){
  chomp;
  my @sp=split(/\t/);

  if($_=~m/^##/){
    if($_=~m/FORMAT/){
      if($filtflag==0){
        $filtflag++;
        $header.="##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n##FORMAT=<ID=AF,Number=1,Type=String,Description=\"Allele Frequency\">\n##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"allele depth: # of supporting ref,alt reads at the site\">\n##FORMAT=<ID=SR,Number=.,Type=Integer,Description=\"strand counts for ref: # of supporting forward,reverse reads for reference allele\">\n##FORMAT=<ID=SA,Number=.,Type=Integer,Description=\"strand counts for alt: # of supporting forward,reverse reads for alterantive allele\">\n##FORMAT=<ID=DP,Number=1,Type=String,Description=\"Total Depth\">\n";
        next;
      }
      else{
        next;
      }
    }
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

  ##down to pay-dirt
  else{
    my $lastthree=&AFlancet(@sp);
    my $out="";
    if($lastthree eq 0){
      next;
    }
    else{
      my @popd=pop(@sp);
      @popd=pop(@sp);
      @popd=pop(@sp);
      $out.=join("\t",@sp[0..$#sp]);
      $out.="\t$lastthree\n";
      @sp=split(/\t/,$out);
      $raw.=$out;

      if($sp[6] eq "PASS"){
        my $tflag=&tumourFilter(@sp);
        my $gflag=&germlineFilter(@sp);
        my $iflag=&indelgermlineFilter(@sp);
        my $siflag=&indelOrSNV(@sp);
        if($siflag==0){
          if(($tflag == 1) && ($gflag == 1)){
            $snv.=$out;
            next;
          }
        }
        else{
          if(($tflag == 1) && ($iflag == 1)){
            $indel.=$out;
            next;
          }
        }
      }
    }
  }
}

my $outDir=dirname($vcf);
my $outName=$outDir . "/" . $id . ".lancet";
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
  if($info[0] eq "."){return(0);}
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
  if($info[0] eq "."){return(0);}
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

sub AFlancet {
  my $gaf="";
  my $taf="";
  my $out="";
  my @tf=split(/\:/,$_[$tcol]);
  my @gf=split(/\:/,$_[$gcol]);
  if(($tf[0] eq "0/0") || ($tf[0] eq ".") || ($gf[0] eq ".")){
    return(0);
  }
  else{

    ##REF, ALT for T, G
    my @tra=split(/\,/,$tf[1]);
    my @gra=split(/\,/,$gf[1]);
    my $tt=$tra[0]+$tra[1];
    my $taff=$tra[1]/$tt;
    my $taf=sprintf("%.3f",$taff);
    my $gt=$gra[0]+$gra[1];
    my $gaff=$gra[1]/$gt;
    my $gaf=sprintf("%.3f",$gaff);

    $out="GT:AD:AF:SR:SA:DP\t";
    if($tcol>$gcol){
      $out.="$gf[0]:$gf[1]:$gaf:" . join(":",@gf[2..$#gf]);
      $out.="\t$tf[0]:$tf[1]:$taf:" . join(":",@tf[2..$#tf]);
    }
    else{
      $out.="$tf[0]:$tf[1]:$taf:" . join(":",@tf[2..$#tf]);
      $out.="\t$gf[0]:$gf[1]:$gaf:" . join(":",@gf[2..$#gf]);
    }
  }
  return($out);
}
