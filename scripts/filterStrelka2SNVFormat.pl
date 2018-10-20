#! perl

use strict;
use warnings;
use Cwd;
use File::Basename;

##impose filters on Strelka2 VCF format
##inputs:
##-ID = Tumour ID (VCF must be somatic, 2 sample, Tumour ID -> other is Germline)
##-DP = Required depth Tumour and Germline
##-MD = Required min depth of ALT in Tumour
##-VCF = VCF input file (in current wd or full path)
##NB hardcoded for 0 ALT in Germline

#per sample:
#DP:FDP:SDP:SUBDP:AU:CU:GU:TU
#Depth : filtered : deletions_spanning : below_toer1_map_quality : A : C : G : T

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
my $filtflag=0;
my $header="";
my $snv="";
my $raw="";

open(VCF,$vcf);
while(<VCF>){
  chomp;
  my @sp=split(/\t/);

  if($_=~m/^##/){
    if($_=~m/FORMAT/){
      if($filtflag==0){
        $filtflag++;
        $header.="##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GENOTYPE\">\n";
        $header.="##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"ALELE DEPTHS\">\n";
        $header.="##FORMAT=<ID=AF,Number=1,Type=Float,Description=\"ALLELE FREQUENCY (ADalt/DP)\">\n";
        $header.=$_ . "\n";
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
    my $line="";
    my $lasttwo=&GTADAFStrelka2(@sp);
    if($lasttwo eq 0){
      next;
    }
    else{
      my @popd=pop(@sp);
      @popd=pop(@sp);
      my $pref="GT:AD:AF:";
      my $sp1=$pref . $sp[-1];
      $sp[-1]=$sp1;
      $line.=join("\t",@sp[0..$#sp]);
      $line.="\t$lasttwo\n";
      $raw.=$line;
      @sp=split(/\t/,$line);
      if($sp[6] eq "PASS"){
        my $tflag=&tumourFilter(@sp);
        my $gflag=&germlineFilter(@sp);
        if(($tflag == 1) && ($gflag == 1)){
          my $siflag=&indelOrSNV(@sp);
          if($siflag==0){
            $snv.=$line;
            next;
          }
        }
      }
      else{
        next;
      }
    }
  }
}

my $outDir=dirname($vcf);
my $outName=$outDir . "/" . $id . ".strelka2";
my $snvName=$outName . ".snv.pass.vcf";
print "Output to $snvName\n";
open(OUT,">$snvName");
print OUT $header;
print OUT $snv;
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
  if(($refalt[0] eq ".") || ($refalt[1] eq ".")){
    return(0);
  }
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

sub findBaseInfo {
  ##which index for base in info?
  my $index="";
  if($_[0] eq "A"){$index=-4}
  if($_[0] eq "C"){$index=-3}
  if($_[0] eq "G"){$index=-2}
  if($_[0] eq "T"){$index=-1}
  return($index);
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

sub GTADAFStrelka2 {
  my $ggt="";
  my $tgt="";
  my @tf=split(/\:/,$_[$tcol]);
  if($tf[0]==0){
    return(0);
  }
  else{
    my @gf=split(/\:/,$_[$gcol]);
    my $ref=&findBaseInfo($_[3]);
    my $alt=&findBaseInfo($_[4]);

    ##REF, ALT for T, G
    my @tr=split(/\,/,$tf[$ref]);
    my @ta=split(/\,/,$tf[$alt]);
    my @gr=split(/\,/,$gf[$ref]);
    my @ga=split(/\,/,$gf[$alt]);

    ##make GTs
    ##"A,T,11,2" -> 0/1:11,2:0.154:13
    ##"C,G,11,0" -> 0/0:11,0:0:11

    ##germline 0/0
    if($ga[0]==0){
      $ggt="0/0:$gr[0],0:0:";
    }
    ##tumour 0/1
    if(($tr[0]>0) && ($ta[0]>0)){
      my $tot=$tr[0]+$ta[0];
      my $fq=$ta[0]/$tot;
      my $fqr=sprintf("%.3f", $fq);
      $tgt="0/1:$tr[0],$ta[0]:$fqr:";
    }
    ##germline 0/1
    if(($gr[0]>0) && ($ga[0]>0)){
      my $tot=$gr[0]+$ga[0];
      my $fq=$ga[0]/$tot;
      my $fqr=sprintf("%.3f", $fq);
      $ggt="0/1:$gr[0],$ga[0]:$fqr:";
    }
  }
  if(($tgt eq "") || ($ggt eq "")){
    return(0);
  }
  elsif($tcol eq -2){
    my $out="$tgt$_[$tcol]\t$ggt$_[$gcol]";
    return($out);
  }
  else{
    my $out="$ggt$_[$gcol]\t$tgt$_[$tcol]";
    return($out);
  }
}
