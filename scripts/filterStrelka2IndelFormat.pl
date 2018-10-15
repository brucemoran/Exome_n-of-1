#! /bin/perl5

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
        $header.="##FORMAT=<ID=GT,Number=1,Type=String,Description=\"GENOTYPE\">\n";
        $header.="##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"ALLELE DEPTHS\">\n";
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
    my $lasttwo=&GTADAFStrelka2Indel(@sp);
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
      @sp=split(/\t/,$line);
      $raw.=$line;
      if($sp[6] eq "PASS"){
        my $tflag=&tumourFilter(@sp);
        my $iflag=&indelgermlineFilter(@sp);
        my $siflag=&indelOrSNV(@sp);
        if(($siflag==1) && ($tflag == 1) && ($iflag == 1)){
            $indel.=$line;
          }
        else{
          next;
        }
      }
    }
  }
}

my $outDir=dirname($vcf);
my $outName=$outDir . "/" . $id . ".strelka2";
my $indelName=$outName . ".indel.pass.vcf";
open(OUT,">$indelName");
print OUT $header;
print OUT $indel;
close OUT;

my $rawName=$outName . ".raw.vcf";
open(OUT,">>$rawName");
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

sub GTADAFStrelka2Indel {

  ##make GT:AD:AF
  my $ggt="";
  my $tgt="";

  my @ti=split(/\:/,$_[$tcol]);
  my @gi=split(/\:/,$_[$gcol]);
  my @tref=split(/\,/,$ti[2]);
  my @talt=split(/\,/,$ti[3]);
  my @gref=split(/\,/,$gi[2]);
  my @galt=split(/\,/,$gi[3]);

  ##make GTs
  ##"A,T,11,2" -> 0/1:11,2:0.154:13
  ##"C,G,11,0" -> 0/0:11,0:0:11

  ##germline 0/0
  if($galt[0]==0){
    my $ggt="0/0:$gref[0],0:0:";
  }
  ##germline 0/1
  if(($gref[0]>0) && ($galt[0]>0)){
    my $tot=$gref[0]+$galt[0];
    my $fq=$galt[0]/$tot;
    my $fqr=sprintf("%.3f", $fq);
    $ggt="0/1:$gref[0],$galt[0]:$fqr:";
  }
  ##germline 1/1
  if(($gref[0]==0) && ($galt[0]!=0)){
    $ggt="1/1:0,$galt[0]:1:";
  }
  ##tumour 0/0
  if($talt[0]==0){
    my $tgt="0/0:$tref[0],0:0:";
  }
  ##tumour 0/1
  if(($tref[0]>0) && ($talt[0]>0)){
    my $tot=$tref[0]+$talt[0];
    my $fq=$talt[0]/$tot;
    my $fqr=sprintf("%.3f", $fq);
    $tgt="0/1:$tref[0],$talt[0]:$fqr:";
  }
  ##tumour 1/1
  if(($tref[0]==0) && ($talt[0]!=0)){
    $tgt="1/1:0,$talt[0]:0:";
  }

  if(($tgt eq "") || ($ggt eq "")){
    return(0);
  }
  if($tcol == -2){
    my $out="$tgt$_[$tcol]\t$ggt$_[$gcol]";
    return($out);
  }
  else{
    my $out="$ggt$_[$gcol]\t$tgt$_[$tcol]";
    return($out);
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
