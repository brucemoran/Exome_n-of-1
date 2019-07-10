#!/usr/bin/env nextflow

params.help = ""

if (params.help) {
  log.info ''
  log.info '----------------------------------------------'
  log.info 'NEXTFLOW BAM QC, TRIM, ALIGN, SOMATIC SNV, CNA'
  log.info '----------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info 'nextflow run exome_n-of-1.simg.nf \
              --sampleCsv "data/sample.csv" \
              --refDir "refs" \
              --includeOrder "tumour_A,tumour_B"'
  log.info ''
  log.info 'Mandatory arguments:'
  log.info '    --sampleCsv      STRING      CSV format, headers: type (either "germline" or "somatic"),sampleID,/path/to/read1.fastq.gz,/path/to/read2.fastq.gz '
  log.info '    --refDir      STRING      dir in which reference data and required indices held; recommended to run associated reference creation NextFlow, DNAseq_references; stuck on GRCh37 for several reasons=('
  log.info '    --sampleMeta      STRING      CSV format, headers: sampleID,meta; use meta for naming in PCGR, CPSR report; sampleID matches from sampleCsv [[[replace with meta entry for sampleCsv]]]'
  log.info ''
  log.info 'Optional arguments:'
  log.info '    --includeOrder      STRING      in final plots, use this ordering of samples (if multiple somatic samples); comma-separated, no spaces'
  log.info '    --germline      STRING      run HaplotypeCaller on germline sample and annotate with CPSR'
  log.info ''
  exit 1
}

/* -2: Global Variables
*/
params.outDir = "$workflow.launchDir/analysis"
params.scriptDir = "$workflow.launchDir/scripts"

/* -1: Install scripts required if not extant
*/
process scrpts {

  publishDir path: "$params.scriptDir", mode: "copy", pattern: "*"

  output:
  file('facets_cna.call.R') into facetscallscript
  file('facets_cna_consensus.call.R') into facetsconcscript
  file('facets_cna_consensus.func.R') into facetsconfscript
  file('filterLancetSomaticFormat.pl') into filterlancetscript
  file('filterMuTect2Format.pl') into filtermutect2script
  file('filterStrelka2IndelFormat.pl') into filterstrelka2iscript
  file('filterStrelka2SNVFormat.pl') into filterstrelka2sscript
  file('MuTect2_contamination.call.R') into mutect2contamscript
  file('SomaticSignatures.call.R') into somaticsignaturesscript
  file('variants_GRanges*.R') into variantsGRangesscript
  file('vcf42.head.txt') into (pcgrVcfHead, completedscrpts)

  script:
  """
  git clone https://github.com/brucemoran/Exome_n-of-1
  mv ./Exome_n-of-1/scripts/* ./
  rm -rf ./Exome_n-of-1

  git clone https://github.com/brucemoran/somaticVariantConsensus
  mv ./somaticVariantConsensus/scripts/* ./
  rm -rf ./somaticVariantConsensus
  """
}
completedscrpts.subscribe{ println "Completed scripts" }

/* 0.00: Input using sample.csv
*/
Channel.fromPath("$params.sampleCsv", type: 'file')
       .splitCsv( header: true )
       .map { row -> [row.type, row.sampleID, file(row.read1), file(row.read2)] }
       .set { bbduking }

/* 0.0: Input trimming
*/
process bbduk {

  publishDir path: "$params.outDir/$sampleID/bbduk", mode: "copy", pattern: "*.txt"

  input:
  set val(type), val(sampleID), file(read1), file(read2) from bbduking

  output:
  file('*') into completed0_0
  set val(type), val(sampleID), file('*.bbduk.R1.fastq.gz'), file('*.bbduk.R2.fastq.gz') into bwa_memming
  set val(type), val(sampleID), file('*.bbduk.R1.fastq.gz'), file('*.bbduk.R2.fastq.gz'), file(read1), file(read2) into fastping
  set val(type), val(sampleID), file(read1), file(read2) into fastqcing

  script:
  """
  {
    sh bbduk.sh ${params.quarter_javamem} \
      in1=$read1 \
      in2=$read2 \
      out1=$sampleID".bbduk.R1.fastq.gz" \
      out2=$sampleID".bbduk.R2.fastq.gz" \
      k=31 \
      mink=5 \
      hdist=1 \
      ktrim=r \
      trimq=20 \
      qtrim=rl \
      maq=20 \
      ref=/usr/local/bbmap/resources/adapters.fa \
      tpe \
      tbo \
      stats=$sampleID".bbduk.adapterstats.txt" \
      overwrite=T
  } 2>&1 | tee > $sampleID".bbduk.runstats.txt"
  """
}

/* 0.1: QC of per, post-bbduk
*/
process fastp {

  publishDir "$params.outDir/$sampleID/fastp", mode: "copy", pattern: "*.html"

  input:
  set val(type), val(sampleID), file(preread1), file(preread2), file(postread1), file(postread2) from fastping

  output:
  file('*.html') into completed0_1
  file('*.json') into fastp_multiqc

  script:
  """
  fastp -w ${task.cpus} -h $sampleID"_pre.fastp.html" -j $sampleID"_pre.fastp.json" --in1 $preread1 --in2 $preread2

  fastp -w ${task.cpus} -h $sampleID"_post.fastp.html" -j $sampleID"_post.fastp.json" --in1 $postread1 --in2 $postread2
  """
}

/* 0.1: QC of per, post-bbduk
*/
process fastqc {

  publishDir "$params.outDir/$sampleID/fastqc", mode: "copy", pattern: "*.html"

  input:
  set val(type), val(sampleID), file(read1), file(read2) from fastqcing

  output:
  file('*.html') into fastqc_multiqc

  script:
  """
  #!/bin/bash
  fastqc -t ${task.cpus} --noextract -o ./ $read1 $read2
  """
}

/* 1.0: Input alignment
*/
process bwamem {

  cache 'deep'

  publishDir path: "$params.outDir/$sampleID/bwa", mode: "copy", pattern: "*.log.txt"

  input:
  set val(type), val(sampleID), file(read1), file(read2) from bwa_memming
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  set file(amb), file(ann), file(bwt), file(pac), file(sa) from Channel.value([params.amb, params.ann, params.bwt, params.pac, params.sa])

  output:
  set val(type), val(sampleID), file('*.bam'), file('*.bai') into (cramming, dup_marking)
  file('*.log.txt') into completedbwa

  script:
  """
  DATE=\$(date +"%Y-%m-%dT%T")
  RGLINE="@RG\\tID:$sampleID\\tPL:ILLUMINA\\tSM:$sampleID\\tDS:$type\\tCN:UCD\\tLB:LANE_X\\tDT:\$DATE"

  {
  bwa mem \
    -t${task.cpus} \
    -M \
    -R \$RGLINE \
    $fasta \
    $read1 $read2 | \
    samtools sort -T "tmp."$sampleID -o $sampleID".sort.bam"

  samtools index $sampleID".sort.bam"

  } 2>&1 | tee > $sampleID".bwa-mem.log.txt"
  """
}

/* 1.01: CRAM alignment and output
* WIP: output iupload shema for ENA/EGA
*/
process cram {

  publishDir path: "$params.outDir/$sampleID/bwa", mode: "copy", pattern: "*.cra*"

  input:
  set val(type), val(sampleID), file(bam), file(bai) from cramming
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])

  output:
  set file('*.cram'), file('*.crai') into completedcram

  script:
  """
  samtools view -hC -T $fasta $sampleID".sort.bam" > $sampleID".sort.cram"
  samtools index $sampleID".sort.cram"
  """
}

/* 1.1: MarkDuplicates
*/
process mrkdup {

  publishDir path: "$params.outDir/$sampleID/picard/markdup", mode: "copy", pattern: "*[!.metrics.txt]"
  publishDir path: "$params.outDir/$sampleID/picard/metrics", mode: "copy", pattern: "*.metrics.txt"

  input:
  set val(type), val(sampleID), file(bam), file(bai) from dup_marking

  output:
  file('*md.metrics.txt') into mrkdup_multiqc
  set val(type), val(sampleID), file('*.md.bam'), file('*.md.bam.bai') into gatk4recaling

  script:
  """
  OUTBAM=\$(echo $bam | sed 's/bam/md.bam/')
  OUTMET=\$(echo $bam | sed 's/bam/md.metrics.txt/')
  {
    picard-tools ${params.quarter_javamem} \
      MarkDuplicates \
      TMP_DIR=./ \
      INPUT=$bam \
      OUTPUT=/dev/stdout \
      COMPRESSION_LEVEL=0 \
      QUIET=TRUE \
      METRICS_FILE=\$OUTMET \
      REMOVE_DUPLICATES=FALSE \
      ASSUME_SORTED=TRUE \
      VALIDATION_STRINGENCY=LENIENT \
      VERBOSITY=ERROR | samtools view -Shb - > \$OUTBAM

  samtools index \$OUTBAM
  } 2>&1 | tee > $sampleID".picard-tools_markDuplicates.log.txt"
  """
}

/* 1.2: GATK4 BestPractices
* as per best practices of GATK4
*/
process gtkrcl {

  label 'c20_60G_cpu_mem'

  publishDir path: "$params.outDir/$sampleID/gatk4/bestpractice", mode: "copy"

  input:
  set val(type), val(sampleID), file(bam), file(bai) from gatk4recaling
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  set file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])
  file(exomeintlist) from Channel.value(params.exomeintlist)

  output:
  file('*.table') into gtkrcl_multiqc
  set val(type), val(sampleID), file('*.bqsr.bam') into (germfiltering, somafiltering)
  set val(type), val(sampleID), file('*.bqsr.bam') into gatk_germ

  script:
  """
  {
    gatk BaseRecalibrator \
    -R $fasta \
    -I $bam \
    --known-sites $dbsnp \
    --use-original-qualities \
    -O ${sampleID}.recal_data.table \
    -L $exomeintlist

  #ApplyBQSR
  OUTBAM=\$(echo $bam | sed 's/bam/bqsr.bam/')
  gatk ApplyBQSR \
    -R $fasta \
    -I $bam \
    --bqsr-recal-file ${sampleID}.recal_data.table \
    --add-output-sam-program-record \
    --use-original-qualities \
    -O \$OUTBAM \
    -L $exomeintlist

  } 2>&1 | tee > $sampleID".GATK4_recal.log.txt"
  """
}

/* 1.21: GATK4 Germline
*/
process gatkgerm {

  label 'c40_120G_cpu_mem'

  publishDir "$params.outDir/$type/$sampleID/gatk4/HC_germline", mode: "copy", pattern: "*"

  input:
  set val(type), val(sampleID), file(bam) from gatk_germ
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  set file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])
  file(exomeintlist) from Channel.value(params.exomeintlist)
  set file(omni), file(otbi), file(kgp1), file(ktbi), file(hpmp), file(htbi) from Channel.value([params.omni, params.otbi, params.kgp1, params.ktbi, params.hpmp, params.htbi])

  output:
  set val(sampleID), file('*.HC.vcf.gz'), file('*.HC.vcf.gz.tbi') into germ_vcf

  script:
  """
  {
    #HaplotypeCaller
    INPUTBAM=$bam
    OUTVCF=\$(echo \$INPUTBAM | sed 's/bam/hc.vcf/')
    gatk --java-options ${params.full_javamem} HaplotypeCaller \
      -R $fasta \
      -I \$INPUTBAM \
      --dont-use-soft-clipped-bases \
      --standard-min-confidence-threshold-for-calling 20 \
      --dbsnp $dbsnp \
      --native-pair-hmm-threads ${task.cpus} \
      -O $sampleID".HC.vcf" \
      -L $exomeintlist

    bgzip $sampleID".HC.vcf"
    tabix $sampleID".HC.vcf.gz"

  } 2>&1 | tee $sampleID".GATK4_HaplotypeCaller-germline.log.txt"

  """
}

/* 1.22: CPSR annotation of GATK4 Germline
*/
process cpsrreport {

  label 'c10_30G_cpu_mem'

  publishDir "$params.outDir/calls/reports", mode: "copy", pattern: "*.html"
  publishDir "$params.outDir/calls/variants/cpsr", mode: "copy", pattern: "*[!.html]"

  input:
  set val(sampleID), file(vcf), file(tbi) from germ_vcf
  each file(cpsr_grch37) from CPSR
  each file(cpsrmeta) from Channel.value([params.cpsrpcgr])

  output:
  file('*') into cpsr_vcfs

  script:
  """
  {
    METAID=\$(grep $sampleID $cpsrmeta | cut -d "," -f 2)
    cpsr.py \
      --no-docker \
      --no_vcf_validate \
      $vcf \
      $cpsr_grch37 \
      ./ \
      grch37 \
      0 \
      $cpsr_grch37"/data/grch37/cpsr_configuration_default.toml" \
      \$METAID

  } 2>&1 | tee $sampleID".cpsr.log.txt"
  """
}

/* 1.3: filter germ channel for all processes subsequent, index bam
*/
process grmflt {

  input:
  set val(type), val(sampleID), file(bam) from germfiltering

  output:
  file(bam) into germbamcombine
  file('*.bam.bai') into germbaicombine
  set val(sampleID), file(bam), file('*.bam.bai') into gmultimetricing

  when:
  type == "germline"

  script:
  """
  samtools index $bam > $bam".bai"
  """
}

/* 1.4: filter somatic channel for all processes subsequent, index bam
*/
process smaflt {

  input:
  set val(type), val(sampleID), file(bam) from somafiltering

  output:
  set val(sampleID), file(bam), file ('*.bam.bai') into (multimetricing, germcombine)

  when:
  type != "germline"

  script:
  """
  samtools index $bam > $bam".bai"
  """
}

/*1.5 combine germline with somatic and unique those outputs
*/
process combinegs {

    input:
    set val(sampleID), file(bam), file (bai) from germcombine
    each file(germlinebam) from germbamcombine
    each file(germlinebai) from germbaicombine

    output:
    set val(sampleID), file(bam), file(bai), stdout, file(germlinebam), file(germlinebai) into ( mutect2somaticing, facetsomaing, msisensoring, mantastrelka2ing, lanceting )
    stdout into vcfGRaID

    """
    echo $germlinebam | perl -ane '@s=split(/\\./); print \$s[0];'
    """
}

/* 2.0: PicardTools metrics suite for MultiQC HTML report
*/
MULTIALL = gmultimetricing.mix(multimetricing)

process mltmet {

  label 'c8_24G_cpu_mem'

  publishDir "$params.outDir/$sampleID/metrics"

  input:
  set val(sampleID), file(bam), file(bai) from MULTIALL
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  file(exomeintlist) from Channel.value(params.exomeintlist)

  output:
  file('*') into all2_0
  file('*.txt') into multimetrics_multiqc

  script:
  """
  {
    picard-tools CollectHsMetrics \
      I=$bam \
      O=$sampleID".hs_metrics.txt" \
      TMP_DIR=./ \
      R=$fasta \
      BAIT_INTERVALS=$exomeintlist  \
      TARGET_INTERVALS=$exomeintlist

    picard-tools CollectAlignmentSummaryMetrics \
      I=$bam \
      O=$sampleID".AlignmentSummaryMetrics.txt" \
      TMP_DIR=./ \
      R=$fasta

    picard-tools CollectMultipleMetrics \
      I=$bam \
      O=$sampleID".CollectMultipleMetrics.txt" \
      TMP_DIR=./ \
      R=$fasta

    picard-tools CollectSequencingArtifactMetrics \
      I=$bam \
      O=$sampleID".artifact_metrics.txt" \
      TMP_DIR=./ \
      R=$fasta

    picard-tools EstimateLibraryComplexity \
      I=$bam \
      O=$sampleID".est_lib_complex_metrics.txt" \
      TMP_DIR=./

    picard-tools CollectInsertSizeMetrics \
      I=$bam \
      O=$sampleID".insert_size_metrics.txt" \
      H=$bam".histogram.pdf" \
      TMP_DIR=./

  } 2>&1 | tee > $sampleID".picard.metrics.log"
  """
}

/*2.1: SCNA with facets CSV snp-pileup
*/
process fctcsv {

  label 'c8_24G_cpu_mem'

  publishDir "$params.outDir/$sampleID/facets"

  input:
  set val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from facetsomaing
  set file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])
  file(facetsR) from facetscallscript

  output:
  file('*.tab') into facets_consensusing
  file('*.pcgr.tsv') into facets_pcgr
  file('*') into facetsoutputR

  script:
  """
  CSVFILE=\$(echo $tumourbam | sed 's/bam/facets.r10.csv/')

  {
    snp-pileup \
      $dbsnp \
      -r 10 \
      -p \
      \$CSVFILE \
      $germlinebam \
      $tumourbam

    Rscript --vanilla $facetsR \$CSVFILE

    echo -e "Chromosome\\tStart\\tEnd\\tSegment_Mean" > $sampleID".cncf-jointsegs.pcgr.tsv"
    tail -n+2 $sampleID".fit_cncf-jointsegs.tab" | awk '{print \$1"\\t"\$10"\\t"\$11"\\t"\$5}' >> $sampleID".cncf-jointsegs.pcgr.tsv"

    tail -n+2 $sampleID".fit_ploidy-purity.tab" > $sampleID".fit_ploidy-purity.pcgr.tsv"

  } 2>&1 | tee > $sampleID".facets_snpp_call.log.txt"
  """
}

/* 2.11: SCNA consensus from facets
*/
process fctcon {

  label 'c8_24G_cpu_mem'

  publishDir "$params.outDir/calls/scna/facets"

  input:
  file(filesn) from facets_consensusing.collect()
  file(callR) from facetsconcscript
  file(funcR) from facetsconfscript
  file(dict) from Channel.value(params.dict)
  file(cosmic) from Channel.value(params.cosmic)

  output:
  file('*') into completed2_11

  script:
  """
  {
  OUTID=\$(basename ${workflow.launchDir})
  Rscript --vanilla $callR \
    $dict \
    $cosmic \
    \$OUTID \
    $funcR
  } 2>&1 | tee > "facets_cons.log.txt"
  """
}

/* 2.2: MuTect2
* NB --germline-resource dollar-sign{dbsnp} removed as no AF causing error
*/
process mutct2 {

  label 'c40_120G_cpu_mem'

  publishDir path: "$params.outDir/$sampleID/mutect2", mode: "copy"
  publishDir path: "$params.outDir/calls/variants/vcf", mode: "copy", pattern: '*raw.vcf'

  input:
  set val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from mutect2somaticing
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  file(exomeintlist) from Channel.value(params.exomeintlist)
  set file(gps), file(gpstbi) from Channel.value([params.gps, params.gpstbi])
  file(filterpl) from filtermutect2script

  output:
  file('*.pass.vcf') into mutect2_veping mode flatten
  file('*.raw.vcf') into mutect2_rawVcf
  file('*') into completedmutect2call
  set val(sampleID), file('*calculatecontamination.table') into contamination

  script:
  """
  {
    gatk --java-options ${params.full_javamem} \
      Mutect2 \
      --native-pair-hmm-threads ${task.cpus} \
      --reference $fasta \
      --input $germlinebam \
      --input $tumourbam \
      --normal-sample $germlineID \
      --tumor-sample $sampleID \
      --output $sampleID".md.recal.mutect2.vcf" \
      -L $exomeintlist

    gatk --java-options ${params.full_javamem} \
      GetPileupSummaries \
      -I $tumourbam \
      -V $gps \
      -O $sampleID".getpileupsummaries.table" \
      -L $exomeintlist

    gatk CalculateContamination \
      -I $sampleID".getpileupsummaries.table" \
      -O $sampleID".calculatecontamination.table"

    gatk --java-options ${params.full_javamem} \
      FilterMutectCalls \
      --contamination-table $sampleID".calculatecontamination.table" \
      --interval-padding 5 \
      --output $sampleID".md.recal.mutect2.FilterMutectCalls.vcf" \
      --unique-alt-read-count 3 \
      --variant $sampleID".md.recal.mutect2.vcf" \
      -L $exomeintlist

    perl $filterpl \
      ID=$sampleID \
      DP=14 \
      MD=2 \
      VCF=$sampleID".md.recal.mutect2.FilterMutectCalls.vcf"

  } 2>&1 | tee > $sampleID".GATK4_mutect2.log.txt"
  """

}

/* 2.31: MuTect2 Contamination
*/
process mutct2_contam {

  echo true

  label 'c40_120G_cpu_mem'

  publishDir path: "$params.outDir/", mode: "copy", pattern: '*issue.table'

  input:
  set val(sampleID), file(contable) from contamination
  file script from mutect2contamscript

  output:
  file('*.table') into completedcontam

  """
  Rscript --vanilla $script $contable $sampleID
  """
}

/* 2.4: Manta output is a pre-req for Strelka2, so call both here
*/
process mntstr {

  label 'c20_60G_cpu_mem'

  publishDir path: "$params.outDir/$sampleID/manta-strelka2"
  publishDir path: "$params.outDir/calls/variants/vcf", mode: "copy", pattern: '*raw.vcf'

  input:
  set val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from mantastrelka2ing
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  set file(exomebedgz), file(exomebedgztbi) from Channel.value([params.exomebedgz, params.exomebedgztbi])
  file(indelscript) from filterstrelka2iscript
  file(snvscript) from filterstrelka2sscript

  output:
  val(sampleID) into completed2_4
  file('*.pass.vcf') into strelka2_merge
  file('*.raw.vcf') into strelka2_rawVcf
  file('manta/*') into completedmantacall

  script:
  """
  {
    configManta.py \
      --normalBam=$germlinebam \
      --tumourBam=$tumourbam \
      --referenceFasta=$fasta \
      --runDir=manta

    manta/runWorkflow.py -m local

    configureStrelkaSomaticWorkflow.py \
      --exome \
      --referenceFasta=$fasta \
      --callRegions $exomebedgz \
      --indelCandidates=manta/results/variants/candidateSmallIndels.vcf.gz \
      --normalBam=$germlinebam \
      --tumorBam=$tumourbam \
      --runDir=strelka2

    strelka2/runWorkflow.py -m local

    TUMOURSNVVCF=\$(echo $tumourbam | sed 's/bam/strelka2.snv.vcf/')
    gunzip -c strelka2/results/variants/somatic.snvs.vcf.gz | \
    perl -ane 'if(\$F[0]=~m/^#/){if(\$_=~m/^#CHROM/){
        \$_=~s/NORMAL/$germlineID/;
        \$_=~s/TUMOR/$sampleID/;
        print \$_;next;}
        else{print \$_;next;}
      }
      else{print \$_;}' > \$TUMOURSNVVCF

    perl $snvscript \
     ID=$sampleID \
     DP=14 \
     MD=2 \
     VCF=\$TUMOURSNVVCF

    TUMOURINDELVCF=\$(echo $tumourbam | sed 's/bam/strelka2.indel.vcf/')
    gunzip -c strelka2/results/variants/somatic.indels.vcf.gz | \
    perl -ane 'if(\$F[0]=~m/^#/){if(\$_=~m/^#CHROM/){
        \$_=~s/NORMAL/$germlineID/;
        \$_=~s/TUMOR/$sampleID/;
        print \$_;next;}
        else{print \$_;next;}}
      else{print \$_;}' > \$TUMOURINDELVCF

    perl $indelscript \
      ID=$sampleID \
      DP=14 \
      MD=2 \
      VCF=\$TUMOURINDELVCF

  } 2>&1 | tee > $sampleID".manta-strelka2.log.txt"
  """
}

/* 2.41: MergeVcf passes from snv indel from Manta/Strelka2
*/
process mntstrmerge {

  label 'c5_15G_cpu_mem'

  input:
  file(vcfs) from strelka2_merge.collect()

  output:
  file('*.snv_indel.pass.vcf') into strelka2_veping

  script:
  """
    gatk MergeVcfs \
      -I 12922_18-A11.strelka2.snv.pass.vcf \
      -I 12922_18-A11.strelka2.indel.pass.vcf \
      -O 12922_18-A11.strelka2.snv_indel.pass.vcf
  """
}

/* 2.5: Lancet
*/
process lancet {

  process.container = "library://bruce.moran/default/variant_callers:lancet.centos7"
  singularity.enabled = true

  publishDir path: "$params.outDir/$sampleID/lancet"
  publishDir path: "$params.outDir/calls/variants/vcf", mode: "copy", pattern: '*raw.vcf'

  input:
  set val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from lanceting
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  file(exomebed) from Channel.value(params.exomebed)
  file(filterLancet) from filterlancetscript

  output:
  file('*.pass.vcf') into lancet_veping mode flatten
  file('*.raw.vcf') into lancet_rawVcf
  file('*') into completedlancetcall

  script:
  """
  TUMOURVCF=\$(echo $tumourbam | sed 's/bam/lancet.vcf/')
  {
    lancet \
      --num-threads ${task.cpus} \
      --ref $fasta \
      --bed $exomebed \
      --tumor $tumourbam \
      --normal $germlinebam | \
      perl -ane 'if(\$F[0]=~m/^\\#CHROM/){
        \$_=~s/TUMOR/$sampleID/;
        \$_=~s/NORMAL/$germlineID/;
        print \$_;}
      else{print \$_;}' > \$TUMOURVCF

    perl $filterLancet \
      ID=$sampleID \
      DP=14 \
      MD=2 \
      VCF=\$TUMOURVCF

  } 2>&1 | tee > $sampleID".lancet.log.txt"

  """
}

/* 3.0: Annotate Vcfs
*/
ALLVCFS = lancet_veping
          .mix( mutect2_veping )
          .mix( strelka2_veping )

process vepann {

  label 'c10_30G_cpu_mem'

  publishDir path: "$params.outDir/calls/variants/vcf", mode: "copy", pattern: '*.vcf'

  input:
  each file(vcf) from ALLVCFS
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])

  output:
  file('*.vcf') into annoVcfs
  file('*.vcf') into (completedvep, runGRanges)

  script:
  """
  VCFANNO=\$(echo $vcf | sed "s/.vcf/.vep.vcf/")

  vep --dir_cache /usr/local/ensembl-vep/cache \
    --offline \
    --assembly GRCh37 \
    --vcf_info_field ANN \
    --symbol \
    --species homo_sapiens \
    --check_existing \
    --cache \
    --merged \
    --fork ${task.cpus} \
    --af_1kg \
    --af_gnomad \
    --vcf \
    --input_file $vcf \
    --output_file \$VCFANNO \
    --format "vcf" \
    --fasta $fasta \
    --hgvs \
    --canonical \
    --ccds \
    --force_overwrite \
    --verbose
  """
}

/* 3.1 RData GRanges from processed VCFs
* take publishDir and check for number of files therein
* each sample has 9 associated (raw,snv,indel per caller)
* NB increment if adding callers!
*/
ALLRAWVEPVCFS = runGRanges
             .mix(lancet_rawVcf)
             .mix(mutect2_rawVcf)
             .mix(strelka2_rawVcf)

process vcfGRa {

  label 'c20_60G_cpu_mem'

  publishDir "$params.outDir/calls/variants/pdf", pattern: '*.pdf'
  publishDir "$params.outDir/calls/variants/vcf", pattern: '*.vcf'
  publishDir "$params.outDir/calls/variants/data", pattern: '*.[*RData,*tab]'

  input:
  file(rawGRangesvcff) from ALLRAWVEPVCFS.collect()
  set file(callR), file(funcR) from variantsGRangesscript
  file(vcfHead) from pcgrVcfHead
  val(vcfGermlineID) from vcfGRaID

  output:
  file('*.ALL.pcgr.all.tab.vcf') into pcgrvcfs
  file('*') into completedvcfGRangesConsensus

  script:
  """
  OUTID=\$(basename ${workflow.launchDir})
  Rscript --vanilla $callR \
    $funcR \
    $vcfGermlineID \
    "snv_indel.pass.vep.vcf" \
    \$OUTID \
    ${params.includeOrder}

  ##header VCF
  for VCF in *.pcgr.all.tab.vcf; do
    cat $vcfHead > 1;
    cat \$VCF >> 1;
    mv 1 \$VCF;
  done
  """
}


/* 3.2 PCGR report
* take all mutations in consensus.tab from pass.vcfs into single VCF for PCGR
*/
process pcgrreport {

  label 'c10_30G_cpu_mem'

  publishDir "$params.outDir/calls/reports", mode: "copy", pattern: "*html"
  publishDir "$params.outDir/calls/variants/pcgr", mode: "copy", pattern: "*[!.html]"

  input:
  file(vcf) from pcgrvcfs
  file(facetstsvs) from facets_pcgr.collect()
  file(pcgr_grch37) from PCGR
  file(pcgrmeta) from Channel.value([params.cpsrpcgr])

  output:
  file('*') into completedPCGR

  script:
  """
  {
  for VCF in *vcf; do
    SAMPLEID=\$(echo \$VCF | cut -d "." -f 1)
    METAID=\$(grep \$SAMPLEID $pcgrmeta | cut -d "," -f 2)
    PLOIDY=\$(cut -f 1 \$SAMPLEID".fit_ploidy-purity.pcgr.tsv")
    PURITY=\$(cut -f 2 \$SAMPLEID".fit_ploidy-purity.pcgr.tsv")
    pcgr.py $pcgr_grch37 \
      ./ \
      grch37 \
      $pcgr_grch37/data/grch37/pcgr_configuration_default.toml \
      \$METAID \
      --input_vcf \$VCF \
      --input_cna \$SAMPLEID".cncf-jointsegs.pcgr.tsv" \
      --tumor_ploidy \$PLOIDY \
      --tumour_purity \$PURITY \
      --no-docker \
      --force_overwrite \
      --no_vcf_validate
  done
  } 2>&1 | tee pcgr.log.txt
  """
}

/* 4.0 Run multiQC to finalise report
*/
process mltiQC {

  label 'c10_30G_cpu_mem'

  publishDir path: "$params.outDir/calls/reports", mode: "copy", pattern: "*html"

  input:
  file(fastps) from fastp_multiqc.collect()
  file(fastqcs) from fastqc_multiqc.collect()
  file(gtkrcls) from gtkrcl_multiqc.collect()
  file(multimetrics) from multimetrics_multiqc.collect()
  file(mrkdups) from mrkdup_multiqc.collect()

  output:
  file('*') into completedmultiqc

  script:
  """
  OUTID=\$(basename ${workflow.launchDir})
  multiqc . -i \$OUTID --tag DNA -f -c /usr/local/multiqc_config_BMB.yaml
  """
}
