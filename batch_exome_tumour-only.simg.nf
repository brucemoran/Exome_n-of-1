#!/usr/bin/env nextflow

params.help = ""

if (params.help) {
  log.info ''
  log.info ' -----------------------------------------------------------'
  log.info '| BATCH TUMOUR-ONLY FASTQ QC, TRIM, ALIGN, SOMATIC SNV, CNA |'
  log.info ' -----------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info 'nextflow run batch_exome_tumour-only.simg.nf \
              --sampleCsv "sample.csv" \
              --refDir "refs"'
  log.info ''
  log.info 'This batch processing script only allows any number of tumour-only fastq pairs to be processed'
  log.info ''
  log.info 'Mandatory arguments:'
  log.info '    --sampleCsv      STRING      CSV format, headers: caseID,sampleID,read1,read2'
  log.info '    --refDir      STRING      dir in which reference data and required indices held; recommended to run associated reference creation NextFlow, DNAseq_references; still stuck on hg19 for several reasons;()'
  log.info ''
  exit 1
}

/* -2: Global Variables
*/
params.runDir = "$workflow.launchDir"
params.outDir = "$params.runDir/tumour-only_analysis"
params.scriptDir = "$params.runDir/scripts"
params.{params.exomebed} = "$params.refDir/exome.bed"

/* Reference data as params
*/
//For BWA use fasta as this is base prefix; also supply dict, fai
params.fasta = Channel.fromPath("$params.refDir/*fasta").getVal()
params.fai = Channel.fromPath("$params.refDir/*fasta.fai").getVal()
params.dict = Channel.fromPath("$params.refDir/*dict").getVal()

params.exomebedintlist = Channel.fromPath("$params.refDir/exome.bed.interval_list").getVal()
params.exomebed = Channel.fromPath("$params.refDir/exome.bed.gz").getVal()
params.exomebedtbi = Channel.fromPath("$params.refDir/exome.bed.gz.tbi").getVal()

params.dbsnp = Channel.fromPath("$params.refDir/dbsnp*.gz").getVal()
params.dbsnptbi = Channel.fromPath("$params.refDir/dbsnp*.tbi").getVal()

params.cosmic = Channel.fromPath("$params.refDir/COSMIC_CGC.bed").getVal()
params.ssrs = Channel.fromPath("$params.refDir/msisensor_microsatellites.list").getVal()
params.gps = Channel.fromPath("$params.refDir/mutect2_GetPileupSummaries.vcf.gz").getVal()
params.gpstbi = Channel.fromPath("$params.refDir/mutect2_GetPileupSummaries.vcf.gz.tbi").getVal()

/* -1: Install scripts required if not extant
*/
process scrpts {

  publishDir "$params.scriptDir", mode: "copy", pattern: "*"

  output:
  file('*') into completedmin1
  file('facets_cna.call.R') into facetscallscript
  file('facets_cna_consensus.call.R') into facetsconcscript
  file('facets_cna_consensus.func.R') into facetsconfscript
  file('filterMuTect2Format.pl') into filtermutect2script
  file('filterPiscesFormat.pl') into filterpiscesscript
  file('filterStrelka2IndelFormat.pl') into filterstrelka2iscript
  file('filterStrelka2SNVFormat.pl') into filterstrelka2sscript
  file('MuTect2_contamination.call.R') into mutect2contamscript
  file('variants_GRanges*.R') into variantsGRangesscript

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
completedmin1.subscribe { println "Scripts, images output to: $params.scriptDir" }

/* 0.000: Input using sample.csv
*/
Channel.fromPath("$params.sampleCsv", type: 'file')
       .splitCsv( header: true )
       .map { row -> [row.caseID, row.sampleID, file(row.read1), file(row.read2)] }
       .set { bbduking }

/* 0.0: Input trimming
*/
process bbduke {

  label 'c10_30G_cpu_mem'

  publishDir "$params.outDir/$caseID/$sampleID/bbduk", mode: "copy", pattern: "*.txt"

  input:
  set val(caseID), val(sampleID), file(read1), file(read2) from bbduking

  output:
  set val(caseID), val(sampleID), file(read1), file(read2) into fastpingpre
  set val(caseID), val(sampleID), file('*.bbduk.R1.fastq.gz'), file('*.bbduk.R2.fastq.gz') into (bwa_memming, fastpingpost)

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

/* 0.1: Input trimming
*/
process fastp {

  label 'c8_24G_cpu_mem'

  publishDir "$params.outDir/$caseID/$sampleID/fastp", mode: "copy", pattern: "*.html"

  input:
  set val(caseID), val(sampleID), file(read1), file(read2) from fastpingpre
  set val(caseID), val(sampleID), file(read1), file(read2) from fastpingpost

  output:
  file('*.json') into fastp_multiqc

  script:
  """
  fastp -w ${task.cpus} -h $sampleID".fastp.html" -j $sampleID".fastp.json" --in1 $read1 --in2 $read2
  """
}

/* 1.0: Input alignment
*/
process bwamem {

  label 'c20_60G_cpu_mem'

  publishDir "$params.outDir/$caseID/$sampleID/bwa", mode: "copy", pattern: "*[!bam]"

  input:
  set val(caseID), val(sampleID), file(read1), file(read2) from bwa_memming
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])

  output:
  set val(caseID), val(sampleID), file('*.bam'), file('*.bai') into dup_marking

  script:
  """
  DATE=\$(date +"%Y-%m-%dT%T")
  RGLINE="@RG\\tID:$sampleID\\tPL:ILLUMINA\\tSM:$sampleID\\tDS:tumour\\tCN:UCD\\tLB:LANE_X\\tDT:\$DATE"

  {
    bwa mem \
    -t ${task.cpus} \
    -M \
    -R \$RGLINE \
    ${params.fasta} \
    $read1 $read2 | \
    samtools sort -T "tmp."$sampleID -o $sampleID".sort.bam"

  samtools index $sampleID".sort.bam"

  samtools view -hC -T ${params.fasta} $sampleID".sort.bam" > $sampleID".sort.cram"
  } 2>&1 | tee > $sampleID".bwa-mem.log.txt"
  """
}

/* 1.1: MarkDuplicates
*/
process mrkdup {

  label 'c40_120G_cpu_mem'

  publishDir "$params.outDir/$caseID/$sampleID/picard/markdup", mode: "copy", pattern: "*.log.txt"
  publishDir "$params.outDir/$caseID/$sampleID/picard/metrics", mode: "copy", pattern: "*.metrics.txt"

  input:
  set val(caseID), val(sampleID), file(bam), file(bai) from dup_marking

  output:
  file('*md.metrics.txt') into mrkdup_multiqc
  set val(caseID), val(sampleID), file('*.md.bam'), file('*.md.bam.bai') into gatk4recaling

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

  label 'c10_30G_cpu_mem'

  publishDir "$params.outDir/$caseID/$sampleID/gatk4", mode: "copy", pattern: "*.txt"

  input:
  set val(caseID), val(sampleID), file(bam), file(bai) from gatk4recaling

  output:
  file('*.table') into gtkrcl_multiqc
  set val(caseID), val(sampleID), file('*.bqsr.bam'), file('*.bqsr.bam.bai') into (multimetricing, mutect2ing, facetsing, msisensoring, mantastrelka2ing, piscesing)

  script:
  """
  {
    gatk BaseRecalibrator \
    -R ${params.fasta} \
    -I $bam \
    --known-sites ${params.dbsnp} \
    --use-original-qualities \
    -O ${sampleID}.recal_data.table \
    -L ${params.{params.exomebed}intlist}

  #ApplyBQSR
  OUTBAM=\$(echo $bam | sed 's/bam/bqsr.bam/')
  gatk ApplyBQSR \
    -R ${params.fasta} \
    -I $bam \
    --bqsr-recal-file ${sampleID}.recal_data.table \
    --add-output-sam-program-record \
    --use-original-qualities \
    -O \$OUTBAM \
    -L ${params.{params.exomebed}intlist}

  samtools index \$OUTBAM
  } 2>&1 | tee > $sampleID".GATK4_recal.log.txt"
  """
}

/* 2.0: Metrics suite, this will produce an HTML report
*/
process mltmet {

  label 'c8_24G_cpu_mem'

  publishDir "$params.outDir/$caseID/$sampleID/metrics"

  input:
  set val(caseID), val(sampleID), file(bam), file(bai) from multimetricing

  output:
  file('*') into all2_0
  file('*.txt') into multimetrics_multiqc
  set file(fa), file(fai) into mutect2_fa

  script:
  """
  {
    picard-tools CollectHsMetrics \
      I=$bam \
      O=$sampleID".hs_metrics.txt" \
      TMP_DIR=./ \
      R=${params.fasta} \
      BAIT_INTERVALS=${params.{params.exomebed}intlist} \
      TARGET_INTERVALS=${params.{params.exomebed}intlist}

    picard-tools CollectAlignmentSummaryMetrics \
      I=$bam \
      O=$sampleID".AlignmentSummaryMetrics.txt" \
      TMP_DIR=./ \
      R=${params.fasta}

    picard-tools CollectMultipleMetrics \
      I=$bam \
      O=$sampleID".CollectMultipleMetrics.txt" \
      TMP_DIR=./ \
      R=${params.fasta}

    picard-tools CollectSequencingArtifactMetrics \
      I=$bam \
      O=$sampleID".artifact_metrics.txt" \
      TMP_DIR=./ \
      R=${params.fasta}

    picard-tools EstimateLibraryComplexity \
      I=$bam \
      O=$sampleID".est_lib_complex_metrics.txt" \
      TMP_DIR=./

    picard-tools CollectInsertSizeMetrics \
      I=$bam \
      O=$sampleID".insert_size_metrics.txt" \
      H=$bam".histogram.pdf" \
      TMP_DIR=./

  } 2>&1 | tee > $sampleID".picard-metrics.log.txt"
  """
}

/*2.1: SCNA with facets CSV snp-pileup
*/
process fctcsv {

  label 'c8_24G_cpu_mem'

  publishDir "$params.outDir/$caseID/$sampleID/facets"

  input:
  set val(caseID), val(sampleID), file(bam), file(bai) from facetsing
  file(facetsR) from facetscallscript

  output:
  file('*') into facetsoutputR

  script:
  """
  CSVFILE=\$(echo $bam | sed 's/bam/facets.r10.csv/')

  {
    snp-pileup \
      -r 10 \
      -p \
      ${params.dbsnp} \
      \$CSVFILE \
      $bam

    Rscript --vanilla ${params.fasta}cetsR \$CSVFILE
  } 2>&1 | tee > $sampleID".facets_snpp_call.log.txt"
  """
}

/* 2.2: MSIsensor
*/
process msisen {

  label 'c10_30G_cpu_mem'

  publishDir "$params.outDir/$sampleID/msisensor"
  publishDir "$params.outDir/calls/msisensor", pattern: '*.txt'

  input:
  set val(caseID), val(sampleID), file(tumourbam), file(tumourbai) from msisensoring

  output:
  val(sampleID) into completed2_2
  file('*') into msisensoroutput

  script:
  """
  msisensor msi \
    -d ${params.ssrs} \
    -t $tumourbam \
    -e ${params.exomebed} \
    -o $sampleID \
    -b ${task.cpus}

  MSI=\$( tail -n1 $sampleID | cut -f 3)
  mv $sampleID $sampleID".MSI-pc_"\$MSI".txt"
  """
}

/* 2.3: MuTect2
* NB --germline-resource dollar-sign{dbsnp} removed as no AF causing error
*/
process mutct2 {

  label 'c40_120G_cpu_mem'

  publishDir "$params.outDir/$sampleID/mutect2"
  publishDir "$params.outDir/calls/variants/vcf", pattern: '*raw.vcf'

  input:
  set val(caseID), val(sampleID), file(tumourbam), file(tumourbai) from mutect2ing
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
      --reference ${params.fasta} \
      --input $tumourbam \
      --tumor-sample $sampleID \
      --output $sampleID".md.recal.mutect2.vcf" \
      -L ${params.{params.exomebed}intlist}

    gatk --java-options ${params.full_javamem} \
      GetPileupSummaries \
      -I $tumourbam \
      -V ${params.gps} \
      -O $sampleID".getpileupsummaries.table" \
      -L ${params.exomebedintlist}

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
      -L ${params.{params.exomebed}intlist}

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

   label 'c10_30G_cpu_mem'

   publishDir "$params.outDir/", pattern: '*issue.table'

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

  label 'c40_120G_cpu_mem'

  publishDir "$params.outDir/$sampleID/manta-strelka2"
  publishDir "$params.outDir/calls/variants/vcf", pattern: '*raw.vcf'

  input:
  set val(caseID), val(sampleID), file(tumourbam), file(tumourbai) from mantastrelka2ing
  file(indelscript) from filterstrelka2iscript
  file(snvscript) from filterstrelka2sscript

  output:
  val(sampleID) into completed2_4
  file('*.pass.vcf') into strelka2_veping mode flatten
  file('*.raw.vcf') into strelka2_rawVcf
  file('manta/*') into completedmantacall

  script:
  """
  {
    configManta.py \
      --tumourBam=$tumourbam \
      --referenceFasta=${params.fasta} \
      --runDir=manta

    manta/runWorkflow.py -m local

    configureStrelkaSomaticWorkflow.py \
      --exome \
      --referenceFasta=${params.fasta} \
      --callRegions ${params.exomebed}gz \
      --indelCandidates=manta/results/variants/candidateSmallIndels.vcf.gz \
      --tumorBam=$tumourbam \
      --runDir=strelka2

    strelka2/runWorkflow.py -m local

    TUMOURSNVVCF=\$(echo $tumourbam | sed 's/bam/strelka2.snv.vcf/')
    gunzip -c strelka2/results/variants/somatic.snvs.vcf.gz | \
    perl -ane 'if(\$F[0]=~m/^#/){if(\$_=~m/^#CHROM/){
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

/* 2.5: Pisces
*/
process pisces {

  label 'c40_120G_cpu_mem'

  publishDir "$params.outDir/$sampleID/pisces"
  publishDir "$params.outDir/calls/variants/vcf", pattern: '*raw.vcf'

  input:
  set val(caseID), val(sampleID), file(tumourbam), file(tumourbai) from piscesing
  file(filterpl) from filterpiscesscript

  output:
  file('*.pass.vcf') into pisces_veping mode flatten
  file('*.raw.vcf') into pisces_rawVcf
  file('*') into completedpiscescall

  script:
  """
  TUMOURVCF=\$(echo $tumourbam | sed 's/bam/vcf/')
  {
    dotnet /app/CreateGenomeSizeFile_5.2.9.122/CreateGenomeSizeFile.dll -g ./ -s "Homo sapiens (hg19)"

    dotnet /app/Pisces_5.2.9.122/Pisces.dll \
      --bam $tumourbam \
      --genomepaths ./ \
      --intervalpaths ${params.exomebedintlist} \
      --filterduplicates true \
      --onlyuseproperpairs true \
      --mindepth 3 \
      --gvcf false \
      --outfolder ./ \
      --maxthreads ${task.cpus}

    perl $filterpl \
      ID=$sampleID \
      DP=14
      MD=2 \
      VCF=\$TUMOURVCF

  } 2>&1 | tee > $sampleID".pisces.log.txt"

  """
}

/* 3.0: Annotate Vcfs
*/
ALLVCFS = pisces_veping
          .mix( mutect2_veping )
          .mix( strelka2_veping )

process vepann {

  label 'c10_30G_cpu_mem'

  publishDir "$params.outDir/calls/variants/vcf", pattern: '*.vcf'

  input:
  each file(vcf) from ALLVCFS

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
    --fasta ${params.fasta} \
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
             .mix(pisces_rawVcf)
             .mix(mutect2_rawVcf)
             .mix(strelka2_rawVcf)

vartypes = Channel.from( "snv", "indel" )

process vcfGRa {

  label 'c20_60G_cpu_mem'

  publishDir "$params.outDir/calls/variants/pdf", pattern: '*.pdf'
  publishDir "$params.outDir/calls/variants/vcf", pattern: '*.vcf'
  publishDir "$params.outDir/calls/variants/data", pattern: '*.{RData,tab}'

  input:
  file(rawGRangesvcff) from ALLRAWVEPVCFS.collect()
  each vartype from vartypes
  set file(callR), file(funcR) from variantsGRangesscript

  output:
  val(vartype) into completed3_1
  file('*') into completedvcfGRangesConsensus

  script:
  """
  OUTID=\$(basename ${params.runDir})
  Rscript --vanilla $callR \
    $funcR \
    "NULL" \
    $vartype".pass.vep.vcf" \
    \$OUTID \
    ${params.includeOrder}
  """
}

/* 4.0 Run multiQC to finalise report
*/
process mltiQC {

  label 'c10_30G_cpu_mem'

  publishDir "$params.outDir", mode: "copy", pattern: "*.html"

  input:
  file('fastp/*') from fastp_multiqc.collect()
  file('mrkdup/*') from mrkdup_multiqc.collect()
  file('gtkrcl/*') from gtkrcl_multiqc.collect()
  file('multimetrics/*') from multimetrics_multiqc.collect()

  output:
  file('*') into completedmultiqc

  script:
  """
  OUTID=\$(basename ${params.runDir})
  multiqc ./ -i \$OUTID --tag DNA -f -c /usr/local/multiqc_config_BMB.yaml
  """
}
