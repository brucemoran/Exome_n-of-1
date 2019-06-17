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
  log.info '    --sampleCsv      STRING      CSV format, headers: caseID,sampleID,meta,read1,read2'
  log.info '    --refDir      STRING      dir in which reference data and required indices held; recommended to run associated reference creation NextFlow, DNAseq_references; still stuck on hg19 for several reasons;()'
  log.info ''
  exit 1
}

/* -2: Global Variables
*/
params.runDir = "$workflow.launchDir"
params.outDir = "$params.runDir/tumour-only_analysis"
params.scriptDir = "$params.runDir/scripts"

/* Reference data as params
*/
params.fasta = Channel.fromPath("$params.refDir/*fasta").getVal()
params.fai = Channel.fromPath("$params.refDir/*fasta.fai").getVal()
params.dict = Channel.fromPath("$params.refDir/*dict").getVal()

params.amb = Channel.fromPath("$params.refDir/*fasta.amb").getVal()
params.ann = Channel.fromPath("$params.refDir/*fasta.ann").getVal()
params.bwt = Channel.fromPath("$params.refDir/*fasta.bwt").getVal()
params.pac = Channel.fromPath("$params.refDir/*fasta.pac").getVal()
params.sa = Channel.fromPath("$params.refDir/*fasta.sa").getVal()

params.exomeintlist = Channel.fromPath("$params.refDir/exome.bed.interval_list").getVal()
params.exomebed = Channel.fromPath("$params.refDir/exome.bed").getVal()
params.exomebedgz = Channel.fromPath("$params.refDir/exome.bed.gz").getVal()
params.exomebedgztbi = Channel.fromPath("$params.refDir/exome.bed.gz.tbi").getVal()

params.dbsnp = Channel.fromPath("$params.refDir/dbsnp*.gz").getVal()
params.dbsnptbi = Channel.fromPath("$params.refDir/dbsnp*.tbi").getVal()
params.omni = Channel.fromPath("$params.refDir/KG_omni*.gz").getVal()
params.otbi = Channel.fromPath("$params.refDir/KG_omni*.gz.tbi").getVal()
params.kgp1 = Channel.fromPath("$params.refDir/KG_phase1*.gz").getVal()
params.ktbi = Channel.fromPath("$params.refDir/KG_phase1*.gz.tbi").getVal()
params.hpmp = Channel.fromPath("$params.refDir/hapmap*.gz").getVal()
params.htbi = Channel.fromPath("$params.refDir/hapmap*.gz.tbi").getVal()
params.mlls = Channel.fromPath("$params.refDir/Mills*.gz").getVal()
params.mtbi = Channel.fromPath("$params.refDir/Mills*.gz.tbi").getVal()

params.cosmic = Channel.fromPath("$params.refDir/COSMIC_CGC.bed").getVal()
params.ssrs = Channel.fromPath("$params.refDir/msisensor_microsatellites.list").getVal()
params.gps = Channel.fromPath("$params.refDir/mutect2_GetPileupSummaries.vcf.gz").getVal()
params.gpstbi = Channel.fromPath("$params.refDir/mutect2_GetPileupSummaries.vcf.gz.tbi").getVal()

Channel.fromPath("$params.refDir/pcgr/", type: 'dir').into { CPSR; PCGR }

/* -1: Install scripts required if not extant
*/
process scrpts {

  publishDir "$params.scriptDir", mode: "copy", pattern: "*"

  output:
  file('facets_cna.call.R') into facetscallscript
  set file('facets_cna_consensus.call.R'), file('facets_cna_consensus.func.R') into facetsconscript
  file('filterMuTect2Format.pl') into filtermutect2script
  file('filterPiscesFormat.pl') into filterpiscesscript
  file('MuTect2_contamination.call.R') into mutect2contamscript
  file('variants_GRanges*.R') into variantsGRangesscript
  file('vcf42.head.txt') into pcgrVcfHead

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

/* 0.000: Input using sample.csv
*/
Channel.fromPath("$params.sampleCsv", type: 'file')
       .splitCsv( header: true )
       .map { row -> [row.caseID, row.sampleID, row.meta, file(row.read1), file(row.read2)] }
       .set { bbduking }

/* 0.0: Input trimming
*/
process bbduke {

  label 'c10_30G_cpu_mem'

  publishDir "$params.outDir/$caseID/$sampleID/bbduk", mode: "copy", pattern: "*.txt"

  input:
  set val(caseID), val(sampleID), val(meta), file(read1), file(read2) from bbduking

  output:
  set val(caseID), val(sampleID), val(meta), file(read1), file(read2) into fastpingpre
  set val(caseID), val(sampleID), val(meta), file('*.bbduk.R1.fastq.gz'), file('*.bbduk.R2.fastq.gz') into (bwa_memming, fastpingpost)

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
  set val(caseID), val(sampleID), val(meta), file(read1), file(read2) from fastpingpre
  set val(caseID), val(sampleID), val(meta), file(read1), file(read2) from fastpingpost

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

  publishDir "$params.outDir/$caseID/$sampleID/bwa", mode: "copy", pattern: "*[cram, log.txt]"

  input:
  set val(caseID), val(sampleID), val(meta), file(read1), file(read2) from bwa_memming
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  set file(amb), file(ann), file(bwt), file(pac), file(sa) from Channel.value([params.amb, params.ann, params.bwt, params.pac, params.sa])

  output:
  set val(caseID), val(sampleID), val(meta), file('*.bam'), file('*.bai') into dup_marking

  script:
  """
  DATE=\$(date +"%Y-%m-%dT%T")
  RGLINE="@RG\\tID:$sampleID\\tPL:ILLUMINA\\tSM:$sampleID\\tDS:tumour\\tCN:UCD\\tLB:LANE_X\\tDT:\$DATE"

  {
    bwa mem \
    -t ${task.cpus} \
    -M \
    -R \$RGLINE \
    $fasta \
    $read1 $read2 | \
    samtools sort -T "tmp."$sampleID -o $sampleID".sort.bam"

  samtools index $sampleID".sort.bam"

  samtools view -hC -T $fasta $sampleID".sort.bam" > $sampleID".sort.cram"
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
  set val(caseID), val(sampleID), val(meta), file(bam), file(bai) from dup_marking

  output:
  file('*md.metrics.txt') into mrkdup_multiqc
  set val(caseID), val(sampleID), val(meta), file('*.md.bam'), file('*.md.bam.bai') into gatk4recaling

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
*/
process gtkrcl {

  label 'c10_30G_cpu_mem'

  publishDir "$params.outDir/$caseID/$sampleID/gatk4", mode: "copy", pattern: "*.txt"

  input:
  set val(caseID), val(sampleID), val(meta), file(bam), file(bai) from gatk4recaling
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  set file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])
  file(exomeintlist) from Channel.value(params.exomeintlist)

  output:
  file('*.table') into gtkrcl_multiqc
  set val(caseID), val(sampleID), val(meta), file('*.bqsr.bam'), file('*.bqsr.bam.bai') into (gatk_germ, multimetricing, mutect2ing, facetsing, msisensoring, mantastrelka2ing, piscesing)

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

  samtools index \$OUTBAM
  } 2>&1 | tee > $sampleID".GATK4_recal.log.txt"
  """
}

/* 1.21: GATK4 Germline
*/
process gatkgerm {

  label 'c40_120G_cpu_mem'

  publishDir "$params.outDir/$caseID/$sampleID/gatk4/HC_germline", mode: "copy", pattern: "*"

  input:
  set val(caseID), val(sampleID), val(meta), file(bam), file(bai) from gatk_germ
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  set file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])
  file(exomeintlist) from Channel.value(params.exomeintlist)
  set file(omni), file(otbi), file(kgp1), file(ktbi), file(hpmp), file(htbi), file(mlls), file(mtbi) from Channel.value([params.omni, params.otbi, params.kgp1, params.ktbi, params.hpmp, params.htbi, params.mlls, params.mtbi])

  output:
  set val(sampleID), val(meta), file('*recal_filt.vcf.gz'), file('*recal_filt.vcf.gz.tbi') into germ_vcf

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
      -O input.recal.vcf \
      -L $exomeintlist

    gatk --java-options ${params.full_javamem} VariantRecalibrator \
      -R $fasta \
      --variant input.recal.vcf \
      --resource:hapmap,known=false,training=true,truth=true,prior=15.0 $hpmp \
      --resource:omni,known=false,training=true,truth=false,prior=12.0 $omni \
      --resource:KG,known=false,training=true,truth=false,prior=10.0 $kgp1 \
      --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $dbsnp \
      -an QD -an DP -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
      --mode SNP \
      --output "SNP.recal.vcf" \
      --tranches-file "SNP.tranches" \
      --rscript-file "SNP.recal.plots.R"

    gatk --java-options ${params.full_javamem} ApplyVQSR \
      -R $fasta \
      --variant "input.recal.vcf" \
      --truth-sensitivity-filter-level 99.0 \
      --tranches-file "SNP.tranches" \
      --recal-file "SNP.recal.vcf" \
      --mode SNP \
      --output "SNP.recal_filt.vcf"

    ##select all indels
    gatk --java-options ${params.full_javamem} SelectVariants \
      -R $fasta \
      --variant input.recal.vcf \
      --select-type-to-include INDEL \
      --output "INDEL.recal.vcf"

    gatk --java-options ${params.full_javamem} MergeVcfs \
     -I "SNP.recal_filt.vcf" \
     -I "INDEL.recal.vcf" \
     -O $sampleID".recal_filt.vcf"

    uniq $sampleID".recal_filt.vcf" > 1; mv 1 $sampleID".recal_filt.vcf"
    bgzip $sampleID".recal_filt.vcf"; tabix $sampleID".recal_filt.vcf.gz"

  } 2>&1 | tee $sampleID".GATK4_HaplotypeCaller-germline.log.txt"

  """
}

/* 1.22: CPSR annotation of GATK4 Germline
*/
process cpsrreport {

  label 'c40_120G_cpu_mem'

  publishDir "$params.outDir/calls/reports", mode: "copy", pattern: "*.html"
  publishDir "$params.outDir/calls/variants/cpsr", mode: "copy", pattern: "*[!.html]"

  input:
  set val(sampleID), val(meta), file(vcf), file(tbi) from germ_vcf
  each file(cpsr_grch37) from CPSR

  output:
  file('*') into cpsr_vcfs

  script:
  """
  {
    cpsr.py \
      --no-docker \
      --no_vcf_validate \
      $vcf \
      $cpsr_grch37 \
      ./ \
      grch37 \
      0 \
      $cpsr_grch37"/data/grch37/cpsr_configuration_default.toml" \
      $meta

  } 2>&1 | tee $sampleID".cpsr.log.txt"
  """
}

/* 2.0: Metrics suite, this will produce an HTML report
*/
process mltmet {

  label 'c8_24G_cpu_mem'

  publishDir "$params.outDir/$caseID/$sampleID/metrics"

  input:
  set val(caseID), val(sampleID), val(meta), file(bam), file(bai) from multimetricing
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
      BAIT_INTERVALS=$exomeintlist \
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

  } 2>&1 | tee > $sampleID".picard-metrics.log.txt"
  """
}

/* 2.2: MSIsensor
*/
process msisen {

  label 'c10_30G_cpu_mem'

  publishDir "$params.outDir/$sampleID/msisensor"
  publishDir "$params.outDir/calls/msisensor", pattern: '*.txt'

  input:
  set val(caseID), val(sampleID), val(meta), file(tumourbam), file(tumourbai) from msisensoring
  file(exomebed) from Channel.value(params.exomebed)
  file(ssrs) from Channel.value(params.ssrs)

  output:
  val(sampleID) into completed2_2
  file('*') into msisensoroutput

  script:
  """
  msisensor msi \
    -d $ssrs \
    -t $tumourbam \
    -e $exomebed \
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
  set val(caseID), val(sampleID), val(meta), file(tumourbam), file(tumourbai) from mutect2ing
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  file(exomeintlist) from Channel.value(params.exomeintlist)
  set file(gps), file(gpstbi) from Channel.value([params.gps, params.gpstbi])
  file(filterpl) from filtermutect2script

  output:
  file('*.pass.vcf') into mutect2_veping mode flatten
  file('*.raw.vcf') into mutect2_rawVcf
  file('*') into completedmutect2call
  set val(sampleID), val(meta), file('*calculatecontamination.table') into contamination

  script:
  """
  {
    gatk --java-options ${params.full_javamem} \
      Mutect2 \
      --native-pair-hmm-threads ${task.cpus} \
      --reference $fasta \
      --input $tumourbam \
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

   label 'c10_30G_cpu_mem'

   publishDir "$params.outDir/", pattern: '*issue.table'

   input:
   set val(sampleID), val(meta), file(contable) from contamination
   file script from mutect2contamscript

   output:
   file('*.table') into completedcontam

   """
   Rscript --vanilla $script $contable $sampleID
   """
}


/* 2.5: Pisces
*/
process pisces {

  label 'c40_120G_cpu_mem'

  publishDir "$params.outDir/$sampleID/pisces"
  publishDir "$params.outDir/calls/variants/vcf", pattern: '*raw.vcf'

  input:
  set val(caseID), val(sampleID), val(meta), file(tumourbam), file(tumourbai) from piscesing
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  file(exomeintlist) from Channel.value(params.exomeintlist)
  file(filterpl) from filterpiscesscript

  output:
  file('*.pass.vcf') into pisces_veping mode flatten
  file('*.raw.vcf') into pisces_rawVcf
  file('*') into completedpiscescall
  file('*.sampleIDmeta.csv') into pcgrsamples

  script:
  """
  TUMOURVCF=\$(echo $tumourbam | sed 's/bam/vcf/')
  {
    dotnet /app/CreateGenomeSizeFile_5.2.9.122/CreateGenomeSizeFile.dll -g ./ -s "Homo sapiens (hg19)"

    dotnet /app/Pisces_5.2.9.122/Pisces.dll \
      --bam $tumourbam \
      --genomepaths ./ \
      --intervalpaths $exomeintlist \
      --filterduplicates true \
      --onlyuseproperpairs true \
      --mindepth 3 \
      --gvcf false \
      --outfolder ./ \
      --maxthreads ${task.cpus}

    perl $filterpl \
      ID=$sampleID \
      DP=14 \
      MD=2 \
      VCF=\$TUMOURVCF

    echo "$sampleID,$meta" > $sampleID".sampleIDmeta.csv"
  } 2>&1 | tee > $sampleID".pisces.log.txt"

  """
}

/* 3.0: Annotate Vcfs
*/
ALLVCFS = pisces_veping
          .mix( mutect2_veping )

process vepann {

  label 'c10_30G_cpu_mem'

  publishDir "$params.outDir/calls/variants/vcf", pattern: '*.[vcf,tsv]'

  input:
  each file(vcf) from ALLVCFS
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])

  output:
  file('*.vcf') into annoVcfs
  file('*.tsv') into annoTsvs
  file('*.vcf') into runGRanges

  script:
  """
  VCFANNO=\$(echo $vcf | sed "s/.vcf/.vep.vcf/")
  TSVANNO=\$(echo $vcf | sed "s/.vcf/.vep.tsv/")
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

  cat \$VCFANNO | perl -ane '@s=split(/\\|/,\$F[7]);
    print "\$F[0]:\$F[1]_\$F[3]>\$F[4]\\t\$s[1]\\t\$s[2]\\t\$s[3]\\t\$s[4]\\t\$s[45]\\t\$s[64]\\t\$F[9]\\n";' > \$TSVANNO
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

vartypes = Channel.from( "snv" )

process vcfGRa {

  label 'c20_60G_cpu_mem'

  publishDir "$params.outDir/calls/variants/pdf", pattern: '*.pdf'
  publishDir "$params.outDir/calls/variants/vcf", pattern: '*.vcf'
  publishDir "$params.outDir/calls/variants/data", pattern: '*.[*RData,*tab]'

  input:
  file(rawGRangesvcff) from ALLRAWVEPVCFS.collect()
  each vartype from vartypes
  set file(callR), file(funcR) from variantsGRangesscript
  file(vcfHead) from pcgrVcfHead

  output:
  val(vartype) into completed3_1
  file('*.ALL.pcgr.all.tab.vcf') into pcgrvcfs
  file('*') into completedvcfGRangesConsensus

  script:
  """
  OUTID=\$(basename ${params.runDir})
  Rscript --vanilla $callR \
    $funcR \
    NULL \
    $vartype".pass.vep.vcf" \
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

  label 'c20_60G_cpu_mem'

  publishDir "$params.outDir/calls/reports", mode: "copy", pattern: "*html"
  publishDir "$params.outDir/calls/variants/pcgr", mode: "copy", pattern: "*[!.html]"

  input:
  file(vcf) from pcgrvcfs
  each file(pcgr_grch37) from PCGR
  each file(pcgrmeta) from pcgrsamples

  output:
  file('*') into completedPCGR

  script:
  """
  {
    sed 's/vcf_tumor_only = false/vcf_tumor_only = true/' $pcgr_grch37/data/grch37/pcgr_configuration_default.toml > pcgr_TO.toml

    SAMPLEID=\$(cut -d "," -f 1 $pcgrmeta | cut -d "." -f 1)
    VCF=\$(ls | grep \$SAMPLEID | grep pcgr.all.tab.vcf)
    METAID=\$(cut -d "," -f 2 $pcgrmeta)
      pcgr.py $pcgr_grch37 \
      ./ \
      grch37 \
      pcgr_TO.toml \
      \$METAID \
      --input_vcf \$VCF \
      --no-docker \
      --force_overwrite \
      --no_vcf_validate
  } 2>&1 | tee pcgr.log.txt
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
