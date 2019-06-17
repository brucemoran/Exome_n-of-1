#!/usr/bin/env nextflow

params.help = ""

if (params.help) {
  log.info ''
  log.info ' -------------------------------------------------------------'
  log.info '| BATCH TUMOUR-NORMAL FASTQ QC, TRIM, ALIGN, SOMATIC SNV, CNA |'
  log.info ' -------------------------------------------------------------'
  log.info ''
  log.info 'Usage: '
  log.info 'nextflow run batch_exome_tumour-normal.simg.nf \
              --sampleCsv "sample.csv" \
              --refDir "refs"'
  log.info ''
  log.info 'This batch processing script only allows one tumour/matched-normal pair; these are processed as per standard, but NB different input format'
  log.info ''
  log.info 'Mandatory arguments:'
  log.info '    --sampleCsv      STRING      CSV format, headers: caseID,soma_sampleID,soma_read1,soma_read2,germ_sampleID,germ_read1,germ_read2'
  log.info '    --refDir      STRING      dir in which reference data and required indices held; recommended to run associated reference creation NextFlow, DNAseq_references; still stuck on hg19 for several reasons;()'
  log.info ''
  exit 1
}

/* -2: Global Variables
*/
params.runDir = "$workflow.launchDir"
params.outDir = "$params.runDir/analysis"
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
  file('*') into completedmin1
  file('facets_cna.call.R') into facetscallscript
  file('facets_cna_consensus.call.R') into facetsconcscript
  file('facets_cna_consensus.func.R') into facetsconfscript
  file('filterLancetSomaticFormat.pl') into filterlancetscript
  file('filterMuTect2Format.pl') into filtermutect2script
  file('filterStrelka2IndelFormat.pl') into filterstrelka2iscript
  file('filterStrelka2SNVFormat.pl') into filterstrelka2sscript
  file('MuTect2_contamination.call.R') into mutect2contamscript
  file('variants_GRanges_consensus_plot_batch*.R') into variantsGRangesscript
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
completedmin1.subscribe { println "Scripts, images output to: $params.scriptDir" }

/* 0.000: Input using sample.csv
*/
Channel.fromPath("$params.sampleCsv", type: 'file')
       .splitCsv( header: true )
       .map { row -> [row.caseID, row.soma_sampleID, file(row.soma_read1), file(row.soma_read2), row.germ_sampleID, file(row.germ_read1), file(row.germ_read2)] }
       .set { split_soma_germ }

/* 0.00: Input using sample.csv
*/
process splt_sg {

  publishDir "$params.outDir/$caseID", mode: "copy", pattern: "*.csv"

  input:
  set val(caseID), val(soma_sampleID), file(soma_read1), file(soma_read2), val(germ_sampleID), file(germ_read1), file(germ_read2) from split_soma_germ

  output:
  file('*.csv') into splitcsv2

  """
  SR1=\$(readlink -e $soma_read1)
  SR2=\$(readlink -e $soma_read2)
  GR1=\$(readlink -e $germ_read1)
  GR2=\$(readlink -e $germ_read2)
  echo "caseID,type,sampleID,read1,read2" > $caseID".csv"
  echo "$caseID,somatic,$soma_sampleID,\$SR1,\$SR2" >> $caseID".csv"
  echo "$caseID,germline,$germ_sampleID,\$GR1,\$GR2" >> $caseID".csv"
  """
}

/*two set channels need to be processed the same; use
*/
splitcsv2.splitCsv( header: true )
         .map { row -> [row.caseID, row.type, row.sampleID, file(row.read1), file(row.read2)] }
         .set { bbduking }

/* 0.0: Input trimming
*/
process bbduke {

  label 'c10_30G_cpu_mem'

  publishDir "$params.outDir/$caseID/$sampleID/bbduk", mode: "copy", pattern: "*.txt"

  input:
  set val(caseID), val(type), val(sampleID), file(read1), file(read2) from bbduking

  output:
  val(sampleID) into completed0_0
  set val(caseID), val(type), val(sampleID), file(read1), file(read2) into fastping
  set val(caseID), val(type), val(sampleID), file('*.bbduk.R1.fastq.gz'), file('*.bbduk.R2.fastq.gz') into bwa_memming

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
completed0_0.subscribe { println "Completed BBDuk: " + it }

process fastp {

  label 'c8_24G_cpu_mem'

  publishDir "$params.outDir/$caseID/$sampleID/fastp", mode: "copy", pattern: "*.html"

  input:
  set val(caseID), val(type), val(sampleID), file(read1), file(read2) from fastping

  output:
  file('*.html') into completed0_1
  file('*.json') into fastp_multiqc

  script:
  """
  fastp -w ${task.cpus} -h $sampleID".fastp.html" -j $sampleID".fastp.json" --in1 $read1 --in2 $read2
  """
}
completed0_1.subscribe { println "Completed Fastp: " + it }

/* 1.0: Input alignment
*/
process bwamem {

  label 'c20_60G_cpu_mem'

  publishDir "$params.outDir/$caseID/$sampleID/bwa", mode: "copy", pattern: "*[!bam]"

  input:
  set val(caseID), val(type), val(sampleID), file(read1), file(read2) from bwa_memming
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  set file(amb), file(ann), file(bwt), file(pac), file(sa) from Channel.value([params.amb, params.ann, params.bwt, params.pac, params.sa])

  output:
  val(sampleID) into completed1_0
  set val(caseID), val(type), val(sampleID), file('*.bam'), file('*.bai') into dup_marking

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

  samtools view -hC -T $fasta $sampleID".sort.bam" > $sampleID".sort.cram"
  } 2>&1 | tee > $sampleID".bwa-mem.log.txt"
  """
}
completed1_0.subscribe { println "Completed bwa-mem: " + it }

/* 1.1: MarkDuplicates
*/
process mrkdup {

  label 'c40_120G_cpu_mem'

  publishDir "$params.outDir/$caseID/$sampleID/picard/markdup", mode: "copy", pattern: "*.log.txt"
  publishDir "$params.outDir/$caseID/$sampleID/picard/metrics", mode: "copy", pattern: "*.metrics.txt"

  input:
  set val(caseID), val(type), val(sampleID), file(bam), file(bai) from dup_marking

  output:
  val(sampleID) into completed1_1
  file('*md.metrics.txt') into mrkdup_multiqc
  set val(caseID), val(type), val(sampleID), file('*.md.bam'), file('*.md.bam.bai') into gatk4recaling

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
completed1_1.subscribe { println "Completed MarkDuplicates: " + it }

/* 1.2: GATK4 BestPractices
* as per best practices of GATK4
*/
process gtkrcl {

  label 'c10_30G_cpu_mem'

  publishDir "$params.outDir/$caseID/$sampleID/gatk4", mode: "copy", pattern: "*.txt"

  input:
  set val(caseID), val(type), val(sampleID), file(bam), file(bai) from gatk4recaling
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  set file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])
  file(exomeintlist) from Channel.value(params.exomeintlist)

  output:
  val(sampleID) into completed1_2
  file('*.table') into gtkrcl_multiqc
  set val(caseID), val(type), val(sampleID), file('*.bqsr.bam') into (germfiltering, somafiltering)
  set val(caseID), val(type), val(sampleID), file(bam) into gatk_germ

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
completed1_2.subscribe { println "Completed GATK4 BaseRecalibration: " + it }

/* 1.21: GATK4 Germline
*/
// process gatkgerm {
//
//   label 'c40_120G_cpu_mem'
//
//   publishDir "$params.outDir/$caseID/$sampleID/gatk4/HC_germline", mode: "copy", pattern: "*"
//
//   input:
//   set val(type), val(sampleID), file(bam) from gatk_germ
//   set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
//   set file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])
//   file(exomeintlist) from Channel.value(params.exomeintlist)
//   set file(omni), file(otbi), file(kgp1), file(ktbi), file(hpmp), file(htbi) from Channel.value([params.omni, params.otbi, params.kgp1, params.ktbi, params.hpmp, params.htbi])
//
//   output:
//   set val(sampleID), file('*.HC.vcf.gz'), file('*.HC.vcf.gz.tbi') into germ_vcf
//
//   script:
//   """
//   {
//     #HaplotypeCaller
//     INPUTBAM=$bam
//     OUTVCF=\$(echo \$INPUTBAM | sed 's/bam/hc.vcf/')
//     gatk --java-options ${params.full_javamem} HaplotypeCaller \
//       -R $fasta \
//       -I \$INPUTBAM \
//       --dont-use-soft-clipped-bases \
//       --standard-min-confidence-threshold-for-calling 20 \
//       --dbsnp $dbsnp \
//       --native-pair-hmm-threads ${task.cpus} \
//       -O $sampleID".HC.vcf" \
//       -L $exomeintlist
//
//     bgzip $sampleID".HC.vcf"
//     tabix $sampleID".HC.vcf.gz"
//
//   } 2>&1 | tee $sampleID".GATK4_HaplotypeCaller-germline.log.txt"
//
//   """
// }
//
// /* 1.22: CPSR annotation of GATK4 Germline
// */
// process cpsrreport {
//
//   label 'c40_120G_cpu_mem'
//
//   publishDir "$params.outDir/calls/reports", mode: "copy", pattern: "*html"
//   publishDir "$params.outDir/calls/variants/pcgr", mode: "copy", pattern: "*[!.html]"
//
//   input:
//   set val(sampleID), file(vcf), file(tbi) from germ_vcf
//   file(cpsr_grch37) from CPSR
//
//   output:
//   file('*') into cpsr_vcfs
//
//   script:
//   """
//   {
//     ##activate the conda env
//     . activate pcgr
//
//     cpsr.py $cpsr_grch37 \
//       ./ \
//       grch37 \
//       $cpsr_grch37"/data/grch37/cpsr_configuration_default.toml" \
//       $sampleID \
//       --input_vcf $vcf \
//       --no-docker
//
//   } 2>&1 | tee $sampleID".cpsr.log.txt"
//   """
// }

/* 1.3: filter germ into a channel, index bam
*/
process grmflt {

  input:
  set val(caseID), val(type), val(sampleID), file(bam) from germfiltering

  output:
  set val(caseID), val(sampleID), file(bam), file ('*.bam.bai') into (gmultimetricing, gmutect2somaticing, gfacetsomaing, gmsisensoring, gmantastrelka2ing, glanceting)
  val(sampleID) into vcfGRaID

  when:
  type == "germline"

  script:
  """
  samtools index $bam > $bam".bai"
  """
}

/* 1.4: filter somatic into a channel, index bam
*/
process smaflt {

  input:
  set val(caseID), val(type), val(sampleID), file(bam) from somafiltering

  output:
  set val(caseID), val(sampleID), file(bam), file ('*.bam.bai') into (multimetricing, mutect2somaticing, facetsomaing, msisensoring, mantastrelka2ing, lanceting)

  when:
  type != "germline"

  script:
  """
  samtools index $bam > $bam".bai"
  """
}

/* join all s, g by caseID
*/
facetsomaing.join(gfacetsomaing).set { facetsing }
msisensoring.join(gmsisensoring).set { msisensor }
mutect2somaticing.join(gmutect2somaticing).set { mutect2somatic }
mantastrelka2ing.join(gmantastrelka2ing).set { mantastrelka2 }
lanceting.join(glanceting).set { lancet }

/* 2.0: Metrics suite, this will produce an HTML report
*/
MULTIALL = gmultimetricing.mix(multimetricing)

process mltmet {

  label 'c8_24G_cpu_mem'

  publishDir "$params.outDir/$caseID/$sampleID/metrics"

  input:
  set val(caseID), val(sampleID), file(bam), file(bai) from MULTIALL
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  set file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])
  file(exomeintlist) from Channel.value(params.exomeintlist)

  output:
  val(sampleID) into completed2_0
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

  } 2>&1 | tee > $sampleID".picard-metrics.log.txt"
  """
}
completed2_0.subscribe { println "Completed running metrics" }

/*2.1: SCNA with facets CSV snp-pileup
*/
process fctcsv {

  label 'c8_24G_cpu_mem'

  publishDir "$params.outDir/$caseID/$sampleID/facets"

  input:
  set val(caseID), val(sampleID), file(bam), file(bai), val(germlineID), file(germlinebam), file(germlinebai) from facetsing
  set file(dbsnp), file(dbsnptbi) from Channel.value([params.dbsnp, params.dbsnptbi])
  file(facetsR) from facetscallscript

  output:
  val(sampleID) into completed2_1
  file('*') into facetsoutputR

  script:
  """
  CSVFILE=\$(echo $bam | sed 's/bam/facets.r10.csv/')

  {
    snp-pileup \
      $dbsnp \
      -r 10 \
      -p \
      \$CSVFILE \
      $germlinebam \
      $bam

    Rscript --vanilla $facetsR \$CSVFILE
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
  set val(caseID), val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from msisensor
  file(exomebed) from Channel.value(params.exomebed)
  file(ssrs) from Channel.value(params.ssrs)

  output:
  file('*') into msisensoroutput

  script:
  """
  msisensor msi \
    -d $ssrs \
    -n $germlinebam \
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
  set val(caseID), val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from mutect2somatic
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  file(exomeintlist) from Channel.value(params.exomeintlist)
  set file(gps), file(gpstbi) from Channel.value([params.gps, params.gpstbi])
  file(filterpl) from filtermutect2script

  output:
  val(sampleID) into completed2_3
  file('*.pass.vcf') into mutect2_veping mode flatten
  file('*.raw.vcf') into mutect2_rawVcf
  file('*') into completedmutect2call
  set val(sampleID), file('*calculatecontamination.table') into contamination

  script:
  """
  GPS=$gps
  if [[ \$GPS =~ "tbi" ]];then
    GPS=\$(echo $gps | sed 's/.tbi//')
  fi
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
      -V \$GPS \
      -O $sampleID".getpileupsummaries.table" \
      -L $exomeintlist

    gatk CalculateContamination \
      -I $sampleID".getpileupsummaries.table" \
      -O $sampleID".calculatecontamination.table"

    if [[ \$(grep NaN $sampleID".calculatecontamination.table" | wc -l) == 1 ]]; then
      USECONT=""
    else
      USECONT="--contamination-table ${sampleID}.calculatecontamination.table"
    fi

    gatk --java-options ${params.full_javamem} \
      FilterMutectCalls \
      --reference $fasta \
      --interval-padding 5 \
      --output $sampleID".md.recal.mutect2.FilterMutectCalls.vcf" \
      --unique-alt-read-count 3 \
      --variant $sampleID".md.recal.mutect2.vcf" \
      -L $exomeintlist \$USECONT

    perl $filterpl \
      ID=$sampleID \
      DP=14 \
      MD=2 \
      VCF=$sampleID".md.recal.mutect2.FilterMutectCalls.vcf"

  } 2>&1 | tee > $sampleID".GATK4_mutect2.log.txt"
  """

}
completed2_3.subscribe { println "Completed Mutect2: " + it }

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
  set val(caseID), val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from mantastrelka2
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  set file(exomebedgz), file(exomebedgztbi) from Channel.value([params.exomebedgz, params.exomebedgztbi])
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
completed2_4.subscribe { println "Completed Manta-Strelka2: " + it }

/* 2.5: Lancet
*/
process lancet {

  label 'c40_120G_cpu_mem'

  publishDir "$params.outDir/$sampleID/lancet"
  publishDir "$params.outDir/calls/variants/vcf", pattern: '*raw.vcf'

  input:
  set val(caseID), val(sampleID), file(tumourbam), file(tumourbai), val(germlineID), file(germlinebam), file(germlinebai) from lancet
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])
  file(exomebed) from Channel.value(params.exomebed)
  file(filterLancet) from filterlancetscript

  output:
  val(sampleID) into completed2_5
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
completed2_5.subscribe { println "Completed lancet: " + it }

/* 3.0: Annotate Vcfs
*/
ALLVCFS = lancet_veping
          .mix( mutect2_veping )
          .mix( strelka2_veping )

process vepann {

  label 'c10_30G_cpu_mem'

  publishDir "$params.outDir/calls/variants/vcf", pattern: '*.vcf'

  input:
  each file(vcf) from ALLVCFS
  set file(fasta), file(fai), file(dict) from Channel.value([params.fasta, params.fai, params.dict])

  output:
  val("VEP") into completed3_0
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
completed3_0.subscribe { println "Completed " + it + " annotation" }

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
  publishDir "$params.outDir/calls/variants/data", pattern: '*.{RData,tab}'

  input:
  file(rawGRangesvcff) from ALLRAWVEPVCFS.collect()
  set file(callR), file(funcR) from variantsGRangesscript
  file(vcfHead) from pcgrVcfHead

  output:
  val(vartype) into completed3_1
  file('*.pcgr.all.tab.vcf') into pcgrvcfs
  file('*') into completedvcfGRangesConsensus

  script:
  """
  OUTID=\$(basename ${params.runDir})
  Rscript --vanilla $callR \
    $funcR \
    ${params.tumourPattern} \
    "snv.pass.vep.vcf" \
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
  file(pcgr_grch37) from PCGR

  output:
  file('*') into completedPCGR

  script:
  """
  ##activate the conda env
  . activate pcgr
  SAMPLEID=\$(echo $vcf | cut -d "." -f 1)
  pcgr.py $pcgr_grch37 \
    ./ \
    grch37 \
    $pcgr_grch37/data/grch37/grch37/pcgr_configuration_default.toml \
    \$SAMPLEID \
    --input_vcf $vcf \
    --no-docker \
    --force_overwrite
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
completedmultiqc.subscribe { println "Completed MultiQC" }
