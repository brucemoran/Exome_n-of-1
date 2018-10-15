#! /bin/bash

##consensus facets CNA on nextflow
DICTFILE=\$(echo ${params.genomefa} | cut -d "." -f 1)".dict"
Rscript --vanilla ${workflow.projectDir}/templates/facets_cna_consensus.caller.R \
  ${params.rLibPath} \
  ${params.refDir}/\$DICTFILE \
  ${params.refDir}/COSMIC_CGC.bed \
  ${params.runID} \
  ${workflow.projectDir}/templates/facets_cna_consensus.func.R
