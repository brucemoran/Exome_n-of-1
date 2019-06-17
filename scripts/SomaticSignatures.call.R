#! R

##inputs
libs <- c("SomaticSignatures", "MutationalPatterns", "GenomicRanges", "bedr", "ggplot2", "rtracklayer", "tidyverse")
libsLoaded <- lapply(libs,function(lib){suppressMessages(library(lib, character.only = TRUE))})
options(scipen=999)
argsIn <- commandArgs(trailingOnly = TRUE)

##input is single commma-sep argument of all VCFs to import
##NB should be cat'd SNV, indel
vcfList <- as.list(dir(pattern=argsIn[1]))
twobitref <- TwoBitFile(argsIn[2])

##make ref, alt SNV
refAlleles <- altAlleles <- c()
sampleNames <- unlist(lapply(vcfList,function(x){
  strsplit(x, "\\.")[[1]][1]
}))

for(x in 1:length(vcfList)){
  vcfList[[x]] <- read.table(vcfList[[x]])
  vcfList[[x]][,3] <- vcfList[[x]][,2]
  refAlleles <- as.vector(vcfList[[x]][,4])
  altAlleles <- as.vector(vcfList[[x]][,5])
  colnames(vcfList[[x]])<-c("chr", "start", "end", colnames(vcfList[[x]][c(4,5,6,7,8,9,10)]))
  dfin <- data.frame(sampleNames=rep(sampleNames[x],
                     times=length(refAlleles)),
                     refAlleles=refAlleles,
                     altAlleles=altAlleles)
  VRangeIN <- VRanges(seqnames=Rle(vcfList[[x]]$chr),
                      ranges=IRanges(vcfList[[x]]$start,end=vcfList[[x]]$end),
                      ref=factor(refAlleles),
                      alt=factor(altAlleles),
                      sampleNames=factor(rep(sampleNames[x])))
  VRangeIN <- unname(VRangeIN)
  if(x==1){
    VRangeUse <- VRangeIN
  }
  if(x>1){
    VRangeUse <- c(VRangeIN, VRangeUse)
  }
}

VRangeUse <- sort(VRangeUse)

##motifs
use_motifs <- mutationContext(VRangeUse, ref=twobitref)
pdf("mutation_spectrum_somaticSignatures.pdf")
  plotMutationSpectrum(use_motifs, "sampleNames")
dev.off()

##decompose
use_mmn <- motifMatrix(use_motifs, group = "sampleNames", normalize = TRUE)
##counts
use_mmf <- motifMatrix(use_motifs, group = "sampleNames", normalize = FALSE)

##no zeros
nzero <- function(x){x[apply(x,1,sum)>0,]}
use_mmnzr <- nzero(use_mmn)

##NMF and decompose to find best number of signatures, plot out
n_sigs = 2:dim(use_mmnzr)[2]
gof_nmf = assessNumberSignatures(use_mmnzr, n_sigs, nReplicates = 5)

##NMF signature def based on number of samples
sigs_nmf = identifySignatures(use_mmnzr, max(n_sigs), nmfDecomposition)

##plots
pdf(paste0("identifySignatures.nSigs_", paste0(max(n_sigs)), "-nReps_5.somaticSignatures.pdf"))
  plotObservedSpectrum(sigs_nmf,colorby=c("alteration"))
dev.off()

pdf(paste0("sampleMap.nSigs_", paste0(max(n_sigs)), "-nReps_5.somaticSignatures.pdf"))
  plotSampleMap(sigs_nmf)
dev.off()

p <- plotSamples(sigs_nmf) +
  theme(legend.position = "right") +
  xlab("Samples") +
  ggtitle("Somatic Signatures") +
  scale_fill_manual(values=rainbow(19)) +
  theme(axis.text.x = element_text(size = 9))
ggsave(filename=paste0("samples_propbar.nSigs-", max(n_sigs), ".SomaticSignatures.pdf"))

##MutationalPatterns

##function to convert context, mutation from use_mmnzr
contextInter <- function(f){
  ul <- unlist(strsplit(strsplit(f," ")[[1]],""))
  paste0(ul[3],"[",ul[1],">",ul[2],"]",ul[5])
}

##create context to match below input from MutationlPatterns COSMIC data
use_mmf_ci <- as_tibble(use_mmf, rownames="tri_con") %>%
              dplyr::mutate(tri_conmu = unlist(lapply(tri_con, contextInter)))
use_mmf_ci_df <- use_mmf_ci %>%
                 dplyr::select(tri_conmu, 2:c(dim(use_mmf_ci)[2]-1)) %>%
                 dplyr::arrange(tri_conmu) %>%
                 column_to_rownames("tri_conmu")

##from https://github.com/mskcc/mutation-signatures/blob/master/signature_analysis.R
##names for signatures
sign_map <- c('Age','APOBEC.1','BRCA1/2','Smoking','Sig_5', 'MMR.1','UV','Sig_8','IGHV_hypermut','POLE','TMZ','Sig_12','APOBEC.2','Sig_14','MMR.2','Sig_16','Sig_17','Sig_18','Sig_19','MMR.3','Sig_21','Aristolochic_acid','Sig_23','Sig_24','Sig_25','MMR.4','Sig_27','Sig_28','Tobacco', 'Sig30')

##read input from cosmic
sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/",
                "signatures_probabilities.txt", sep = "")
cancer_signatures <- read.table(sp_url, sep = "\t", header = TRUE)
new_order <- match(row.names(use_mmf_ci_df), cancer_signatures$Somatic.Mutation.Type)
cancer_signatures <- cancer_signatures[as.vector(new_order),]
row.names(cancer_signatures) <- cancer_signatures$Somatic.Mutation.Type
cosmic_signatures <- as.matrix(cancer_signatures[,4:33])
colnames(cosmic_signatures) <- sign_map
cos_sim_samples_signatures <- cos_sim_matrix(use_mmf_ci_df, cosmic_signatures)

##order clusters
hclust_cosmic <- cluster_signatures(cosmic_signatures, method = "average")
# store signatures in new order
cosmic_order <- colnames(cosmic_signatures)[hclust_cosmic$order]
# fit
fit_res <- fit_to_signatures(use_mmf_ci_df, cosmic_signatures)

##plots
pdf("cosine_sim.heatmap.COSMIC.MutationalPatterns.pdf")
  plot_cosine_heatmap(cos_sim_samples_signatures,
                      col_order = cosmic_order,
                      cluster_rows = TRUE)
dev.off()

pdf("contribution_heatmap.COSMIC.MutationalPatterns.pdf")
  plot_contribution_heatmap(fit_res$contribution,
                          cluster_samples = TRUE,
                          method = "complete")
dev.off()

##save + exit
save(VRangeUse, use_mmn, use_mmf, sigs_nmf, use_mmf_ci_df, cosmic_signatures, cosmic_order, fit_res, file="signatures.RData")

quit(save = "no", status = 0)
