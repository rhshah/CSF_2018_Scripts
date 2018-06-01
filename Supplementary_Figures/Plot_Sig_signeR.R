library(signeR)
library(VariantAnnotation)
# BSgenome, equivalent to the one used on the variant call
library(BSgenome.Hsapiens.UCSC.hg19)
library(rtracklayer)

vcfobj <- readVcf("../NewRuns/signatures/s_IM_GBM_75_RZ_HC.vcf", "hg19")
mut <- genCountMatrixFromVcf(BSgenome.Hsapiens.UCSC.hg19, vcfobj)

target_regions <- import(con="/Users/shahr2/git/BIC-variants_pipeline/targets/IMPACT410_hg19/IMPACT410_hg19_baits.bed", format="bed")
opp <- genOpportunityFromGenome(BSgenome.Hsapiens.UCSC.hg19,
                                target_regions, nsamples=nrow(mut))
signatures <- signeR(M=mut, Opport=opp)
pdf(file = "../NewRuns/signatures/s_IM_GBM_75_RZ_HC_sig.pdf",
    width = 10,
    height = 2)
SignPlot(signatures$SignExposures)
dev.off()
vcfobj <- readVcf("../NewRuns/signatures/s_IM_GBM_75_RZ_HC.vcf", "hg19")
mut <- genCountMatrixFromVcf(BSgenome.Hsapiens.UCSC.hg19, vcfobj)

target_regions <- import(con="/Users/shahr2/git/BIC-variants_pipeline/targets/IMPACT410_hg19/IMPACT410_hg19_baits.bed", format="bed")
opp <- genOpportunityFromGenome(BSgenome.Hsapiens.UCSC.hg19,
                                target_regions, nsamples=nrow(mut))
signatures <- signeR(M=mut, Opport=opp)
pdf(file = "../NewRuns/signatures/s_IM_GBM_75_RZ_HC_sig.pdf",
    width = 10,
    height = 2)
SignPlot(signatures$SignExposures)
dev.off()
#
vcfobj <- readVcf("../NewRuns/signatures/s_IM_GBM_75_CSFa_HC.vcf", "hg19")
mut <- genCountMatrixFromVcf(BSgenome.Hsapiens.UCSC.hg19, vcfobj)

target_regions <- import(con="/Users/shahr2/git/BIC-variants_pipeline/targets/IMPACT410_hg19/IMPACT410_hg19_baits.bed", format="bed")
opp <- genOpportunityFromGenome(BSgenome.Hsapiens.UCSC.hg19,
                                target_regions, nsamples=nrow(mut))
signatures <- signeR(M=mut, Opport=opp)
pdf(file = "../NewRuns/signatures/s_IM_GBM_75_CSFa_HC_sig.pdf",
    width = 10,
    height = 2)
SignPlot(signatures$SignExposures)
dev.off()
#
vcfobj <- readVcf("../NewRuns/signatures/s_IM_GBM_75_CSFc_HC.vcf", "hg19")
mut <- genCountMatrixFromVcf(BSgenome.Hsapiens.UCSC.hg19, vcfobj)

target_regions <- import(con="/Users/shahr2/git/BIC-variants_pipeline/targets/IMPACT410_hg19/IMPACT410_hg19_baits.bed", format="bed")
opp <- genOpportunityFromGenome(BSgenome.Hsapiens.UCSC.hg19,
                                target_regions, nsamples=nrow(mut))
signatures <- signeR(M=mut, Opport=opp)
pdf(file = "../NewRuns/signatures/s_IM_GBM_75_CSFc_HC_sig.pdf",
    width = 10,
    height = 2)
SignPlot(signatures$SignExposures)
dev.off()
#
vcfobj <- readVcf("../NewRuns/signatures/s_IM_GBM_97_CSFa_HC.vcf", "hg19")
mut <- genCountMatrixFromVcf(BSgenome.Hsapiens.UCSC.hg19, vcfobj)

target_regions <- import(con="/Users/shahr2/git/BIC-variants_pipeline/targets/IMPACT410_hg19/IMPACT410_hg19_baits.bed", format="bed")
opp <- genOpportunityFromGenome(BSgenome.Hsapiens.UCSC.hg19,
                                target_regions, nsamples=nrow(mut))
signatures <- signeR(M=mut, Opport=opp)
pdf(file = "../NewRuns/signatures/s_IM_GBM_97_CSFa_HC_sig.pdf",
    width = 10,
    height = 2)
SignPlot(signatures$SignExposures)
dev.off()

#
vcfobj <- readVcf("../NewRuns/signatures/s_IM_GBM_111_RZ_HC.vcf", "hg19")
mut <- genCountMatrixFromVcf(BSgenome.Hsapiens.UCSC.hg19, vcfobj)

target_regions <- import(con="/Users/shahr2/git/BIC-variants_pipeline/targets/IMPACT410_hg19/IMPACT410_hg19_baits.bed", format="bed")
opp <- genOpportunityFromGenome(BSgenome.Hsapiens.UCSC.hg19,
                                target_regions, nsamples=nrow(mut))
signatures <- signeR(M=mut, Opport=opp)
pdf(file = "../NewRuns/signatures/s_IM_GBM_111_RZ_HC_sig.pdf",
    width = 10,
    height = 2)
SignPlot(signatures$SignExposures)
dev.off()
#
vcfobj <- readVcf("../NewRuns/signatures/s_IM_GBM_111_CSF_HC.vcf", "hg19")
mut <- genCountMatrixFromVcf(BSgenome.Hsapiens.UCSC.hg19, vcfobj)

target_regions <- import(con="/Users/shahr2/git/BIC-variants_pipeline/targets/IMPACT410_hg19/IMPACT410_hg19_baits.bed", format="bed")
opp <- genOpportunityFromGenome(BSgenome.Hsapiens.UCSC.hg19,
                                target_regions, nsamples=nrow(mut))
signatures <- signeR(M=mut, Opport=opp)
pdf(file = "../NewRuns/signatures/s_IM_GBM_111_CSF_HC_sig.pdf",
    width = 10,
    height = 2)
SignPlot(signatures$SignExposures)
dev.off()
#
vcfobj <- readVcf("../NewRuns/signatures/s_IM_GBM_63_RZ_HC.vcf", "hg19")
mut <- genCountMatrixFromVcf(BSgenome.Hsapiens.UCSC.hg19, vcfobj)

target_regions <- import(con="/Users/shahr2/git/BIC-variants_pipeline/targets/IMPACT410_hg19/IMPACT410_hg19_baits.bed", format="bed")
opp <- genOpportunityFromGenome(BSgenome.Hsapiens.UCSC.hg19,
                                target_regions, nsamples=nrow(mut))
signatures <- signeR(M=mut, Opport=opp)
pdf(file = "../NewRuns/signatures/s_IM_GBM_63_RZ_HC_sig.pdf",
    width = 10,
    height = 2)
SignPlot(signatures$SignExposures)
dev.off()
#
vcfobj <- readVcf("../NewRuns/signatures/s_IM_GBM_63_R_HC.vcf", "hg19")
mut <- genCountMatrixFromVcf(BSgenome.Hsapiens.UCSC.hg19, vcfobj)

target_regions <- import(con="/Users/shahr2/git/BIC-variants_pipeline/targets/IMPACT410_hg19/IMPACT410_hg19_baits.bed", format="bed")
opp <- genOpportunityFromGenome(BSgenome.Hsapiens.UCSC.hg19,
                                target_regions, nsamples=nrow(mut))
signatures <- signeR(M=mut, Opport=opp)
pdf(file = "../NewRuns/signatures/s_IM_GBM_63_R_HC_sig.pdf",
    width = 10,
    height = 2)
SignPlot(signatures$SignExposures)
dev.off()
#
vcfobj <- readVcf("../NewRuns/signatures/s_IM_GBM_63_CSFa_HC.vcf", "hg19")
mut <- genCountMatrixFromVcf(BSgenome.Hsapiens.UCSC.hg19, vcfobj)

target_regions <- import(con="/Users/shahr2/git/BIC-variants_pipeline/targets/IMPACT410_hg19/IMPACT410_hg19_baits.bed", format="bed")
opp <- genOpportunityFromGenome(BSgenome.Hsapiens.UCSC.hg19,
                                target_regions, nsamples=nrow(mut))
signatures <- signeR(M=mut, Opport=opp)
pdf(file = "../NewRuns/signatures/s_IM_GBM_63_CSFa_HC_sig.pdf",
    width = 10,
    height = 2)
SignPlot(signatures$SignExposures)
dev.off()
#
vcfobj <- readVcf("NewRuns/signatures/s_IM_GBM_63_CSFb_HC.vcf", "hg19")
mut <- genCountMatrixFromVcf(BSgenome.Hsapiens.UCSC.hg19, vcfobj)

target_regions <- import(con="/Users/shahr2/git/BIC-variants_pipeline/targets/IMPACT410_hg19/IMPACT410_hg19_baits.bed", format="bed")
opp <- genOpportunityFromGenome(BSgenome.Hsapiens.UCSC.hg19,
                                target_regions, nsamples=nrow(mut))
signatures <- signeR(M=mut, Opport=opp)
pdf(file = "../NewRuns/signatures/s_IM_GBM_63_CSFb_HC_sig.pdf",
    width = 10,
    height = 2)
SignPlot(signatures$SignExposures)
dev.off()
#
vcfobj <- readVcf("../NewRuns/signatures/s_IM_GBM_21_R_HC.vcf", "hg19")
mut <- genCountMatrixFromVcf(BSgenome.Hsapiens.UCSC.hg19, vcfobj)

target_regions <- import(con="/Users/shahr2/git/BIC-variants_pipeline/targets/IMPACT410_hg19/IMPACT410_hg19_baits.bed", format="bed")
opp <- genOpportunityFromGenome(BSgenome.Hsapiens.UCSC.hg19,
                                target_regions, nsamples=nrow(mut))
signatures <- signeR(M=mut, Opport=opp)
pdf(file = "../NewRuns/signatures/s_IM_GBM_21_R_HC_sig.pdf",
    width = 10,
    height = 2)
SignPlot(signatures$SignExposures)
dev.off()
#
vcfobj <- readVcf("../NewRuns/signatures/s_IM_GBM_21_CSF_HC.vcf", "hg19")
mut <- genCountMatrixFromVcf(BSgenome.Hsapiens.UCSC.hg19, vcfobj)

target_regions <- import(con="/Users/shahr2/git/BIC-variants_pipeline/targets/IMPACT410_hg19/IMPACT410_hg19_baits.bed", format="bed")
opp <- genOpportunityFromGenome(BSgenome.Hsapiens.UCSC.hg19,
                                target_regions, nsamples=nrow(mut))
signatures <- signeR(M=mut, Opport=opp)
pdf(file = "../NewRuns/signatures/s_IM_GBM_21_CSF_HC_sig.pdf",
    width = 10,
    height = 2)
SignPlot(signatures$SignExposures)
dev.off()
#
vcfobj <- readVcf("../NewRuns/signatures/s_IM_GBM_21_PZ_HC.vcf", "hg19")
mut <- genCountMatrixFromVcf(BSgenome.Hsapiens.UCSC.hg19, vcfobj)

target_regions <- import(con="/Users/shahr2/git/BIC-variants_pipeline/targets/IMPACT410_hg19/IMPACT410_hg19_baits.bed", format="bed")
opp <- genOpportunityFromGenome(BSgenome.Hsapiens.UCSC.hg19,
                                target_regions, nsamples=nrow(mut))
signatures <- signeR(M=mut, Opport=opp)
pdf(file = "../NewRuns/signatures/s_IM_GBM_21_PZ_HC_sig.pdf",
    width = 10,
    height = 2)
SignPlot(signatures$SignExposures)
dev.off()

