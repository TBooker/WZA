#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

library(vegan)
library(pegas)

## This script implements an RDA analysis on simulated data
## The input is the a BGZipped VCF file and the output is a dataframe of the SNPs
## and the number of standard deviations their Z scores from the single constained axis.

vcf_genotypes_raw <- read.vcf(args[1],
                      from =1,
                       to = 1e6)

positions <- VCFloci(args[1])

print("Extracting genotype matrix - part 1")

vcf_genotypes <- as.matrix(vcf_genotypes_raw)
colnames( vcf_genotypes ) <- positions$POS
vcf_genotypes[vcf_genotypes=="0|0"] <- 0
vcf_genotypes[vcf_genotypes=="1|0"] <- 1
vcf_genotypes[vcf_genotypes=="0|1"] <- 1
vcf_genotypes[vcf_genotypes=="1|1"] <- 2

vcf_genotypes[vcf_genotypes=="2|1"] <- 0
vcf_genotypes[vcf_genotypes=="1|2"] <- 0
vcf_genotypes[vcf_genotypes=="2|0"] <- 0
vcf_genotypes[vcf_genotypes=="0|2"] <- 0
vcf_genotypes[vcf_genotypes=="2|2"] <- 0

print("Extracting genotype matrix - part 2")

genotype_matrix <- apply(vcf_genotypes, 2, as.numeric)

envs <- read.table(args[2])

names(envs) <- c("env")

print("Building RDA model")

model <- rda(genotype_matrix ~ ., data = envs, scale = T)

print("Extracting loadings")
load.rda <- scores( model, choices = c(1),
                    display = "species")

output <- data.frame(loading = load.rda,
                     pos = row.names(load.rda))

write.csv(output, file = args[3], row.names = FALSE)
