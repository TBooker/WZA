#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

library(vegan)
library(pegas)

## This script implements an RDA analysis on simulated data
## The input is the a BGZipped VCF file and the output is a dataframe of the SNPs
## and the number of standard deviations their Z scores from the single constained axis.
data <- as.matrix( read.table(args[1]) ) 

nPops <- dim(data)[2]/2

extract_pops <- seq(1,nPops*2,2)
allele_freqs <- t( data[,extract_pops]/100 )

envs_raw <- c(t(read.table(args[2])))
envs <- data.frame(env = envs_raw)

model <- rda(allele_freqs ~ ., data = envs, scale = T)

print("Extracting loadings")
load.rda <- scores( model, choices = c(1),
                    display = "species")

output <- data.frame(loading = load.rda,
                     pos = row.names(load.rda))

write.csv(output, file = args[3], row.names = FALSE)
