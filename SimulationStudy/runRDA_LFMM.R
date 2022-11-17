#!/usr/bin/env Rscript
args = commandArgs(trailingOnly = TRUE)

library(vegan)
library(pegas)
library(LEA)

## This script implements an RDA analysis on simulated data
## The input is the BayPass input file and the output is a dataframe of the SNPs
## and the number of standard deviations their Z scores from the single constrained axis.

data <- as.matrix( read.table(args[1])  ) 

nHaps = data[1,1]+data[1,2]
select_cols = seq(1,ncol(data),2)
data_f = data[,select_cols]/nHaps

allele_freqs <- t( data_f )

envs_raw <- c(t(read.table(args[2]) ))
envs <- data.frame(env = envs_raw)

model <- rda(allele_freqs ~ ., data = envs, scale = T)

print("Extracting loadings")

load.rda <- scores( model, choices = c(1),
                    display = "species")

#output <- data.frame(loading = load.rda,
#                     pos = row.names(load.rda))

#write.csv(output, file = "~/work/GEA/simulations/SimulationStudy/TEMP/XXX_RDA.txt", row.names = FALSE)


#install.packages("devtools")
#devtools::install_github("bcm-uga/LEA")


# Fitting an LFMM with K = 3 factors
mod <- lfmm2(input = allele_freqs, env = envs, K = 3)

pv <- lfmm2.test(object = mod,
                input = allele_freqs,
                env = envs,
                linear = TRUE)


library(ggplot2)

gea <- read.csv(args[3])

gea$LEA <- pv$pvalues
gea$RDA <- load.rda

output_name = paste0(strsplit(args[3],".csv")[[1]],"_singleSNPtests.csv")

write.csv(gea, 
          file = output_name, 
          row.names = FALSE)

