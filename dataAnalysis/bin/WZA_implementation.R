## This code assumes that you have a dataframe (called GEA) 
## Each row corresponds to a single SNP
##
## The dataframe should have at least three columns:
## 
## GEA$pVal - the p-value for a particular SNP - could be spearman's rho or whatever you happen to habe
## GEA$pbar_qbar - the average minor allele frequency across populations (pbar) multiplied by the average major allele frequency  
## GEA$gene - the name of the gene or genomic region that a SNP corresponds to

## Convert one-sided p-values to Z-scores
GEA$z <- qnorm(GEA$pVal , lower.tail = F)
## Calculate the numerator of the Weighted-Z score
weiZ_num <- tapply( GEA$pbar_qbar * GEA$z, GEA$gene, sum )
## Calculate the denominator of the Weighted-Z score
weiZ_den <- sqrt(tapply( GEA$pbar_qbar^2, GEA$gene, sum ))
## Bring data together into a new DF
Z_df <- data.frame( gene = names(weiZ_num), weiZ = weiZ_num/weiZ_den)
