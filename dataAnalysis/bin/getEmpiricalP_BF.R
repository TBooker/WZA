rm (list = ls ())

# Read in the SNP table  (of Bayes factors)

#all_pine_bf <- read.table ("../Yeaman_Data/var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window.all_assoc_pval_bf", T, comment.char = "&")

## This function converts Bayes factors to empirical Ps. Takes ranks, reverses them and converts to empirical p

BF_to_empiricalP <- function(x){
  rank(-1*x)/length(x)
}

## Add the empirical ps to the data frame
#all_pine_bf$DD_0_ep <- BF_to_empiricalP(all_pine_bf$DD_0)

## Write the results to a file

#write.table(all_pine_bf,file = 'var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window.all_assoc_pval_bf_eP.tsv',  quote=FALSE, sep='\t', row.names = F)



# Read in the SNP table  (of p-values)

all_pine_bf <- read.table ("../../Yeaman_Data/doi_10.5061_dryad.0t407__v1/yeaman_pine_spruce_convergence_archive2/var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window.all_assoc_pval_raw", T, comment.char = "&")
all_pine_bf <- read.table ("../../Yeaman_Data/doi_10.5061_dryad.0t407__v1/temp", T, comment.char = "&")

quantile( na.omit(all_pine_bf$DD_0), 0.99)

## This function converts Bayes factors to empirical Ps. Takes ranks, reverses them and converts to empirical p

pVal_to_empiricalP <- function(x){
  rank(x)/length(x)
}

## Add the empirical ps to the data frame
all_pine_bf$DD_0_ep <- pVal_to_empiricalP(all_pine_bf$DD_0)

## Write the results to a file

write.table(all_pine_bf,file = 'var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window.all_assoc_pval_eP.tsv',  quote=FALSE, sep='\t', row.names = F)


