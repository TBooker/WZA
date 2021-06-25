rm (list = ls ())

all_pine_bf <- read.table ("~/UBC/Convergence/Yeaman_Data/var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window.all_assoc_pval_bf", T, comment.char = "&")
quantile( na.omit(all_pine_bf$DD_0), 0.99)
BF_to_empiricalP <- function(x){
  rank(-1*x)/length(x)
}

all_pine_bf$DD_0_ep <- BF_to_empiricalP(all_pine_bf$DD_0)

write.table(all_pine_bf,file = '~/UBC/Convergence/Yeaman_Data/var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window.all_assoc_pval_bf_eP.tsv',  , quote=FALSE, sep='\t', row.names = F)
