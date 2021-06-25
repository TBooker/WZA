head -n1 ../../../../Convergence/Yeaman_Data/var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window.all_assoc_pval_bf_eP.tsv > $1.tsv

grep $1 ../../../../Convergence/Yeaman_Data/var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window.all_assoc_pval_bf_eP.tsv >> $1.tsv

#../../../../Convergence/Yeaman_Data/var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window.all_assoc_pval_bf_eP.tsv
