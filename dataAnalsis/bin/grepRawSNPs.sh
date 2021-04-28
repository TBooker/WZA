while read t; do
	head -n1 var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window.all_assoc_pval_eP.tsv > $t.tsv

	grep $t var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window.all_assoc_pval_eP.tsv >> $t.tsv

  echo "$t"
done <ScaffoldsToFocusOn.txt


#../../../../Convergence/Yeaman_Data/var_out_GATK3_allhet_pine688_ALL.summary.ALL.annots.sorted.GOOD.window.all_assoc_pval_bf_eP.tsv
