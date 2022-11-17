for s in s0.0136 s0.003
do
  cd $s

parallel "python3 ../../bin/summariseWZA_AUC_singleFile.py \
        --result {1}_output_n{2}_i20_SNP_number/ \
        --singleSNP {1}_output_n{2}_i20/ \
        --LA 0.005 \
        --output {1}_output_n{2}_i20_SNP_number.csv" ::: BC_Map cline trunc ::: 10 20 40
cd ../
done
exit 0

## Run the WZA on windows of varying numbers of SNPs
for s in s0.0136 s0.003
do
  cd $s
  for f in BC_Map cline trunc
    do
      for d in 10 20 40
        do
        n=20
        for i in $(ls ${f}_output_n${d}_i${n}/*wBP.csv)
          do
            echo $i
            id=$(echo $i | cut -d"/" -f2 | cut -d"_" -f1)
            echo $id
            parallel "python3 ../../bin/general_WZA_script_dev.py \
            --correlations $i \
            --summary_stat geno_k_tau_p_value \
            --window gene \
            --output ${f}_output_n${d}_i${n}_SNP_number/${id}_WZA_{}snps.csv \
            --sep ',' \
            --MAF maf \
            --snp_window {} \
            --retain LA" ::: 2 4 8 10
            done
          done
        done
        cd ../

      done
