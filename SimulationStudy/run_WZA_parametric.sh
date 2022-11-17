## Run the WZA on each GEA method
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
            parallel -j2 "python3 ../../bin/general_WZA_script_dev.py \
                          --correlations $i \
                          --summary_stat  {} \
                          --window gene \
                          --output ${f}_output_n${d}_i${n}_parametric/${id}_WZA_{}_unweighted.csv \
                          --MAF maf \
                          --retain LA  \
                          --sep ',' " ::: geno_k_tau_p_value
            done
          done
        done
        cd ../

      done
