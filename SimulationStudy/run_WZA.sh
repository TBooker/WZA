
# # Merge BayPass results with the actual results
#
# echo "BC map - mutMap"
# parallel "python3 ../bin/x_y_cbind.py --x s{2}/BC_Map_output_mutMap_n{4}_i{3}/{1}_{2}_1.12Loci.directionalSelection_d40n50_mutMap_singleSNPtests.csv \
#  --y s{2}/bayPass_results/BC_Map_output_mutMap_n{4}_i{3}/{1}_{2}_1.12Loci.directionalSelection_d40n50_mutMap.bayPass_summary_betai_reg.out \
#  --output s{2}/BC_Map_output_mutMap_n{4}_i{3}/{1}_{2}_1.12Loci.directionalSelection_d{4}n{3}_mutMap_singleSNPtests_wBP.csv" ::: $(seq 1 30) ::: 0.003 0.0136 ::: 5 10 20 ::: 10 20 40
#
#  echo "cline - mutMap"
# parallel "python3 ../bin/x_y_cbind.py --x s{2}/cline_output_mutMap_n{4}_i{3}/{1}_{2}_3.12Loci.directionalSelection_d40n50_mutMap_singleSNPtests.csv \
#   --y s{2}/bayPass_results/cline_output_mutMap_n{4}_i{3}/{1}_{2}_3.12Loci.directionalSelection_d40n50_mutMap.bayPass_summary_betai_reg.out \
#   --output s{2}/cline_output_mutMap_n{4}_i{3}/{1}_{2}_3.12Loci.directionalSelection_d{4}n{3}_mutMap_singleSNPtests_wBP.csv" ::: $(seq 1 30) ::: 0.003 0.0136 ::: 5 10 20 ::: 10 20 40
#
#   echo "trunc - mutMap"
# parallel "python3 ../bin/x_y_cbind.py --x s{2}/trunc_output_mutMap_n{4}_i{3}/{1}_{2}_2.12Loci.directionalSelection_d40n50_mutMap_singleSNPtests.csv \
#    --y s{2}/bayPass_results/trunc_output_mutMap_n{4}_i{3}/{1}_{2}_2.12Loci.directionalSelection_d40n50_mutMap.bayPass_summary_betai_reg.out \
#    --output s{2}/trunc_output_mutMap_n{4}_i{3}/{1}_{2}_2.12Loci.directionalSelection_d{4}n{3}_mutMap_singleSNPtests_wBP.csv" ::: $(seq 1 30) ::: 0.003 0.0136 ::: 5 10 20 ::: 10 20 40
#
#
#    echo "BC map - uniform"
# parallel "python3 ../bin/x_y_cbind.py --x s{2}/BC_Map_output_n{4}_i{3}/{1}_{2}_1.12Loci.directionalSelection_d40n50_mutMap_singleSNPtests.csv \
#   --y s{2}/bayPass_results/BC_Map_output_n{4}_i{3}/{1}_{2}_1.12Loci.directionalSelection_d40n50_mutMap.bayPass_summary_betai_reg.out \
#   --output s{2}/BC_Map_output_n{4}_i{3}/{1}_{2}_1.12Loci.directionalSelection_d{4}n{3}_singleSNPtests_wBP.csv" ::: $(seq 1 30) ::: 0.003 0.0136 ::: 5 10 20 ::: 10 20 40
#
#   echo "cline - uniform"
# parallel "python3 ../bin/x_y_cbind.py --x s{2}/cline_output_n{4}_i{3}/{1}_{2}_3.12Loci.directionalSelection_d40n50_mutMap_singleSNPtests.csv \
#   --y s{2}/bayPass_results/cline_output_n{4}_i{3}/{1}_{2}_3.12Loci.directionalSelection_d40n50_mutMap.bayPass_summary_betai_reg.out \
#   --output s{2}/cline_output_n{4}_i{3}/{1}_{2}_3.12Loci.directionalSelection_d{4}n{3}_singleSNPtests_wBP.csv" ::: $(seq 1 30) ::: 0.003 0.0136 ::: 5 10 20 ::: 10 20 40
#
#   echo "trunc - uniform"
# parallel "python3 ../bin/x_y_cbind.py --x s{2}/trunc_output_n{4}_i{3}/{1}_{2}_2.12Loci.directionalSelection_d40n50_mutMap_singleSNPtests.csv \
#   --y s{2}/bayPass_results/trunc_output_n{4}_i{3}/{1}_{2}_2.12Loci.directionalSelection_d40n50_mutMap.bayPass_summary_betai_reg.out \
#   --output s{2}/trunc_output_n{4}_i{3}/{1}_{2}_2.12Loci.directionalSelection_d{4}n{3}_singleSNPtests_wBP.csv" ::: $(seq 1 30) ::: 0.003 0.0136 ::: 5 10 20 ::: 10 20 40
#



## Run the WZA on each GEA method
for s in s0.003 s0.0136
do
  cd $s
  for f in BC_Map cline trunc
    do
      for d in 10 20 40
        do
          for n in 5 10 20
            do
              for i in $(ls ${f}_output_n${d}_i${n}/*wBP.csv)
                do
                  echo $i
                  id=$(echo $i | cut -d"/" -f2 | cut -d"_" -f1)
                  echo $id
                  parallel -j2 "python3 ../../bin/general_WZA_script.py \
                                --correlations $i \
                                --summary_stat  {} \
                                --window gene \
                                --output ${f}_output_n${d}_i${n}/${id}_WZA_{}.csv \
                                --MAF maf \
                                --retain LA  \
                                --sep ',' " ::: geno_k_tau_p_value LEA &
                  parallel -j2 "python3 ../../bin/general_WZA_script.py \
                                --correlations $i \
                                --summary_stat  {} \
                                --window gene \
                                --output ${f}_output_n${d}_i${n}/${id}_WZA_{}.csv \
                                --MAF maf \
                                --retain LA  \
                                --large_i_small_p \
                                --sep ',' " ::: RDA 'BF(dB)'
                  done
                  for i in $(ls ${f}_output_mutMap_n${d}_i${n}/*wBP.csv)
                    do
                      echo $i
                      id=$(echo $i | cut -d"/" -f2 | cut -d"_" -f1)
                      echo $id
                      parallel -j2 "python3 ../../bin/general_WZA_script.py \
                                    --correlations $i \
                                    --summary_stat  {} \
                                    --window gene \
                                    --output ${f}_output_mutMap_n${d}_i${n}/${id}_WZA_{}.csv \
                                    --MAF maf \
                                    --retain LA  \
                                    --sep ',' " ::: geno_k_tau_p_value LEA
                      parallel -j2 "python3 ../../bin/general_WZA_script.py \
                                    --correlations $i \
                                    --summary_stat  {} \
                                    --window gene \
                                    --output ${f}_output_mutMap_n${d}_i${n}/${id}_WZA_{}.csv \
                                    --MAF maf \
                                    --retain LA  \
                                    --large_i_small_p \
                                    --sep ',' " ::: RDA 'BF(dB)'
                      done
                done
           done
    done
    cd ../
done



        #
        # mkdir ${f}_output_n${n}_WZA
        #  for i in $(ls ${f}_output_n${n}/*singleSNPtests.csv)
        #   do
        #       id=$(echo $i | cut -d"/" -f2 | cut -d"_" -f1)
        #       echo $id
        #       # parallel -j2 "python3 ../../bin/general_WZA_script.py \
        #       #               --correlations $i \
        #       #               --summary_stat  {} \
        #       #               --window gene \
        #       #               --output ${f}_output_n${n}_WZA/${id}_WZA_{}_n{$n}.csv \
        #       #               --MAF maf \
        #       #               --retain LA  \
        #       #               --sep ',' " ::: geno_k_tau_p_value LEA
        #       parallel -j2 "python3 ../../bin/general_WZA_script.py \
        #                     --correlations $i \
        #                     --summary_stat  {} \
        #                     --window gene \
        #                     --output ${f}_output_n${n}_WZA/${id}_WZA_{}_n{$n}.csv \
        #                     --MAF maf \
        #                     --retain LA  \
        #                     --large_i_small_p \
        #                     --sep ',' " ::: RDA
        #   done

      # for i in $(ls ${f}_output/*wBP.csv)
      #   do
      #     echo $i
      #     id=$(echo $i | cut -d"/" -f2 | cut -d"_" -f1)
      #     echo $id
      #     parallel -j2 "python3 ../../bin/general_WZA_script.py \
      #                   --correlations $i \
      #                   --summary_stat  {} \
      #                   --window gene \
      #                   --output ${f}_output_WZA/${id}_WZA_{}.csv \
      #                   --MAF maf \
      #                   --retain LA  \
      #                   --sep ',' " ::: geno_k_tau_p_value LEA RDA 'BF(dB)'
      #     parallel -j2 "python3 ../../bin/general_WZA_script.py \
      #                   --correlations $i \
      #                   --summary_stat  {} \
      #                   --window gene \
      #                   --output ${f}_output_WZA/${id}_WZA_{}.csv \
      #                   --MAF maf \
      #                   --retain LA  \
      #                   --large_i_small_p \
      #                   --sep ',' " ::: RDA 'BF(dB)'
      #     done
