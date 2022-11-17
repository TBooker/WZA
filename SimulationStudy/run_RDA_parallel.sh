for s in 0.003 0.0136
  do
    cd s${s}
## BC Map
# n40 - 50 individuals
#    parallel "Rscript ../runRDA_LFMM.R BC_Map_bayPass/{1}_${s}_1.12Loci.directionalSelection_d40n50.bayPass.txt BC_Map_bayPass/{1}_${s}_1.12Loci.directionalSelection_d40n50.bayPass.pc1 BC_Map_output/{1}_${s}_1.12Loci.directionalSelection_d40n50.csv" ::: $(seq 1 30)
# n10 n20 - 50 individuals
#    parallel "Rscript ../runRDA_LFMM.R BC_Map_output_n{2}_bayPass/{1}_${s}_1.12Loci.directionalSelection_d{2}n50.bayPass.txt BC_Map_output_n{2}/{1}_${s}_1.12Loci.directionalSelection_d{2}n50.bayPass.pc1 BC_Map_output_n{2}/{1}_${s}_1.12Loci.directionalSelection_d{2}n50.csv" ::: $(seq 1 30) ::: 10 20
# n10 n20 n40 - 5 10 20 individuals
parallel "Rscript ../runRDA_LFMM.R BC_Map_output_mutMap_n{2}_i{3}/{1}_${s}_1.12Loci.directionalSelection_d40n50_mutMap.bayPass.txt BC_Map_output_mutMap_n{2}_i{3}/{1}_${s}_1.12Loci.directionalSelection_d40n50_mutMap.bayPass.pc1 BC_Map_output_mutMap_n{2}_i{3}/{1}_${s}_1.12Loci.directionalSelection_d40n50_mutMap.csv" ::: $(seq 1 30) ::: 10 20 40 ::: 5 10 20
parallel "Rscript ../runRDA_LFMM.R cline_output_mutMap_n{2}_i{3}/{1}_${s}_3.12Loci.directionalSelection_d40n50_mutMap.bayPass.txt cline_output_mutMap_n{2}_i{3}/{1}_${s}_3.12Loci.directionalSelection_d40n50_mutMap.bayPass.pc1 cline_output_mutMap_n{2}_i{3}/{1}_${s}_3.12Loci.directionalSelection_d40n50_mutMap.csv" ::: $(seq 1 30) ::: 10 20 40 ::: 5 10 20
parallel "Rscript ../runRDA_LFMM.R trunc_output_mutMap_n{2}_i{3}/{1}_${s}_2.12Loci.directionalSelection_d40n50_mutMap.bayPass.txt trunc_output_mutMap_n{2}_i{3}/{1}_${s}_2.12Loci.directionalSelection_d40n50_mutMap.bayPass.pc1 trunc_output_mutMap_n{2}_i{3}/{1}_${s}_2.12Loci.directionalSelection_d40n50_mutMap.csv" ::: $(seq 1 30) ::: 10 20 40 ::: 5 10 20

parallel "Rscript ../runRDA_LFMM.R BC_Map_output_n{2}_i{3}/{1}_${s}_1.12Loci.directionalSelection_d40n50_mutMap.bayPass.txt BC_Map_output_n{2}_i{3}/{1}_${s}_1.12Loci.directionalSelection_d40n50_mutMap.bayPass.pc1 BC_Map_output_n{2}_i{3}/{1}_${s}_1.12Loci.directionalSelection_d40n50_mutMap.csv" ::: $(seq 1 30) ::: 10 20 40 ::: 5 10 20
parallel "Rscript ../runRDA_LFMM.R cline_output_n{2}_i{3}/{1}_${s}_3.12Loci.directionalSelection_d40n50_mutMap.bayPass.txt cline_output_n{2}_i{3}/{1}_${s}_3.12Loci.directionalSelection_d40n50_mutMap.bayPass.pc1 cline_output_n{2}_i{3}/{1}_${s}_3.12Loci.directionalSelection_d40n50_mutMap.csv" ::: $(seq 1 30) ::: 10 20 40 ::: 5 10 20
parallel "Rscript ../runRDA_LFMM.R trunc_output_n{2}_i{3}/{1}_${s}_2.12Loci.directionalSelection_d40n50_mutMap.bayPass.txt trunc_output_n{2}_i{3}/{1}_${s}_2.12Loci.directionalSelection_d40n50_mutMap.bayPass.pc1 trunc_output_n{2}_i{3}/{1}_${s}_2.12Loci.directionalSelection_d40n50_mutMap.csv" ::: $(seq 1 30) ::: 10 20 40 ::: 5 10 20
cd ../
done

exit 0




## Cline
# n40 - 50 individuals
#    parallel "Rscript ../runRDA_LFMM.R cline_bayPass/{1}_${s}_3.12Loci.directionalSelection_d40n50.bayPass.txt cline_bayPass/{1}_${s}_3.12Loci.directionalSelection_d40n50.bayPass.pc1 cline_output/{1}_${s}_3.12Loci.directionalSelection_d40n50.csv" ::: $(seq 1 30)
# n10 n20 - 50 individuals
    parallel "Rscript ../runRDA_LFMM.R cline_output_n{2}_bayPass/{1}_${s}_3.12Loci.directionalSelection_d{2}n50.bayPass.txt cline_output_n{2}/{1}_${s}_3.12Loci.directionalSelection_d{2}n50.bayPass.pc1 cline_output_n{2}/{1}_${s}_3.12Loci.directionalSelection_d{2}n50.csv" ::: $(seq 1 30) ::: 10 20 40 ::: 5 10 20

## Trunc
# n40 - 50 individuals
#    parallel "Rscript ../runRDA_LFMM.R trunc_bayPass/{1}_${s}_2.12Loci.directionalSelection_d40n50.bayPass.txt trunc_bayPass/{1}_${s}_2.12Loci.directionalSelection_d40n50.bayPass.pc1 trunc_output/{1}_${s}_2.12Loci.directionalSelection_d40n50.csv" ::: $(seq 1 30)
# n10 n20 - 50 individuals
  parallel "Rscript ../runRDA_LFMM.R trunc_output_n{2}_i{3}/{1}_${s}_2.12Loci.directionalSelection_d{2}n50.bayPass.txt trunc_output_n{2}/{1}_${s}_2.12Loci.directionalSelection_d{2}n50.bayPass.pc1 trunc_output_n{2}/{1}_${s}_2.12Loci.directionalSelection_d{2}n50.csv" ::: $(seq 1 30) ::: 10 20 40 ::: 5 10 20
  parallel "Rscript ../runRDA_LFMM.R trunc_output_n{2}_i{3}/{1}_${s}_2.12Loci.directionalSelection_d{2}n50.bayPass.txt trunc_output_n{2}/{1}_${s}_2.12Loci.directionalSelection_d{2}n50.bayPass.pc1 trunc_output_n{2}/{1}_${s}_2.12Loci.directionalSelection_d{2}n50.csv" ::: $(seq 1 30) ::: 10 20 40 ::: 5 10 20

#  parallel "python3 ../../bin/x_y_cbind.py --x BC_Map_output/{}_${s}_1.12Loci.directionalSelection_d40n50_singleSNPtests.csv --y BC_Map_bayPass_analysis/{}_${s}_1.12Loci.directionalSelection_d40n50.bayPass_summary_betai_reg.out --output BC_Map_output/{}_${s}_1.12Loci.directionalSelection_d40n50_singleSNPtests_wBP.csv" ::: $(seq 1 30)
#  parallel "python3 ../../bin/x_y_cbind.py --x cline_output/{}_${s}_3.12Loci.directionalSelection_d40n50_singleSNPtests.csv --y cline_bayPass_analysis/{}_${s}_3.12Loci.directionalSelection_d40n50.bayPass_summary_betai_reg.out --output cline_output/{}_${s}_3.12Loci.directionalSelection_d40n50_singleSNPtests_wBP.csv" ::: $(seq 1 30)
#  parallel "python3 ../../bin/x_y_cbind.py --x trunc_output/{}_${s}_2.12Loci.directionalSelection_d40n50_singleSNPtests.csv --y trunc_bayPass_analysis/{}_${s}_2.12Loci.directionalSelection_d40n50.bayPass_summary_betai_reg.out --output trunc_output/{}_${s}_2.12Loci.directionalSelection_d40n50_singleSNPtests_wBP.csv" ::: $(seq 1 30)
    cd ../
done

exit 0
cd s0.0136

parallel "python3 ../../bin/WZA.py --csv BC_Map_output_n{1}_i{2}/ --output  BC_Map_output_n{1}_i{2}.wza.csv --directional --rda BC_Map_output_n{1}_i{2}_RDA/" ::: 10 20 40 ::: 5 10 20

parallel "python3 ../../bin/summariseWZA_AUC.py --wza BC_Map_output_n{1}_i{2}.wza.csv  --LA 0.005 --output BC_Map_output_n{1}_i{2}.wza" ::: 10 20 40 ::: 5 10 20

cd ../

exit 0

cd s0.003
#run_RDA_parallel.sh
parallel "mkdir BC_Map_output_n{1}_i{2}_RDA" ::: 10 20 40 ::: 5 10 20

parallel "Rscript ../runRDA.R BC_Map_output_n{2}_i{3}/{1}_0.003_1.12Loci.directionalSelection_d{2}n{3}.bayPass.txt BC_Map_output_n{2}_i{3}/{1}_0.003_1.12Loci.directionalSelection_d20n5.bayPass.pc1 BC_Map_output_n{2}_i{3}_RDA/{1}_0.003_1.12Loci.directionalSelection_d{2}n{3}.rda.txt" ::: $(seq 1 30) ::: 10 20 40 ::: 5 10 20

parallel "python3 ../../bin/WZA.py --csv BC_Map_output_n{1}_i{2}/ --output  BC_Map_output_n{1}_i{2}.wza.csv --directional --rda BC_Map_output_n{1}_i{2}_RDA/" ::: 10 20 40 ::: 5 10 20

cd ../

exit 0

cd s0.0136
#run_RDA_parallel.sh
parallel "mkdir BC_Map_output_n{1}_i{2}_RDA" ::: 10 20 40 ::: 5 10 20

parallel "Rscript ../runRDA.R BC_Map_output_n{2}_i{3}/{1}_0.0136_1.12Loci.directionalSelection_d{2}n{3}.bayPass.txt BC_Map_output_n{2}_i{3}/{1}_0.0136_1.12Loci.directionalSelection_d{2}n{3}.bayPass.pc1 BC_Map_output_n{2}_i{3}_RDA/{1}_0.0136_1.12Loci.directionalSelection_d{2}n{3}.rda.txt" ::: $(seq 1 30) ::: 10 20 40 ::: 5 10 20

#parallel "python3 ../../bin/WZA.py --csv BC_Map_output_n{1}_i{2}/ --output  BC_Map_output_n{1}_i{2}.wza.csv --directional --rda BC_Map_output_n{1}_i{2}_RDA/" ::: 10 20 40 ::: 5 10 20
cd ../
