cd s0.003

## Perform WZA on each of the maps
  parallel "python3 ../../bin/WZA.py --csv {}_output/ --output {}.wza.csv --bayPass {}_bayPass_analysis/ --rda {}_RDA/ --directional" ::: BC_Map cline trunc

## Summarise the results using the AUC
  parallel "python3 ../../bin/summariseWZA_AUC.py --wza {}.wza.csv --BayPass --RDA  --LA 0.005 --output {}.wza" ::: BC_Map cline trunc

cd ../

cd s0.0136

## Perform WZA on each of the maps
  parallel "python3 ../../bin/WZA.py --csv {}_output/ --output {}.wza.csv --bayPass {}_bayPass_analysis/ --rda {}_RDA/ --directional" ::: BC_Map cline trunc

## Summarise the results using the AUC
  parallel "python3 ../../bin/summariseWZA_AUC.py --wza {}.wza.csv --BayPass --RDA  --LA 0.005 --output {}.wza" ::: BC_Map cline trunc

cd ../
