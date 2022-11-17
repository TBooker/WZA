




cd s0.003

## Perform WZA on each of the maps
  parallel "python3 ../../bin/WZA.py --csv BC_Map_output_n{1}_i{2}/ --output  BC_Map_output_n{1}_i{2}.wza.csv --directional" ::: 10 20 40 ::: 5 10 20

## Summarise the results using the AUC
  parallel "python3 ../../bin/summariseWZA_AUC.py --wza BC_Map_output_n{1}_i{2}.wza.csv  --LA 0.005 --output BC_Map_output_n{1}_i{2}.wza" ::: 10 20 40 ::: 5 10 20

cd ../


cd s0.0136

## Perform WZA on each of the maps
  parallel "python3 ../../bin/WZA.py --csv BC_Map_output_n{1}_i{2}/ --output  BC_Map_output_n{1}_i{2}.wza.csv --directional" ::: 10 20 40 ::: 5 10 20

## Summarise the results using the AUC
  parallel "python3 ../../bin/summariseWZA_AUC.py --wza BC_Map_output_n{1}_i{2}.wza.csv  --LA 0.005 --output BC_Map_output_n{1}_i{2}.wza" ::: 10 20 40 ::: 5 10 20

cd ../



cd s0.003

## Perform WZA on each of the maps
  parallel "python3 ../../bin/WZA.py --csv {1}_output_n{2}/ --output {1}_n{2}.wza.csv --directional" ::: BC_Map cline trunc ::: 10 20

## Summarise the results using the AUC
  parallel "python3 ../../bin/summariseWZA_AUC.py --wza {1}_n{2}.wza.csv  --LA 0.005 --output {1}_n{2}.wza" ::: BC_Map cline trunc ::: 10 20

cd ../


cd s0.0136

## Perform WZA on each of the maps
parallel "python3 ../../bin/WZA.py --csv {1}_output_n{2}/ --output {1}_n{2}.wza.csv --directional" ::: BC_Map cline trunc ::: 10 20

## Summarise the results using the AUC
parallel "python3 ../../bin/summariseWZA_AUC.py --wza {1}_n{2}.wza.csv  --LA 0.005 --output {1}_n{2}.wza" ::: BC_Map cline trunc ::: 10 20

cd ../
