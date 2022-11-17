
## Let's start by performing GEA on each of the simulated datasets


for map in  "cline" "trunc" #"BC_Map"
do
  echo $map
## Set the optima files
  if [[ "$map" == "BC_Map" ]]; then
      optima_file="/home/booker/work/GEA/simulations/slim_configs/BC_Map_environments.14x14.txt"
  elif [[ "$map" == "cline" ]]; then
    optima_file="/home/booker/work/GEA/simulations/slim_configs/oneDCline_environments.14x14.txt"
  elif [[ "$map" == "trunc" ]]; then
      optima_file="/home/booker/work/GEA/simulations/slim_configs/trunc_environments.14x14.txt"
  else
    exit 0
  fi

  if [[ "$map" == "BC_Map" ]]; then
      map_no=1
  elif [[ "$map" == "cline" ]]; then
    map_no=3
  elif [[ "$map" == "trunc" ]]; then
      map_no=2
  else
    exit 0
  fi


  echo $optima_fileq

  for selCoeff in 0.003 0.0136
    do
      cd s$selCoeff
      for n in 10 20 40
        do
          for i in 5 10 20
            do
#        if [[ "$map" == "BC_Map" ]]; then
#          continue
#        fi
        mkdir ${map}_output_mutMap_n${n}_i${i}
        echo s$selCoeff

        parallel -j3 "python3 /home/booker/work/GEA/simulations/bin/parseSimulations.py \
                        --trees ${map}/{}_${selCoeff}_${map_no}.12Loci.directionalSelection.trees \
                        --optima $optima_file \
                        --nPops $n \
                        --nInds $i \
                        --bayPass \
                        --directional \
                        --s $selCoeff \
                         --output ${map}_output_mutMap_n${n}_i${i}/{}_${selCoeff}_${map_no}.12Loci.directionalSelection_d40n50_mutMap \
                         --mutMap" ::: $(seq 1 30)

                                 mkdir ${map}_output_mutMap_n${n}_i${i}
                                 echo s$selCoeff

         mkdir ${map}_output_n${n}_i${i}

         parallel -j3 "python3 /home/booker/work/GEA/simulations/bin/parseSimulations.py \
                         --trees ${map}/{}_${selCoeff}_${map_no}.12Loci.directionalSelection.trees \
                         --optima $optima_file \
                         --nPops $n \
                         --nInds $i \
                         --bayPass \
                         --directional \
                         --s $selCoeff \
                          --output ${map}_output_n${n}_i${i}/{}_${selCoeff}_${map_no}.12Loci.directionalSelection_d40n50_mutMap" ::: $(seq 1 30)
        done
      done
      cd ../
    done
  done
exit 0



## d=40, n=50

for map in  "BC_Map" #"cline" "trunc"
do
  echo $map
## Set the optima files
  if [[ "$map" == "BC_Map" ]]; then
      optima_file="/home/booker/work/GEA/simulations/slim_configs/BC_Map_environments.14x14.txt"
  elif [[ "$map" == "cline" ]]; then
    optima_file="/home/booker/work/GEA/simulations/slim_configs/oneDCline_environments.14x14.txt"
  elif [[ "$map" == "trunc" ]]; then
      optima_file="/home/booker/work/GEA/simulations/slim_configs/trunc_environments.14x14.txt"
  else
    exit 0
  fi

  if [[ "$map" == "BC_Map" ]]; then
      map_no=1
  elif [[ "$map" == "cline" ]]; then
    map_no=3
  elif [[ "$map" == "trunc" ]]; then
      map_no=2
  else
    exit 0
  fi


  echo $optima_file

    for selCoeff in 0.003 0.0136
      do
#        if [[ "$map" == "BC_Map" ]]; then
#          continue
#        fi
        cd s$selCoeff
        mkdir ${map}_output_mutMap
        echo s$selCoeff

        parallel -j2 "python3 /home/booker/work/GEA/simulations/bin/parseSimulations.py \
                        --trees ${map}/{}_${selCoeff}_${map_no}.12Loci.directionalSelection.trees \
                        --optima $optima_file \
                        --nPops 40 \
                        --nInds 50 \
                        --bayPass \
                        --directional \
                        --s $selCoeff \
                         --output ${map}_output_mutMap/{}_${selCoeff}_${map_no}.12Loci.directionalSelection_d40n50_mutMap \
                         --mutMap" ::: $(seq 1 30)
        cd ../
    done
  done
exit 0



## d=40, n=50

for map in  "cline" "trunc" "BC_Map"
do
  echo $map
## Set the optima files
  if [[ "$map" == "BC_Map" ]]; then
      optima_file="/home/booker/work/GEA/simulations/slim_configs/BC_Map_environments.14x14.txt"
  elif [[ "$map" == "cline" ]]; then
    optima_file="/home/booker/work/GEA/simulations/slim_configs/oneDCline_environments.14x14.txt"
  elif [[ "$map" == "trunc" ]]; then
      optima_file="/home/booker/work/GEA/simulations/slim_configs/trunc_environments.14x14.txt"
  else
    exit 0
  fi

  if [[ "$map" == "BC_Map" ]]; then
      map_no=1
  elif [[ "$map" == "cline" ]]; then
    map_no=3
  elif [[ "$map" == "trunc" ]]; then
      map_no=2
  else
    exit 0
  fi


  echo $optima_file

  for interval in 100 1000
    do
    for selCoeff in 0.003 0.0136
      do
        if [[ "$map" == "BC_Map" ]]; then
          continue
        fi
        cd s$selCoeff
        mkdir ${map}_output_i${interval}
        echo s$selCoeff

        parallel -j2 "python3 /home/booker/work/GEA/simulations/bin/parseSimulations.py \
                        --trees ${map}/{}_${selCoeff}_${map_no}.12Loci.directionalSelection.trees \
                        --optima $optima_file \
                        --nPops 40 \
                        --nInds 50 \
                        --bayPass \
                        --directional \
                        --intervals $interval
                        --s $selCoeff \
                        --VCF \
                         --output ${map}_output_i${interval}/{}_${selCoeff}_${map_no}.12Loci.directionalSelection_d40n50_i${interval}" ::: $(seq 1 30)
        mkdir ${map}_bayPass
        mv ${map}_output/*bayPass.txt ${map}_bayPass/
        mkdir ${map}_VCF
        mv ${map}_VCF/*vcf ${map}_VCF
        cd ../
      done
    done
  done
exit 0




## d=20, 10, 5
for inds in 5 10 20
  do
  for pops in 40 20 10
    do
      for map in  "BC_Map"
#      for map in  "cline" "trunc" "BC_Map"
      do
        echo $inds $pops $map
      ## Set the optima files
        if [[ "$map" == "BC_Map" ]]; then
            optima_file="/home/booker/work/GEA/simulations/slim_configs/BC_Map_environments.14x14.txt"
        elif [[ "$map" == "cline" ]]; then
          optima_file="/home/booker/work/GEA/simulations/slim_configs/oneDCline_environments.14x14.txt"
        elif [[ "$map" == "trunc" ]]; then
            optima_file="/home/booker/work/GEA/simulations/slim_configs/trunc_environments.14x14.txt"
        else
          exit 0
        fi

        if [[ "$map" == "BC_Map" ]]; then
            map_no=1
        elif [[ "$map" == "cline" ]]; then
          map_no=3
        elif [[ "$map" == "trunc" ]]; then
            map_no=2
        else
          exit 0
        fi


        echo $optima_file

        for selCoeff in 0.003 0.0136
          do
            if [[ "$map" != "BC_Map" ]]; then
              continue
            fi
            cd s$selCoeff
            mkdir ${map}_output_n${pops}_i${inds}
            echo s$selCoeff

            parallel -j2 "python3 /home/booker/work/GEA/simulations/bin/parseSimulations.py \
                            --trees ${map}/{}_${selCoeff}_${map_no}.12Loci.directionalSelection.trees \
                            --optima $optima_file \
                            --nPops $pops \
                            --nInds $inds \
                            --bayPass \
                            --directional \
                            --s $selCoeff \
                             --output ${map}_output_n${pops}_i${inds}/{}_${selCoeff}_${map_no}.12Loci.directionalSelection_d${pops}n${inds}" ::: $(seq 1 30)
#             mkdir ${map}_output_n${pops}_bayPass
#             mv ${map}_output_n${pops}/*bayPass.txt ${map}_output_n${pops}_bayPass/

            cd ../
          done


      done
  done
done
exit 0



## d=40, n=50

for map in  "cline" "trunc" #"BC_Map"
do
  echo $map
## Set the optima files
  if [[ "$map" == "BC_Map" ]]; then
      optima_file="/home/booker/work/GEA/simulations/slim_configs/BC_Map_environments.14x14.txt"
  elif [[ "$map" == "cline" ]]; then
    optima_file="/home/booker/work/GEA/simulations/slim_configs/oneDCline_environments.14x14.txt"
  elif [[ "$map" == "trunc" ]]; then
      optima_file="/home/booker/work/GEA/simulations/slim_configs/trunc_environments.14x14.txt"
  else
    exit 0
  fi

  if [[ "$map" == "BC_Map" ]]; then
      map_no=1
  elif [[ "$map" == "cline" ]]; then
    map_no=3
  elif [[ "$map" == "trunc" ]]; then
      map_no=2
  else
    exit 0
  fi


  echo $optima_file

  for selCoeff in 0.003 0.0136
    do
      if [[ "$map" == "BC_Map" ]]; then
        continue
      fi
      cd s$selCoeff
  #    mkdir ${map}_output
      echo s$selCoeff

      parallel -j2 "python3 /home/booker/work/GEA/simulations/bin/parseSimulations.py \
                      --trees ${map}/{}_${selCoeff}_${map_no}.12Loci.directionalSelection.trees \
                      --optima $optima_file \
                      --nPops 40 \
                      --nInds 50 \
                      --bayPass \
                      --directional \
                      --s $selCoeff \
                      --VCF \
                       --output ${map}_output/{}_${selCoeff}_${map_no}.12Loci.directionalSelection_d40n50" ::: $(seq 1 30)
      mkdir ${map}_bayPass
      mv ${map}_output/*bayPass.txt ${map}_bayPass/
      mkdir ${map}_VCF
      mv ${map}_VCF/*vcf ${map}_VCF
      cd ../
    done

    cd Vs192
      mkdir ${map}_output
      echo Vs192
      parallel -j2 "python3 /home/booker/work/GEA/simulations/bin/parseSimulations.py \
                      --trees ${map}/{}_0.5_192.trees \
                      --optima $optima_file \
                      --nPops 40 \
                      --nInds 50 \
                      --bayPass \
                      --VCF \
                       --output ${map}_output/{}_0.5_192_d40n50" ::: $(seq 1 30)
      mkdir ${map}_bayPass
      mv ${map}_output/*bayPass.txt ${map}_bayPass/
#      cd ../
      mkdir ${map}_VCF
      mv ${map}_VCF/*vcf ${map}_VCF
#      cd ../
    cd ../

done



cd neutral
mkdir 2D_SteppingStone_output
echo 2D_SteppingStone
parallel -j2 "python3 /home/booker/work/GEA/simulations/bin/parseSimulations.py \
                  --trees 2D_SteppingStone/{}_0.0_1.12Loci.directionalSelection.trees \
                  --optima /home/booker/work/GEA/simulations/slim_configs/BC_Map_environments.14x14.txt \
                  --nPops 40 \
                  --nInds 50 \
                  --output 2D_SteppingStone_output/{}_0.0_1.12Loci.directionalSelection" ::: $(seq 1 30)
mkdir Island_output
echo Island
parallel -j2 "python3 /home/booker/work/GEA/simulations/bin/parseSimulations.py \
                  --trees Island/{}_0.0.12Loci.directionalSelection.trees \
                  --optima /home/booker/work/GEA/simulations/slim_configs/BC_Map_environments.14x14.txt \
                  --nPops 40 \
                  --nInds 50 \
                  --output Island_output/{}_0.0.12Loci.directionalSelection" ::: $(seq 1 30)
cd ../
