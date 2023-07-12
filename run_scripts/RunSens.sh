#!/bin/bash


mkdir -p output/CUPID_sens/bayesian_$1/
mkdir -p output/CUPID_sens/bayesian_$1/limits/

#./bin/CUPID_Sens -m B -n 10000 -s 1 -b 1e-4 -f 0 -e bayesian_baseline_simple_shape > log/bayesian_baseline_simple_shape.log &
#./bin/CUPID_Sens -m B -n 10000 -s 1 -b 1e-4 -f 1 -e bayesian_baseline_full_shape -F 0 -L 100 -G 0
#./bin/CUPID_Sens -m B -n 10000 -s 1 -b 6e-5 -f 0 -e bayesian_optimist_simple_shape > log/bayesian_optimist_simple_shape.log &
#./bin/CUPID_Sens -m B -n 10000 -s 1 -b 6e-5 -f 1 -e bayesian_optimist_full_shape > log/bayesian_optimist_full_shape.log &
#./bin/CUPID_Sens -m B -n 10000 -s 1 -b 2e-5 -f 0 -e bayesian_reach_simple_shape > log/bayesian_reach_simple_shape.log &
#./bin/CUPID_Sens -m B -n 10000 -s 1 -b 2e-5 -f 1 -e bayesian_reach_full_shape > log/bayesian_reach_full_shape.log &

export step=1000

for i in {0..9}
do

    export first=$((x=$step,y=$i,x*y))
    export last=$((x=$first,y=$step,x+y))
    echo "-m B -n 10000 -s 1 -b ${2} -f 1 -e bayesian_${1} -F ${first} -L ${last} -G ${i}"
    qsub Sensitivity.pbs -v "CONFIG=-m B -n 10000 -s 1 -b ${2} -f ${3} -e bayesian_${1} -F ${first} -L ${last} -G ${i}" 

done
