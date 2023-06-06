#!/bin/bash



./bin/CUPID_Sens -m B -n 10000 -s 1 -b 1e-4 -f 0 -e bayesian_baseline_simple_shape > log/bayesian_baseline_simple_shape.log &
./bin/CUPID_Sens -m B -n 10000 -s 1 -b 1e-4 -f 1 -e bayesian_baseline_full_shape > log/bayesian_baseline_full_shape.log &
./bin/CUPID_Sens -m B -n 10000 -s 1 -b 6e-5 -f 0 -e bayesian_optimist_simple_shape > log/bayesian_optimist_simple_shape.log &
./bin/CUPID_Sens -m B -n 10000 -s 1 -b 6e-5 -f 1 -e bayesian_optimist_full_shape > log/bayesian_optimist_full_shape.log &
./bin/CUPID_Sens -m B -n 10000 -s 1 -b 2e-5 -f 0 -e bayesian_reach_simple_shape > log/bayesian_reach_simple_shape.log &
./bin/CUPID_Sens -m B -n 10000 -s 1 -b 2e-5 -f 1 -e bayesian_reach_full_shape > log/bayesian_reach_full_simple_shape.log
