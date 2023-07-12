#!/bin/bash

mkdir -p output/CUPID_sens/frequentist_baseline_simple_shape/
mkdir -p output/CUPID_sens/frequentist_baseline_full_shape/
mkdir -p output/CUPID_sens/frequentist_optimist_simple_shape/
mkdir -p output/CUPID_sens/frequentist_optimist_full_shape/
mkdir -p output/CUPID_sens/frequentist_reach_simple_shape/
mkdir -p output/CUPID_sens/frequentist_reach_full_shape/


./bin/CUPID_Sens -m F -n 100000 -s 200 -b 1e-4 -f 0 -e frequentist_baseline_simple_shape 2&>1 | tee log/frequentist_baseline_simple_shape.log &
./bin/CUPID_Sens -m F -n 100000 -s 200 -b 1e-4 -f 1 -e frequentist_baseline_full_shape > log/frequentist_baseline_full_shape.log &
./bin/CUPID_Sens -m F -n 100000 -s 200 -b 6e-5 -f 0 -e frequentist_optimist_simple_shape > log/frequentist_optimist_simple_shape.log &
./bin/CUPID_Sens -m F -n 100000 -s 200 -b 6e-5 -f 1 -e frequentist_optimist_full_shape > log/frequentist_optimist_full_shape.log &
./bin/CUPID_Sens -m F -n 100000 -s 200 -b 2e-5 -f 0 -e frequentist_reach_simple_shape > log/frequentist_reach_simple_shape.log &
./bin/CUPID_Sens -m F -n 100000 -s 200 -b 2e-5 -f 1 -e frequentist_reach_full_shape > log/frequentist_reach_full_simple_shape.log &
