#!/bin/sh#!/bin/bash

# set the number of cycles
num_cycles=$(ls Gas/inputs/ | wc -l)
num_cycles=$((num_cycles-1))

echo $num_cycles

# loop over the cycles
for ((i=1; i<=$num_cycles; i++))
do
    # run the program with the current cycle variable
    ./NVE_NVT.exe 1 Gas $i
done

