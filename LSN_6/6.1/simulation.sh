#!/bin/bash

# check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 [U|M|C|CHI] [Metro|Gibbs]"
    exit 1
fi

python params_manager.py "$1" "$2"
./simulationTemperatures.sh "$1" "$2"