#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 [U|M|C|CHI] [Metro|Gibbs]"
    exit 1
fi

# Execute the Python script to manage parameters
python params_manager.py "$1" "$2"

# Execute the simulation script with specified parameters
./simulationTemperatures.sh "$1" "$2"