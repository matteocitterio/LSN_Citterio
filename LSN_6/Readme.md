## How to use my code

- To compile the code: `make`
- To make the `*.sh` files executables, hit `chmod +x NAME_OF_THE_.sh_FILE`
- To clean the folder from '*.o' and '*.exe' files, hit `make clean`

To launch a simulation, just run `./simulation.sh`. This simple bash script contains:

``` bash
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
```

Where `parameters_manager.py` takes care of changing the input parameters in order to simulate the system at different temperatures (and external field in the magnitude case). 
Likewise, `simulationTemperatures.sh` is a simple bash routine that runs the simulation in a for loop accepting the input files created with the previous python script.
The program accepts as command line inputs the quantity you want to simulate ([U|M|C|CHI]) and the algorithm you want to use ([Metro|Gibbs]) and builds a folder with all the inputs given and the results obtained.
