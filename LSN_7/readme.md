## How to use my code

The three exercises 7.1, 7.2 and 7.3 are all implemented in the folder named `7.2`. The following guidelines apply:

- To compile the program, hit `make`
- To remove any '*.o' and '*.exe' files, hit `make clean`
- To make any `*.sh` file executable, please use `chmod +x name_of_the_.sh_file`

To run a simulation with my code make sure that all the executables are updated, then use `simulation.sh` to simulate the system.
It takes as input the Phase name, i.e one between `Gas`, `Liquid`, `Solid`. If you want to change other parameters of the simulation:

```
#!/bin/sh#!/bin/bash

vi input.PHASE_NAME
```
and change the desired parameters.
