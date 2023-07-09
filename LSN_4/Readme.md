## How to use my code

The two exercises are both implemented in the same folder, i.e. `./4.1/`.

- To compile the program hit `make`
- File `MD_MC.cpp` contains all the code necessary for computing pressure as described by the expression in the `.ipynb` file.

**Please note** that i adapted the code so that it accepts command line inputs as this came in handy for parameters tuning during the equilibration phase.
A bash script, `equilibration_routine.sh` has been created in order to launch in sequence the program `./NVE_NVT.exe` with different input parameters. 
Arguments:
   - bool `equilibration`: flag used if the program is running to thermalize the system
   - string `input`: string [Liquid| Gas | Solid]
   - int `FileNumber`: int used in `equilibration_routine.sh` to retrieve the proper set of parameters used in the thermalization.

If no equilibration is needed, or simply to run the program as usual, just launch `./NVE_NVT.exe 0 Solid 0`. Namely set to 0 `equilibration` and `FileNumber` and choose a phase; please note that in such case you need to make sure `restart` is set to 1 in `input.in` and have a valid configuration.

- File `input_manager.py` creates a list of files with the desired temperature parameters needed in the equilibration routine. It takes as command line argument a string `input` as the program above.
-  File `equilibration_routine.sh` runs the equilibration with different temperatures parameters as created in `input_manager.py`
