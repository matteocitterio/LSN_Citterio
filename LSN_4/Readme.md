Exercise 4.1 and 4.2 could be found in folder ./4.1/

- File `MD_MC.cpp` contains all the code necessary for computing pressure as described by the expression in the `.ipynb` file.
- File `MD_MC.cpp` has been modified to accepts command line arguments that have been used in the `equilibration_routine.sh` routine. Arguments:
   - bool `equilibration`: flag used if the program is running to thermalize the system
   - string `input`: string [Liquid| Gas | Solid]
   - int `FileNumber`: int used in `equilibration_routine.sh` to retrieve the proper set of parameters used in the thermalization.

to run the program as before, just launch `./NVE_NVT.exe 0 Solid 0`, namely set to 0 `equilibration` and `FileNumber` and choose a phase.

- File `input_manager.py` creates a list of files with the desired temperature parameters needed in the equilibration routine. It takes as command line argument a string `input` as the program above.
-  File `equilibration_routine.sh` runs the equilibration with different temperatures parameters as created in `input_manager.py`
