import os
import sys
import numpy as np

'''
This program builds a set of `input.in` files so that we can cycle over the parameters and extract the best parameters needed for the thermalization.
Input: 
- Please provide the name of the phase, one choice among: ['Liquid', 'Solid', 'Gas'], this will grant access to the specific folder.
'''

acceptedArgs = ['Liquid', 'Solid', 'Gas']                                               # Admissable command line inputs
phase_name = sys.argv[1]                                                                # Retrieve the phase name from the command line argument

# Check if the provided phase name is one of the accepted values
if phase_name in acceptedArgs:
    dir_name ='./'+ phase_name +"/inputs"
else: 
    # Raise an error if the phase name is not accepted
    raise KeyError('Arg not accepted, list of accepted args:', acceptedArgs)

# Check if the directory exists, if not, create it
if not os.path.exists(dir_name):
    print('Making a new directory')
    os.makedirs(dir_name)


# Open the file in read mode
with open('input.in', 'r') as f:
    # read all lines and store them in a list
    lines = f.readlines()

# Define the cycling variable
cycling_var = list(np.arange(.92, 1.02, 0.01))

# Open the file `filenames.txt` inside the `dir_name` directory to store the names of the generated output files
with open(dir_name+'/'+'filenames.txt', 'w') as filenames:

    # Loop over the cycling variable and replace the value 1.1 in the lines
    for i, cv in enumerate(cycling_var):
        # replace the 4th line (index 3) with the cycling variable value
        lines[2] = str(cv) + '\n'

        # create the output file name with the cycling variable value
        output_file = dir_name+'/'+'input_{}.in'.format(cv)

        # write the updated lines to the output file
        with open(output_file, 'w') as f_out:
            f_out.writelines(lines)

        # Print a message indicating the output file name
        print('Cycle {}: output file saved as {}'.format(i, output_file))
        # Write the output file name to `filenames.txt`
        filenames.write(output_file+'\n')
