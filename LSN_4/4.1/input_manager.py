import os
import sys
import numpy as np

acceptedArgs = ['Liquid', 'Solid', 'Gas']
phase_name = sys.argv[1]
if phase_name in acceptedArgs:
    dir_name ='./'+ phase_name +"/inputs"
else: 
    raise KeyError('Arg not accepted, list of accepted args:', acceptedArgs)
if not os.path.exists(dir_name):

    print('Making a new directory')
    os.makedirs(dir_name)


# open the file in read mode
with open('input.in', 'r') as f:
    # read all lines and store them in a list
    lines = f.readlines()

# define the cycling variable
cycling_var = list(np.arange(.92, 1.02, 0.01))

with open(dir_name+'/'+'filenames.txt', 'w') as filenames:

    # loop over the cycling variable and replace the value 1.1 in the lines
    for i, cv in enumerate(cycling_var):
        # replace the 4th line (index 3) with the cycling variable value
        lines[2] = str(cv) + '\n'

        # create the output file name with the cycling variable value
        output_file = dir_name+'/'+'input_{}.in'.format(cv)

        # write the updated lines to the output file
        with open(output_file, 'w') as f_out:
            f_out.writelines(lines)

        # print a message
        print('Cycle {}: output file saved as {}'.format(i, output_file))
        filenames.write(output_file+'\n')
