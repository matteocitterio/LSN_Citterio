import os
import shutil
import sys
import numpy as np

'''
This python program takes care of changing the input parameters in order to simulate the system at different temperatures (and external field in the magnitude case)
'''

# List of accepted command line arguments
acceptedArgs = ['U', 'M', 'C','CHI']
acceptedArgs2 = ['Metro', 'Gibbs']

# Retrieve command line arguments for quantity and method
quantity = sys.argv[1]
method = sys.argv[2]

# Command line input management
if len(sys.argv) < 3:
    raise KeyError('Program ', sys.argv[0], ' usage: ', acceptedArgs, ' ', acceptedArgs2)

# Set the directory path based on quantity and method
if quantity in acceptedArgs:
    dir_name ='./'+ quantity +"/inputs/"+method
else: 
    raise KeyError('Arg not accepted, list of accepted args:', acceptedArgs)

# Check if the specified method is in the accepted list
if method in acceptedArgs2:
    print('Using ', method, ' algorithm')
else:
    raise KeyError('Arg not accepted, list of accepted args:', acceptedArgs2)

# If the directory doesn't exist, create a new one. Otherwise, delete old files.
if not os.path.exists(dir_name):
    print('Making a new directory with default params')
    os.makedirs(dir_name)
else:
    # remove files starting with 'input_'
    for file in os.listdir(dir_name):
        if file.startswith('input_'):
            os.remove(os.path.join(dir_name, file))

# Copy the default parameter file
src_file = './input.dat'
dst_file = os.path.join(dir_name, quantity + '.dat')
if not os.path.exists(dst_file):
    shutil.copy(src_file, dst_file)

# Open the file in read mode
with open(dir_name+'/'+quantity+'.dat', 'r') as f:
    # Read all lines and store them in a list
    lines = f.readlines()

# Define the cycling variable
temperatures = list(np.linspace(0.2, 3, 58))

with open(dir_name+'/'+'filenames.txt', 'w') as filenames:

    # Loop over the cycling variable and replace the value 1.1 in the lines
    for i, cv in enumerate(temperatures):
        # Replace the 2nd line (index 1) with the cycling variable value
        lines[1] = str(cv) + '\n'
        
        # Setting the algorithm
        if method == "Metro":
            lines[5] = str(1) +'\n'
        else: 
            lines[5] = str(0) +'\n'

        # If we are doing the magnitude, add the external field
        if quantity == 'M':
            lines[4] = str(0.02) +'\n'
        else:
            lines[4] = str(0.0) + '\n'

        # Create the output file name with the cycling variable value
        output_file = dir_name+'/'+'input_{}.dat'.format(cv)

        # Write the updated lines to the output file
        with open(output_file, 'w') as f_out:
            f_out.writelines(lines)

        # Print a message
        print('Cycle {}: output file saved as {}'.format(i, output_file))
        filenames.write(output_file+'\n')

# Remove the copied default parameter file
os.remove(dst_file)