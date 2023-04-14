import os
import shutil
import sys
import numpy as np

acceptedArgs = ['U', 'M', 'C','CHI']
acceptedArgs2 = ['Metro', 'Gibbs']
quantity = sys.argv[1]
method = sys.argv[2]

#Doing a bit of command line input managment
if len(sys.argv) < 3:

    raise KeyError('Program ', sys.argv[0], ' usage: ', acceptedArgs, ' ', acceptedArgs2)

if quantity in acceptedArgs:
    dir_name ='./'+ quantity +"/inputs/"+method
else: 
    raise KeyError('Arg not accepted, list of accepted args:', acceptedArgs)

if method in acceptedArgs2:
    print('Using ', method, ' algorithm')
else:
    raise KeyError('Arg not accepted, list of accepted args:', acceptedArgs2)

#If it doesnt exist create a new directory, else delete all the old files
if not os.path.exists(dir_name):
    print('Making a new directory with default params')
    os.makedirs(dir_name)
else:
    # remove files starting with 'input_'
    for file in os.listdir(dir_name):
        if file.startswith('input_'):
            os.remove(os.path.join(dir_name, file))

#copy the default parameter file
src_file = './input.dat'
dst_file = os.path.join(dir_name, quantity + '.dat')
if not os.path.exists(dst_file):
    shutil.copy(src_file, dst_file)

# open the file in read mode
with open(dir_name+'/'+quantity+'.dat', 'r') as f:
    # read all lines and store them in a list
    lines = f.readlines()

# define the cycling variable
temperatures = list(np.linspace(0.2, 3, 58))

with open(dir_name+'/'+'filenames.txt', 'w') as filenames:

    # loop over the cycling variable and replace the value 1.1 in the lines
    for i, cv in enumerate(temperatures):
        # replace the 4th line (index 3) with the cycling variable value
        lines[1] = str(cv) + '\n'
        
        #setting the algorithm
        if method == "Metro":
            lines[5] = str(1) +'\n'
        else: 
            lines[5] = str(0) +'\n'

        #if we are doing the magnitude, add the external field
        if quantity == 'M':
            lines[4] = str(0.02) +'\n'
        else:
            lines[4] = str(0.0) + '\n'

        # create the output file name with the cycling variable value
        output_file = dir_name+'/'+'input_{}.dat'.format(cv)

        # write the updated lines to the output file
        with open(output_file, 'w') as f_out:
            f_out.writelines(lines)

        # print a message
        print('Cycle {}: output file saved as {}'.format(i, output_file))
        filenames.write(output_file+'\n')

os.remove(dst_file)