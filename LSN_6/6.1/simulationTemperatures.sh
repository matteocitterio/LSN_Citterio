#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 2 ]; then
    echo "Usage: $0 [U|M|C|CHI] [Metro|Gibbs]"
    exit 1
fi

# Set the variable from the command line argument
case $1 in
    U) input_base_dir="U/inputs"; results_base_dir="U/results"; output_file="output.ene.0"; output_prefix="output.ene_";;
    C) input_base_dir="C/inputs"; results_base_dir="C/results"; output_file="output.heat.0"; output_prefix="output.heat_";;
    M) input_base_dir="M/inputs"; results_base_dir="M/results"; output_file="output.mag.0"; output_prefix="output.mag_";;
    CHI) input_base_dir="CHI/inputs"; results_base_dir="CHI/results"; output_file="output.chi.0"; output_prefix="output.chi_";;
    *) echo "Invalid argument. Usage: $0 [U|M|C|CHI] [Metro|Gibbs]"; exit 1;;
esac

# Set the variable from the second command line argument
case $2 in
    Metro) input_dir="$input_base_dir/Metro"; results_dir="$results_base_dir/Metro";;
    Gibbs) input_dir="$input_base_dir/Gibbs"; results_dir="$results_base_dir/Gibbs";;
    *) echo "Invalid method. Usage: $0 [U|M|C|CHI] [Metro|Gibbs]"; exit 1;;
esac


# Count the number of cycles
num_cycles=$(ls "$input_dir" | wc -l)

echo "Number of simulations: $num_cycles"

# Create the results base directory if it doesn't exist
if [ ! -d "$results_base_dir" ]; then
    mkdir "$results_base_dir"
fi

# Clear the results directory if it exists, otherwise create it
if [ -d "$results_dir" ]; then
    rm -rf "$results_dir"/*
else
    mkdir "$results_dir"
fi

# Loop over the cycles
for ((i=2; i<=$num_cycles; i++))
do
    # Copy one of the files from the input directory to 'input.dat'
    file=$(ls "$input_dir" | sed -n "${i}p")
    cp "$input_dir/$file" input.dat
    
    # Run the program with the current cycle variable
    ./Monte_Carlo_ISING_1D.exe >/dev/null
    

    # Extract the number from the input file name and format it with two decimal places
    num=$(echo "$file" | sed 's/input_//;s/\.dat//')
    num=$(printf "%.2f" $num | sed 's/0*$//;s/\.$/.00/')

    echo "Doing simulation at temp $num"

    # Copy the output file to a file with the formatted number in the results directory
    cp "$output_file" "$results_dir/$output_prefix$num.dat"
    
done

