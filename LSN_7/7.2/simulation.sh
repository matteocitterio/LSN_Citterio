#!/bin/sh#!/bin/bash

# Check if the input argument is provided
if [ $# -eq 0 ]; then
    echo "Error: Please provide an argument (Gas|Liquid|Solid)."
    exit 1
fi

# Validate the input argument
valid_materials=("Gas" "Liquid" "Solid")
material=$1
echo $material

if [[ ! " ${valid_materials[@]} " =~ " ${material} " ]]; then
    echo "Error: Invalid material. Please choose from Gas, Liquid, or Solid."
    exit 1
fi

# Set the folder name based on the material
if [ "$material" == "Gas" ]; then
    folder_name="gas"
elif [ "$material" == "Liquid" ]; then
    folder_name="liquid"
else
    folder_name="solid"
fi

# Remove files
rm "${material}/instant.epot"
rm "${material}/instant.pres"

# Generate input file name
input_file="input.${folder_name}"

# Copy input file
cp "${input_file}" "input.in"

# Run the program
./NVE_NVT.exe

