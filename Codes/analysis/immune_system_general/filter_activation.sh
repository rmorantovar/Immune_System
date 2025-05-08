#!/bin/bash

# Set the input and output files
path_folder = "~/"
input_file= path_folder + "input.bin"
output_file="output.bin"

# Set the column number to use as the condition
column_num=1

# Set the value to use as the condition
value="foo"

# Iterate over each row in the input file
while read row; do
  # Split the row into columns
  columns=($row)

  # Check the value of the specified column
  if [ "${columns[$column_num-1]}" != "$value" ]; then
    # If the value does not match the condition, write the row to the output file
    echo "$row" >> "$output_file"
  fi
done < "$input_file"
