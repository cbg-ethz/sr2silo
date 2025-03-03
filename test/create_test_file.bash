#!/bin/bash

# Define output file
output_file="test_file2.ndjson"

# Remove file if it already exists
test -f "$output_file" && rm "$output_file"

# Generate JSON entries from 2000 down to 1
for ((i=2000; i>=1; i--)); do
    echo "{\"N_length\":$i}" >> "$output_file"
done

# Print completion message
echo "File '$output_file' created successfully."
