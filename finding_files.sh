#!/bin/bash

# Check if the main directory is provided as an argument
if [[ $# -eq 0 ]]; then
    echo "Usage: $0 <main_directory>"
    exit 1
fi

# Main directory is the first argument
main_dir="$1"

# Initialize counters and lists
count_with_both=0
missing=()

# Loop through all subdirectories in the main directory
for subdir in "$main_dir"/*/; do
    # Check if step5 folder and pathway.pdb exist inside the subdirectory
    if [[ -d "${subdir}step5" && -f "${subdir}step5/pathway.pdb" ]]; then
        ((count_with_both++))
    else
        missing+=("$subdir")
    fi
done

# Output the count of subdirectories meeting the criteria
echo "Number of subdirectories with 'step5' folder and 'pathway.pdb' file: $count_with_both"

# Output the subdirectories that do not meet the criteria
#if [[ ${#missing[@]} -gt 0 ]]; then
#    echo "Subdirectories missing either 'step5' folder or 'pathway.pdb' file:"
#    for miss in "${missing[@]}"; do
#        echo "$miss"
#    done
#else
#    echo "All subdirectories have 'step5' folder and 'pathway.pdb' file."
#fi

