#!/bin/bash

# Define folder names to delete
FOLDERS_TO_DELETE=("step1" "step2" "step3" "step4" "step5")

# Define exact file names to delete
FILES_TO_DELETE=("1_CA.pdb" "2_CB.pdb" "INPUT_STRUCTURE_1" "INPUT_STRUCTURE_2" \
"PDB_INFO" "REFERENCE_FOR_ALIGNMENT" "INPUT_INFO_STRUCTURES")

# Define file patterns to delete (e.g., *.log)
FILE_PATTERNS=("*.txt")

# Main directory to process (use current directory if not specified)
MAIN_DIR=successful-eco

echo "Starting cleanup in: $MAIN_DIR"

# Loop through folders and delete them
for folder in "${FOLDERS_TO_DELETE[@]}"; do
    echo "Searching for folders: $folder"
    find "$MAIN_DIR" -type d -name "$folder" -exec rm -rf {} + -print
done

# Loop through files and delete them
for file in "${FILES_TO_DELETE[@]}"; do
    echo "Searching for files: $file"
    find "$MAIN_DIR" -type f -name "$file" -exec rm -f {} + -print
done

# Loop through file patterns and delete them
for pattern in "${FILE_PATTERNS[@]}"; do
    echo "Searching for files matching pattern: $pattern"
    find "$MAIN_DIR" -type f -name "$pattern" -exec rm -f {} + -print
done

echo "Cleanup complete!"

