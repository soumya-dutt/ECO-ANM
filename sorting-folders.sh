#!/bin/bash

# Define the main directory and output directory
MAIN_DIR="pdbs_clean31"
OUTPUT_DIR="successful-eco"

# Create the output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR"

# Iterate through all subdirectories in the main directory
for SUBDIR in "$MAIN_DIR"/*/; do
  # Define the path to the eco.log file
  LOG_FILE="$SUBDIR/eco.log"
  
  # Check if the eco.log file exists
  if [ -f "$LOG_FILE" ]; then
    # Read the first word of the file
    FIRST_WORD=$(head -n 1 "$LOG_FILE" | awk '{print $1}')
    
    # Check if the first word is SUCCESS
    if [ "$FIRST_WORD" == "SUCCESS" ]; then
      # Copy the subdirectory to the output directory
      cp -r "$SUBDIR" "$OUTPUT_DIR"
    fi
  fi
done

echo "Subdirectories with SUCCESS in eco.log have been copied to $OUTPUT_DIR."
