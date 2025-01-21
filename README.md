# Use md.sh to run the ANM-highthroughput-pipeline.py script.
### the script itself creates all the necessery input files to run ANM-pathway and also creates all the necessery directory tree.
The bin-files contains all the binary files to run ANM-patway.

One thing that needs to be made sure that the prepare_input_strucutre_files.tcl file is copied into all the subdirectories under successful_eco
run this command to copy it into all the subdirectories - 

for dir in successful_eco/*/; do cp prepare_input_strucutre_files.tcl "$dir"; done

VMD is needed for step0
Biopython is needed in the environment for getting pdb information.

Helper scripts are also available. 
cleanup.sh deletes all the things created by ANM-pathway (except for the initial input files)

finding_files.sh helps to identify any sub dir where ANM failed
