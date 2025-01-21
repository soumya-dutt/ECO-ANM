# Prepares two PDBs with CAs from user input.
# This script is intended for VMD.
#
# User input is stored in a file called 'INPUT_INFO_STRUCTURES'
# Format:
# -----------------------------------------------------------------
# (Input PDB file for end point 1) (Input PDB file for end point 2)
# (Chain id) (Residue numbers to include for this chain, VMD syntax) [One chain per line]
# -----------------------------------------------------------------
#
# Tasks performed by this script
# 1. Read the 'INPUT_INFO_STRUCTURES' file
# 2. Extract C-alphas for the given selections and align them with C-alphas of structure 1 as reference
# 3. Write two PDB files with aligned coordinates called '1_CA.pdb' and '2_CA.pdb'.
# 4. Write two myxyz files called 'INPUT_STRUCTURE_1' and 'INPUT_STRUCTURE_2'. 
#    Also write the 'REFERENCE_FOR_ALIGNMENT' file which is same as 'INPUT_STRUCTURE_1'.
# 5. Write the 'PDB_INFO' file from '1_CA.pdb'
# 6. Write a log file with some output.
###################################################################################################################


########## Tcl procedures ###########

# Procedure for writing a myxyz file from a VMD atom selection
proc Write_myxyz_from_vmd_atom_selection {sel filename_myxyz} {
    set coords [$sel get {x y z}]
    set out [open "$filename_myxyz" w]
    foreach atom_coord $coords {
        puts $out $atom_coord
    }
    close $out
}


# Write the PDB_INFO file
# Format:
# ---------------------------------
#  Number of atoms (%d)
#  Each line will have the follwoing fields separated by white space
#  (Atom_serial_number) (Atom_name) (Three_letter_residue_name) (Chain_identifier) (Residue_serial_number) (Element_symbol)
# ---------------------------------
proc Write_PDB_INFO {file_name_PDB} {
    mol load pdb $file_name_PDB
    set sel [atomselect top "all"]

    set NUM_ATOMS [$sel num]

    set SERIAL_NUMBERS [$sel get serial]
    set ATOM_NAMES [$sel get name]
    set RESIDUE_NAMES [$sel get resname]
    set CHAIN_IDS [$sel get chain]
    set RESIDUE_NUMBERS [$sel get resid]
    set ELEMENT_SYMBOLS [$sel get element]

    set out [open "PDB_INFO" w]
    puts $out $NUM_ATOMS
    for {set i 0} {$i < $NUM_ATOMS} {incr i} {
        puts $out "[lindex $SERIAL_NUMBERS $i] [lindex $ATOM_NAMES $i] [lindex $RESIDUE_NAMES $i] [lindex $CHAIN_IDS $i] [lindex $RESIDUE_NUMBERS $i] [lindex $ELEMENT_SYMBOLS $i]"
    }
    close $out
}


########## End of Tcl procedures ##########

###################################################################################################################

# Read 'INPUT_INFO_STRUCTURES'. Extract PDB file names, chain ids and atom selections.
set in [open "INPUT_INFO_STRUCTURES" r]
set file_data [read $in]
close $in

set data [split $file_data "\n"]
set num_elems [llength $data]
set num_lines [expr {$num_elems - 1}]

for {set i 0} {$i < $num_lines} {incr i} {
    set temp_list ""
    set item [lindex $data $i]
    if {$i == 0} {
        set each_line [split $item " "]
	set input_pdb_filename_1 [lindex $each_line 0]
        set input_pdb_filename_2 [lindex $each_line 1]
    } else {
	set each_line [split $item " "]
	set num_elements_one_line [llength $each_line]
	lappend chain_ids [lindex $each_line 0]
	for {set j 1} {$j < $num_elements_one_line} {incr j} {
	    lappend temp_list [lindex $each_line $j]
	}
	lappend atom_selections $temp_list
    }   
}

set num_chains [llength $atom_selections]
#set final_selection "protein and "
set final_selection ""
for {set i 0} {$i < $num_chains} {incr i} {
    set final_selection [concat $final_selection "(chain " [lindex $chain_ids $i] " and name CA and resid " [lindex $atom_selections $i] ") "]
    if {$i != [expr {$num_chains - 1}]} {
        set final_selection [concat $final_selection " or "]
    } 
}


# Read two input PDB files.
mol load pdb $input_pdb_filename_1
set sel_1 [atomselect top "$final_selection"]
mol load pdb $input_pdb_filename_2
set sel_2 [atomselect top "$final_selection"]


# Number of atoms in two selections
set num_atoms_1 [$sel_1 num]
set num_atoms_2 [$sel_2 num]
if {$num_atoms_1 != $num_atoms_2} {
    puts "ERROR: Different number of atoms extracted from two structures based on input atom selections."
    exit
} else {
    set num_atoms $num_atoms_1
}


# Calculate the initial RMSD without alignment.
set initial_rmsd [measure rmsd $sel_1 $sel_2]


# Align
set transformation_matrix [measure fit $sel_2 $sel_1]
$sel_2 move $transformation_matrix


# Calculate the final RMSD after alignment
set rmsd_after_alignment [measure rmsd $sel_1 $sel_2]


# Write two PDB files with C-alphas
$sel_1 writepdb 1_CA.pdb
$sel_2 writepdb 2_CB.pdb


# Write myxyz files: 'INPUT_STRUCTURE_1', 'INPUT_STRUCTURE_1' and 'REFERENCE_FOR_ALIGNMENT' (which is same as 'INPUT_STRUCTURE_1')
set filename_myxyz INPUT_STRUCTURE_1
Write_myxyz_from_vmd_atom_selection $sel_1 $filename_myxyz
set filename_myxyz INPUT_STRUCTURE_2
Write_myxyz_from_vmd_atom_selection $sel_2 $filename_myxyz
set filename_myxyz REFERENCE_FOR_ALIGNMENT
Write_myxyz_from_vmd_atom_selection $sel_1 $filename_myxyz


# Write the 'PDB_INFO' file
set file_name_PDB 1_CA.pdb
Write_PDB_INFO $file_name_PDB


# Write a log file called 'structure_prep.log'
set out [open "structure_prep.log" w]
puts $out "Input filename for end point 1: $input_pdb_filename_1"
puts $out "Input filename for end point 2: $input_pdb_filename_2"
puts $out "String used for VMD atom selection: $final_selection"
puts $out "Number of atoms: $num_atoms"
puts $out "Intial RMSD: $initial_rmsd"
puts $out "RMSD after alignment: $rmsd_after_alignment"
close $out


exit











 
