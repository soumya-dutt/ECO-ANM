import os
import subprocess
import glob
import time
from Bio import PDB

def create_directories(base_dir):
    directories = ["step1", "step2", "step3/step3-1", "step3/step3-2", "step4", "step5"]
    for directory in directories:
        os.makedirs(os.path.join(base_dir, directory), exist_ok=True)

def parse_pdb_info(pdb_file):
    parser = PDB.PDBParser(QUIET=True)
    structure = parser.get_structure("protein", pdb_file)

    chain_residues = {}
    for chain in structure.get_chains():
        residues = [res.get_id()[1] for res in chain.get_residues() if PDB.is_aa(res)]
        if residues:
            chain_residues[chain.id] = f"{min(residues)} to {max(residues)}"

    ca_count = sum(1 for atom in structure.get_atoms() if atom.get_id() == "CA")
    return chain_residues, ca_count

def write_input_info_structures(base_dir, pdb1, pdb2, chains_residues):
    with open(os.path.join(base_dir, "INPUT_INFO_STRUCTURES"), "w") as f:
        f.write(f"{pdb1} {pdb2}\n")
        for chain, residue_range in chains_residues.items():
            f.write(f"{chain} {residue_range}\n")

def write_file_step1(base_dir, num_particles):
    content = f"""NUM_PARTICLES {num_particles}
CUT_OFFS 15.0 15.0
FORCE_CONSTANTS 0.1 0.1
CUSP_TOLERANCE 1.0e-5
NUM_IMAGES_INTERPOLATE 100"""
    with open(os.path.join(base_dir, "step1/INPUT_FIND_STRUCTURE_ON_CUSP"), "w") as f:
        f.write(content)

def write_file_step2(base_dir, num_particles):
    content = f"""NUM_PARTICLES {num_particles}
CUT_OFFS 15.0 15.0
FORCE_CONSTANTS 0.1 0.1
CUSP_TOLERANCE 1.0e-4
NUM_IMAGES_INTERPOLATE 100
STEP_SIZES 0.5 0.5
NUM_ITER_CHECK_STEP_SIZES 20
STEP_SIZE_REDUCTION_FACTORS 0.5 0.5
ENERGY_TOL_CONVERGENCE 0.001
NUM_ITERATIONS 10000"""
    with open(os.path.join(base_dir, "step2/INPUT_MINIMIZE_ON_CUSP"), "w") as f:
        f.write(content)

def write_file_step3(base_dir, num_particles, surface_index, folder):
    content = f"""NUM_PARTICLES {num_particles}
SURFACE_INDEX {surface_index}
CUT_OFF 15.0
FORCE_CONSTANT 0.1
STEP_SIZE 0.04
RMSD_PATHWAY 0.2
RMSD_DIFF_TOL 0.01
ENERGY_FROM_REFERENCE_TOL 0.0001
MAX_NUM_ITER 1000000"""
    with open(os.path.join(base_dir, folder, "INPUT_SLIDE_ONE_SURFACE"), "w") as f:
        f.write(content)

def write_file_step4(base_dir, num_particles):
    content = f"""NUM_PARTICLES {num_particles}
CUT_OFFS 15.0 15.0
FORCE_CONSTANTS 0.1 0.1"""
    with open(os.path.join(base_dir, "step4/INPUT_COLLECT_ENERGY"), "w") as f:
        f.write(content)

def run_command(command, cwd=None, timeout=None):
    try:
        subprocess.run(command, shell=True, cwd=cwd, check=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        print(f"Command '{command}' timed out after {timeout} seconds.")
    except subprocess.CalledProcessError as e:
        print(f"Command '{command}' failed with return code {e.returncode}.")

def count_files(pattern):
    return len(glob.glob(pattern))

def copy_files(source_pattern, destination):
    for file in glob.glob(source_pattern):
        subprocess.run(["cp", file, destination], check=True)

def check_and_copy(file_path, destination):
    if os.path.exists(file_path):
        subprocess.run(["cp", file_path, destination], check=True)
    else:
        print(f"Warning: File '{file_path}' does not exist. Skipping copy.")

def copy_all_initial_files(base_dir):
    initial_files = ["INPUT_STRUCTURE_1", "INPUT_STRUCTURE_2", "PDB_INFO", "REFERENCE_FOR_ALIGNMENT"]
    destinations = ["step1", "step2", "step3/step3-1", "step3/step3-2", "step4", "step5"]
    for dest in destinations:
        for file in initial_files:
            file_path = os.path.join(base_dir, file)
            check_and_copy(file_path, os.path.join(base_dir, dest))

def handle_step3(base_dir, num_particles, bin_files_dir):
    for i, folder in enumerate(["step3/step3-1", "step3/step3-2"], start=1):
        write_file_step3(base_dir, num_particles, i, folder)
        binary_file = os.path.join(bin_files_dir, "3_desc_one_surface_v2")
        check_and_copy(binary_file, os.path.join(base_dir, folder))
        copy_files(os.path.join(base_dir, "step2/minimized_struct_on_cusp"), os.path.join(base_dir, folder))
        run_command(f"./3_desc_one_surface_v2", cwd=os.path.join(base_dir, folder), timeout=300)

        # Count the number of output files and write num_structures_written files
        output_pattern = os.path.join(base_dir, folder, f"OUT_COORDS_SURFACE_{i}_*")
        num_structures = count_files(output_pattern)
        with open(os.path.join(base_dir, folder, f"num_structures_written_{i}"), "w") as f:
            f.write(str(num_structures))

def process_subdirectory(sub_dir, bin_files_dir):
    pdb_files = glob.glob(os.path.join(sub_dir, "*.pdb"))
    if len(pdb_files) != 2:
        print(f"Skipping {sub_dir} as it does not contain exactly 2 PDB files.")
        return

    pdb1, pdb2 = pdb_files
    chains_residues1, ca_count1 = parse_pdb_info(pdb1)
    chains_residues2, ca_count2 = parse_pdb_info(pdb2)

    chains_residues = {**chains_residues1, **chains_residues2}  # Merge chain info from both PDB files
    create_directories(sub_dir)
    write_input_info_structures(sub_dir, os.path.basename(pdb1), os.path.basename(pdb2), chains_residues)
    run_command(f"vmd -dispdev none -e prepare_input_strucutre_files.tcl", cwd=sub_dir)

    copy_all_initial_files(sub_dir)

    num_particles = ca_count1

    # Step 1
    write_file_step1(sub_dir, num_particles)
    binary_file = os.path.join(bin_files_dir, "1_locate_struct_on_cusp")
    check_and_copy(binary_file, os.path.join(sub_dir, "step1/"))
    run_command("./1_locate_struct_on_cusp", cwd=os.path.join(sub_dir, "step1"))

    # Step 2
    write_file_step2(sub_dir, num_particles)
    binary_file = os.path.join(bin_files_dir, "2_find_min_on_cusp_v3")
    check_and_copy(binary_file, os.path.join(sub_dir, "step2/"))
    copy_files(os.path.join(sub_dir, "step1/initial_struct_on_cusp"), os.path.join(sub_dir, "step2/"))
    run_command("./2_find_min_on_cusp_v3", cwd=os.path.join(sub_dir, "step2"))

    # Step 3
    handle_step3(sub_dir, num_particles, bin_files_dir)

    # Step 4
    write_file_step4(sub_dir, num_particles)
    binary_file = os.path.join(bin_files_dir, "4_collec_ener")
    check_and_copy(binary_file, os.path.join(sub_dir, "step4/"))
    copy_files(os.path.join(sub_dir, "step2/minimized_struct_on_cusp"), os.path.join(sub_dir, "step4/"))
    copy_files(os.path.join(sub_dir, "step3/step3-1/OUT_COORDS_SURFACE_1_*"), os.path.join(sub_dir, "step4/"))
    copy_files(os.path.join(sub_dir, "step3/step3-1/num_structures_written_1"), os.path.join(sub_dir, "step4/"))
    copy_files(os.path.join(sub_dir, "step3/step3-2/OUT_COORDS_SURFACE_2_*"), os.path.join(sub_dir, "step4/"))
    copy_files(os.path.join(sub_dir, "step3/step3-2/num_structures_written_2"), os.path.join(sub_dir, "step4/"))
    run_command("./4_collec_ener", cwd=os.path.join(sub_dir, "step4"))

    # Step 5
    binary_file = os.path.join(bin_files_dir, "5_makepathway")
    check_and_copy(binary_file, os.path.join(sub_dir, "step5/"))
    copy_files(os.path.join(sub_dir, "step4/COORDS_IMAGE_*"), os.path.join(sub_dir, "step5/"))
    num_images = count_files(os.path.join(sub_dir, "step5/COORDS_IMAGE_*"))
    run_command(f"./5_makepathway -n {num_particles} -i {num_images} -a 0", cwd=os.path.join(sub_dir, "step5"))

def main():
    main_directory = "successful_eco"  # Replace with your main directory path
    bin_files_dir = os.path.join(main_directory, "../bin-files")

    sub_dirs = [os.path.join(main_directory, d) for d in os.listdir(main_directory) if os.path.isdir(os.path.join(main_directory, d))]
    for sub_dir in sub_dirs:
        process_subdirectory(sub_dir, bin_files_dir)

if __name__ == "__main__":
    main()
