# src/cadock.py

# Import necessary libraries
import pubchempy as pcp
from rdkit import Chem
from rdkit.Chem import AllChem
import os
from Bio.PDB import PDBList
import pandas as pd
import numpy as np
from pymol import cmd
import subprocess
import sys

def generate_minimized_pdb(smiles, pdb_filename):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Invalid SMILES string: {smiles}")
        return False
    mol = Chem.AddHs(mol)
    try:
        AllChem.EmbedMolecule(mol, randomSeed=42)
    except:
        print(f"Failed to generate 3D coordinates for SMILES: {smiles}")
        return False
    try:
        AllChem.UFFOptimizeMolecule(mol, maxIters=200)
    except:
        print(f"Energy minimization failed for SMILES: {smiles}")
        return False
    try:
        Chem.SanitizeMol(mol)
    except:
        print(f"Sanitization failed for SMILES: {smiles}")
        return False
    Chem.MolToPDBFile(mol, pdb_filename)
    print(f"Minimized molecule saved as {pdb_filename}")
    return True

def download_pdb(pdb_id, download_dir):
    if not os.path.exists(download_dir):
        os.makedirs(download_dir)
    pdbl = PDBList()
    pdb_file_path = pdbl.retrieve_pdb_file(pdb_id, file_format='pdb', pdir=download_dir)
    return pdb_file_path

def remove_hetatm(input_pdb, output_pdb):
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            if not line.startswith('HETATM'):
                outfile.write(line)

def convert_pdb_to_pdbqt_receptor(input_pdb, output_pdbqt):
    subprocess.run(['obabel', '-i', 'pdb', input_pdb, '-o', 'pdbqt', '-O', output_pdbqt, '-xr', '-xn', '-xp'], check=True)

def convert_pdb_to_pdbqt_ligand(input_pdb, output_pdbqt):
    subprocess.run(['obabel', '-i', 'pdb', input_pdb, '-o', 'pdbqt', '-O', output_pdbqt, '-h'], check=True)

def run_command_with_output(command, log_file):
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    with open(log_file, 'w') as log:
        for line in process.stdout:
            sys.stdout.write(line)
            log.write(line)
            sys.stdout.flush()
    return process.wait()

def perform_docking(smiles_list, PDB_ID):
    folder_name = 'docking_results'
    receptor_name = PDB_ID

    # Create the folder
    if not os.path.exists(folder_name):
        os.mkdir(folder_name)

    print(f"Receptor Name: {receptor_name}")
    print(f"Number of ligands: {len(smiles_list)}")

    # Generate and Pre-process Ligands
    valid_smiles = []
    for i, smiles in enumerate(smiles_list):
        pdb_filename = f'{folder_name}/ligand_{i+1}.pdb'
        if generate_minimized_pdb(smiles, pdb_filename):
            valid_smiles.append(smiles)
        else:
            print(f"Skipping invalid SMILES: {smiles}")

    print(f"Number of valid SMILES processed: {len(valid_smiles)}")

    # Download and Pre-process Receptor
    downloaded_pdb_path = download_pdb(PDB_ID, folder_name)
    os.rename(downloaded_pdb_path, f'{folder_name}/{receptor_name}_dirty.pdb')

    # Remove HETATM from PDB file
    remove_hetatm(f'{folder_name}/{receptor_name}_dirty.pdb', f'{folder_name}/{receptor_name}.pdb')

    # Define Box (p2rank)
    subprocess.run(['./p2rank_2.4.2/prank', 'predict', '-f', f'{folder_name}/{receptor_name}.pdb'], check=True)

    df = pd.read_csv(f'p2rank_2.4.2/test_output/predict_{receptor_name}/{receptor_name}.pdb_predictions.csv')
    center_x, center_y, center_z = float(df['   center_x'].iloc[0]), float(df['   center_y'].iloc[0]), float(df['   center_z'].iloc[0])

    pred = pd.read_csv(f'p2rank_2.4.2/test_output/predict_{receptor_name}/{receptor_name}.pdb_residues.csv')
    pocket1 = pred[pred[' pocket'] == 1]
    resi = '+'.join([str(i) for i in pocket1[' residue_label']])

    cmd.load(f'{folder_name}/{receptor_name}.pdb')
    cmd.select('pocket1', f'resi {resi}')
    cmd.show('cartoon')
    alpha_carbons = []
    cmd.iterate_state(1, 'pocket1 and name CA', 'alpha_carbons.append([x, y, z])', space={'alpha_carbons': alpha_carbons})

    alpha_carbons = np.array(alpha_carbons)
    min_coords = np.min(alpha_carbons, axis=0)
    max_coords = np.max(alpha_carbons, axis=0)
    cube_size = max_coords - min_coords
    Size_x, Size_y, Size_z = cube_size

    print(center_x, center_y, center_z, Size_x, Size_y, Size_z)

    # Docking
    receptor_pdb = f"{folder_name}/{receptor_name}.pdb"
    receptor_pdbqt = f"{folder_name}/{receptor_name}.pdbqt"

    print("Converting receptor to PDBQT format...")
    convert_pdb_to_pdbqt_receptor(receptor_pdb, receptor_pdbqt)
    print("Receptor conversion complete.")

    results_file = f"{folder_name}/docking_results.txt"
    with open(results_file, 'w') as f:
        f.write("SMILES,Docking Score\n")

    for i, smiles in enumerate(smiles_list):
        print(f"\nProcessing ligand {i+1} of {len(smiles_list)}")
        print(f"SMILES: {smiles}")

        ligand_pdb = f"{folder_name}/ligand_{i+1}.pdb"
        ligand_pdbqt = f"{folder_name}/ligand_{i+1}.pdbqt"

        print("Converting ligand to PDBQT format...")
        convert_pdb_to_pdbqt_ligand(ligand_pdb, ligand_pdbqt)
        print("Ligand conversion complete.")

        output = f"{receptor_name}_ligand_{i+1}.pdbqt"
        log_file = f"{folder_name}/vina_log_{i+1}.txt"
        vina_command = [
            'vina',
            '--receptor', receptor_pdbqt,
            '--ligand', ligand_pdbqt,
            '--out', f'{folder_name}/{output}',
            '--center_x', str(center_x),
            '--center_y', str(center_y),
            '--center_z', str(center_z),
            '--size_x', str(Size_x),
            '--size_y', str(Size_y),
            '--size_z', str(Size_z)
        ]

        print("Starting Vina docking...")
        exit_code = run_command_with_output(vina_command, log_file)

        if exit_code == 0:
            print("Vina docking completed successfully.")
            with open(log_file, 'r') as log:
                score = "N/A"
                for line in log:
                    if line.startswith('   1'):
                        score = line.split()[1]
                        break
            print(f"Best docking score: {score}")
        else:
            print(f"Error running Vina for ligand {i+1}. Check the log file for details.")
            score = "Error"

        with open(results_file, 'a') as f:
            f.write(f"{smiles},{score}\n")

    print(f"\nDocking complete. Results have been saved to {results_file}")
    print(f"Individual log files for each ligand are saved in the {folder_name} directory.")

    # Visualization and Alignment are skipped as they require interactive display capabilities
    # which are not suitable for Colab's backend execution without additional setup
