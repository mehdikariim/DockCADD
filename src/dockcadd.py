
import os
import re
import subprocess
import pandas as pd

from Bio.PDB import PDBList
from pdbfixer import PDBFixer
from openmm.app import PDBFile

from rdkit import Chem
from rdkit.Chem import AllChem

########################################
# Helper Functions
########################################

def run_command_with_live_output(command, log_file):
    """
    Runs a command in a subprocess, writing output to both console and a log file.
    """
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
    with open(log_file, 'w') as lf:
        for line in process.stdout:
            print(line, end="")  # Print to console
            lf.write(line)
    return process.wait()


def keep_only_chain_A_with_fallback(input_pdb, output_pdb):
    """
    Extracts lines for chain A from input_pdb.
    Writes only ATOM/HETATM lines (and TER if chain A) to output_pdb.
    If no chain A lines are found, copies the full file.
    """
    chain_a_count = 0
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            chain_id = line[21] if len(line) >= 22 else None
            if line.startswith(("ATOM", "HETATM")):
                if chain_id == 'A':
                    outfile.write(line)
                    chain_a_count += 1
            elif line.startswith("TER"):
                if chain_id == 'A':
                    outfile.write(line)
            elif line.startswith("END"):
                outfile.write(line)
    if chain_a_count == 0:
        print("[WARN] No chain A lines found; using full PDB instead.")
        with open(input_pdb, 'r') as inf, open(output_pdb, 'w') as outf:
            outf.write(inf.read())


def remove_hetatm(input_pdb, output_pdb):
    """
    Removes all HETATM lines from input_pdb.
    This effectively removes any co-ligand present in the receptor.
    """
    with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
        for line in infile:
            if not line.startswith("HETATM"):
                outfile.write(line)


def fix_with_pdbfixer(pdb_in, pdb_out):
    """
    Uses PDBFixer to add missing residues, atoms, and hydrogens (pH 7).
    """
    fixer = PDBFixer(filename=pdb_in)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(7.0)
    with open(pdb_out, 'w') as f:
        PDBFile.writeFile(fixer.topology, fixer.positions, f)


def convert_pdb_to_pdbqt_receptor(input_pdb, output_pdbqt):
    """
    Converts receptor PDB to PDBQT using OpenBabel.
    """
    cmd = [
        "obabel", "-i", "pdb", input_pdb,
        "-o", "pdbqt", "-O", output_pdbqt,
        "-xr", "-xn", "-xp"
    ]
    subprocess.run(cmd, check=True)


def convert_pdb_to_pdbqt_ligand(input_pdb, output_pdbqt):
    """
    Converts ligand PDB to PDBQT using OpenBabel.
    """
    cmd = [
        "obabel", "-i", "pdb", input_pdb,
        "-o", "pdbqt", "-O", output_pdbqt,
        "-h"
    ]
    subprocess.run(cmd, check=True)


########################################
# p2rank Pocket Prediction
########################################

def run_p2rank_and_get_center(receptor_pdb, pdb_id):
    """
    Runs p2rank on the receptor PDB and returns the top pocket center (x,y,z).
    """
    p2rank_exec = os.path.join(os.getcwd(), "p2rank_2.4.2", "prank")
    if not os.path.isfile(p2rank_exec):
        raise FileNotFoundError(f"p2rank not found at {p2rank_exec}")
    if not os.access(p2rank_exec, os.X_OK):
        os.chmod(p2rank_exec, 0o755)
    cmd = [p2rank_exec, "predict", "-f", receptor_pdb]
    log_file = f"p2rank_{pdb_id}.log"
    ret = run_command_with_live_output(cmd, log_file)
    if ret != 0:
        raise RuntimeError("p2rank prediction failed.")
    base_name = os.path.splitext(os.path.basename(receptor_pdb))[0]
    predictions_csv = f"p2rank_2.4.2/test_output/predict_{base_name}/{base_name}.pdb_predictions.csv"
    df = pd.read_csv(predictions_csv, skipinitialspace=True)
    df.columns = [c.strip().lower() for c in df.columns]
    cx = float(df["center_x"].iloc[0])
    cy = float(df["center_y"].iloc[0])
    cz = float(df["center_z"].iloc[0])
    print(f"[p2rank] Pocket center: ({cx}, {cy}, {cz})")
    return (cx, cy, cz)


########################################
# Ligand Preparation (SMILES and/or SDF)
########################################

def generate_multiple_conformers(mol, num_confs=3):
    """
    Generates multiple 3D conformers for a molecule using RDKit.
    """
    mol = Chem.AddHs(mol)
    cids = AllChem.EmbedMultipleConfs(mol, numConfs=num_confs, randomSeed=42)
    for cid in cids:
        if AllChem.MMFFHasAllMoleculeParams(mol):
            AllChem.MMFFOptimizeMolecule(mol, confId=cid, maxIters=200)
        else:
            AllChem.UFFOptimizeMolecule(mol, confId=cid, maxIters=200)
    return mol


def prepare_ligands(smiles_list=None, sdf_file=None, num_confs=3, out_dir="ligand_prep"):
    """
    Prepares ligand files from a list of SMILES and/or an SDF file.
    For each valid molecule, generates multiple conformers and writes each as a separate PDB.
    Returns a list of (pdb_filepath, label).
    """
    if not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    results = []

    def write_confs_to_pdb(mol, base_name):
        conf_ids = [conf.GetId() for conf in mol.GetConformers()]
        out_paths = []
        for i, cid in enumerate(conf_ids, start=1):
            tmp_mol = Chem.Mol(mol, False, cid)
            pdb_name = f"{base_name}_conf{i}.pdb"
            pdb_path = os.path.join(out_dir, pdb_name)
            Chem.MolToPDBFile(tmp_mol, pdb_path)
            out_paths.append((pdb_path, pdb_name.replace(".pdb", "")))
        return out_paths

    # Process SMILES if provided
    if smiles_list:
        for idx, smi in enumerate(smiles_list, start=1):
            mol = Chem.MolFromSmiles(smi)
            if not mol:
                print(f"[LigandPrep] Warning: invalid SMILES skipped: {smi}")
                continue
            mol3d = generate_multiple_conformers(mol, num_confs=num_confs)
            base_name = f"lig_{idx}"
            results.extend(write_confs_to_pdb(mol3d, base_name))
    
    # Process SDF if provided
    if sdf_file and os.path.isfile(sdf_file):
        suppl = Chem.SDMolSupplier(sdf_file, removeHs=False)
        mol_count = 0
        for i, mol in enumerate(suppl):
            if mol is None:
                print(f"[LigandPrep] Warning: skipping invalid SDF record {i}.")
                continue
            mol_count += 1
            name = mol.GetProp("_Name") if mol.HasProp("_Name") else f"sdf_{mol_count}"
            mol3d = generate_multiple_conformers(mol, num_confs=num_confs)
            base_name = f"{name}_{mol_count}"
            results.extend(write_confs_to_pdb(mol3d, base_name))
    
    if not results:
        print("[LigandPrep] No valid ligands found.")
    else:
        print(f"[LigandPrep] Prepared {len(results)} ligand conformers.")
    return results


########################################
# Main Docking Workflow
########################################

def perform_docking(smiles_list=None, sdf_file=None, pdb_id="5ZMA", num_confs=3, docking_folder="docking_results"):
    """
    Main docking workflow:
      1) Prepare receptor: download PDB, extract chain A (fallback if needed), remove HETATM (co-ligand), and fix with PDBFixer.
      2) Run p2rank to get pocket center (for a 20x20x20 box).
      3) Convert receptor to PDBQT.
      4) Prepare ligands from SMILES and/or SDF (generate multiple conformers).
      5) For each ligand conformer, convert to PDBQT, dock with AutoDock Vina,
         convert best pose to PDB, and merge with receptor to form the final complex.
      6) Write docking scores to a CSV file.
    """
    if not (smiles_list or (sdf_file and os.path.isfile(sdf_file))):
        print("[ERROR] No valid ligand input provided (neither SMILES nor SDF).")
        return

    if not os.path.exists(docking_folder):
        os.makedirs(docking_folder, exist_ok=True)

    print(f"[Docking] Preparing receptor {pdb_id} ...")
    pdbl = PDBList()
    raw_file = pdbl.retrieve_pdb_file(pdb_id, file_format="pdb", pdir=docking_folder)
    raw_pdb = os.path.join(docking_folder, f"{pdb_id}_raw.pdb")
    os.rename(raw_file, raw_pdb)
    chainA_file = os.path.join(docking_folder, f"{pdb_id}_chainA_tmp.pdb")
    keep_only_chain_A_with_fallback(raw_pdb, chainA_file)
    
    # Remove HETATM lines to eliminate co-ligands from the receptor
    receptor_no_het = os.path.join(docking_folder, f"{pdb_id}_chainA_nohet.pdb")
    remove_hetatm(chainA_file, receptor_no_het)
    
    receptor_pdb = os.path.join(docking_folder, f"{pdb_id}_prepared.pdb")
    fix_with_pdbfixer(receptor_no_het, receptor_pdb)

    # p2rank pocket prediction
    cx, cy, cz = run_p2rank_and_get_center(receptor_pdb, pdb_id)
    box_size = 20.0
    print(f"[Docking] Using docking box (20x20x20) centered at ({cx}, {cy}, {cz})")

    # Convert receptor to PDBQT
    receptor_pdbqt = os.path.join(docking_folder, f"{pdb_id}_prepared.pdbqt")
    convert_pdb_to_pdbqt_receptor(receptor_pdb, receptor_pdbqt)

    # Ligand preparation (from SMILES and/or SDF)
    lig_out_dir = os.path.join(docking_folder, "ligands")
    ligand_list = prepare_ligands(smiles_list=smiles_list, sdf_file=sdf_file, num_confs=num_confs, out_dir=lig_out_dir)
    if not ligand_list:
        print("[Docking] No valid ligands to dock. Exiting.")
        return

    results_csv = os.path.join(docking_folder, "docking_results.csv")
    with open(results_csv, "w") as rf:
        rf.write("LigandLabel,Score\n")

    for i, (lig_pdb, label) in enumerate(ligand_list, start=1):
        print(f"\n[Docking] Processing ligand conformer {label} ...")
        ligand_pdbqt = os.path.join(docking_folder, f"{label}.pdbqt")
        convert_pdb_to_pdbqt_ligand(lig_pdb, ligand_pdbqt)

        out_pdbqt = os.path.join(docking_folder, f"{label}_out.pdbqt")
        log_file = os.path.join(docking_folder, f"{label}_vina.log")

        vina_cmd = [
            "vina",
            "--receptor", receptor_pdbqt,
            "--ligand", ligand_pdbqt,
            "--out", out_pdbqt,
            "--center_x", str(cx),
            "--center_y", str(cy),
            "--center_z", str(cz),
            "--size_x", str(box_size),
            "--size_y", str(box_size),
            "--size_z", str(box_size),
            "--num_modes", "10"
        ]
        ret_code = run_command_with_live_output(vina_cmd, log_file)
        best_score = "N/A"
        if ret_code == 0:
            with open(log_file, "r") as lg:
                for line in lg:
                    if re.match(r"^\s*1\s+", line):
                        parts = line.split()
                        if len(parts) >= 2:
                            best_score = parts[1]
                        break
        else:
            best_score = "ERROR"
        print(f"[Docking] Best score for {label}: {best_score}")
        with open(results_csv, "a") as rf:
            rf.write(f"{label},{best_score}\n")

        # Convert best pose from PDBQT to PDB
        docked_pdb = os.path.join(docking_folder, f"{label}_docked.pdb")
        subprocess.run([
            "obabel", "-ipdbqt", out_pdbqt,
            "-opdb", "-O", docked_pdb, "-d"
        ], check=True)

        # Merge receptor and docked ligand to form final complex
        final_complex = os.path.join(docking_folder, f"{label}_complex.pdb")
        with open(final_complex, "w") as fc:
            with open(receptor_pdb, "r") as recf:
                for line in recf:
                    if line.startswith("END"):
                        continue
                    fc.write(line)
            fc.write("TER\n")
            with open(docked_pdb, "r") as ligf:
                for line in ligf:
                    fc.write(line)
            fc.write("END\n")
        print(f"[Docking] Final complex saved as: {final_complex}")

    print(f"\n[DONE] Docking complete. Results saved in {results_csv}")
