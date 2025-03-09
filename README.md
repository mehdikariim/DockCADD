# DockCADD v2

**DockCADD v2** is a streamlined and automated computational framework designed to facilitate molecular docking and drug discovery. It requires minimal input from users and utilizes advanced tools to provide accurate docking results. 

This new version of **DockCADD** is a lightweight, integrated workflow for structure-based drug design. It automatically prepares a receptor (extracting only a specified chain, with a fallback if the chain isn’t found), repairs missing residues/atoms via [PDBFixer](https://github.com/openmm/pdbfixer), predicts the binding pocket with [p2rank](https://github.com/rdk/p2rank), and docks ligands using [AutoDock Vina](https://github.com/ccsb-scripps/AutoDock-Vina). Ligands can be provided as a list of SMILES strings or as an SDF file, and multiple conformers are generated for each molecule using [RDKit](https://www.rdkit.org/). 

## Features

- **Receptor Preparation**:  
  - Downloads a protein structure from the PDB.
  - Extracts only chain A (or falls back to the full structure if chain A isn’t found).
  - Repairs missing residues, atoms, and adds hydrogens using PDBFixer.

- **Pocket Prediction**:  
  - Uses p2rank to predict the binding pocket and extract its center.

- **Ligand Preparation**:  
  - Accepts ligands as a list of SMILES strings and/or an SDF file.
  - Generates multiple 3D conformers per ligand with RDKit.
  - Writes each conformer to a separate PDB file.

- **Docking**:  
  - Converts the receptor and ligand files to PDBQT format using OpenBabel.
  - Docks each ligand conformer with AutoDock Vina.
  - Parses the best docking pose and merges it with the receptor to generate a final complex.

- **Visualization**:  
  - Includes an optional PyMOL visualization function to generate a static PNG snapshot of a final complex.

## Installation

1. **Clone the Repository:**

   ```bash
   git clone https://github.com/mehdikariim/DockCADD-v2.git
   cd DockCADD-v2

2. **Run the Setup Script:**

This script installs all system packages (e.g., PyMOL, OpenBabel, Java), AutoDock Vina, p2rank, and the required Python libraries (including PDBFixer, OpenMM, and RDKit).

     ```bash
     bash scripts/setup.sh

2. **Run the Setup Script:**

3. **Usage:**
You can use the provided Python package to perform docking. Below are two example usage scenarios:

Example 1: **Docking Using a List of SMILES**
    ```bash
    from src.dockcadd import perform_docking

    # Define your list of ligand SMILES and target receptor PDB ID
    smiles_list = ["CCOc1ccc(CC(=O)NC)cc1", "CCCC(=O)NCC1=CC=CC=C1"]
    pdb_id = "5ZMA"

    # Run docking (generates 3 conformers per ligand by default)
    perform_docking(smiles_list=smiles_list, sdf_file=None, pdb_id=pdb_id, num_confs=3, docking_folder="docking_results")

Example 2: **Docking Using an SDF File**
    ```bash
    from src.dockcadd import perform_docking
    
    # Provide the path to your SDF file containing ligands
    sdf_file = "path/to/your_ligands.sdf"
    pdb_id = "5ZMA"
    
    # Run docking using the SDF file (3 conformers per ligand)
    perform_docking(smiles_list=None, sdf_file=sdf_file, pdb_id=pdb_id, num_confs=3, docking_folder="docking_results")


# License
This project is licensed under the MIT License.

# Acknowledgments
AutoDock Vina: AutoDock Vina 1.2.5
p2rank: p2rank 2.4.2
PDBFixer & OpenMM: PDBFixer
RDKit: RDKit


