# DockCADD

**DockCADD** is a streamlined and automated computational framework designed to facilitate molecular docking and drug discovery. It requires minimal input from users and utilizes advanced tools to provide accurate docking results. 

## Features

- **Automation**: Reduces user intervention, making the process fast and reproducible.
- **Efficiency**: Processes large compound libraries effectively.
- **Cost-Effective**: Uses open-source tools, making it accessible for academic and industrial research.
- **Integration**: Combines RDKit, P2Rank, and AutoDock Vina for robust molecular docking.

## Core Workflow

1. **Inputs**:
   - **SMILES**: Simplified Molecular Input Line Entry System representing the chemical structure of ligands.
   - **PDB ID**: Protein Data Bank Identifier specifying the target protein structure.

2. **Workflow**:
   - **Ligand Preparation**: Generates 3D conformations, performs energy minimization, and converts ligands to PDBQT format.
   - **Protein Preparation**: Retrieves and cleans the protein structure, identifies binding pockets, and defines a docking grid box.
   - **Docking**: Utilizes AutoDock Vina for virtual screening, outputting docking scores and top binding poses.

3. **Output**:
   - **Docking Scores**: Quantify binding affinity.
   - **Ranked Ligands**: Based on docking performance.
   - **Additional Data**: Binding poses and interaction details.

## Installation

### Using Google Colab

1. **Open Google Colab:**

   Go to [Google Colab](https://colab.research.google.com/) and create a new notebook.

2. **Clone the Repository and Run Setup:**

   In a new cell, run:

   ```python
   # Clone the Cadock repository
   !git clone https://github.com/mehdikariim/DockCADD.git
   %cd DockCADD

   # Run the setup script
   !bash scripts/setup.sh
   from src.cadock import perform_docking

#Next
   # Define your SMILES list and PDB ID
   smiles_list = ["SMILES 1", "SMILES 2"]  # Replace with your SMILES
   PDB_ID = "PDB ID"  # Replace with your desired PDB ID

#Next
   # Perform docking
   perform_docking(smiles_list, PDB_ID)

#Next
   # Visualize results
   from src.cadock import visualize_results
   visualize_results(smiles_list, PDB_ID, 'docking_results')
