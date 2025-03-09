#!/bin/bash

############################################
# scripts/setup.sh
# Installs all system and Python dependencies needed for DockcaddV2.
############################################

# 1) System Updates and Packages
sudo apt-get update -y
sudo apt-get install -y pymol openbabel wget tar openjdk-11-jdk

# 2) Python Libraries
echo "Upgrading pip and installing Python libraries..."
pip install --upgrade pip

# Install dependencies from requirements.txt if available
if [ -f "requirements.txt" ]; then
    pip install -r requirements.txt
fi

# Additional packages (if not listed in requirements.txt):
pip install git+https://github.com/openmm/pdbfixer.git
pip install openmm
pip install rdkit-pypi

echo "Python dependencies installed successfully."

# 3) AutoDock Vina Installation
VINA_URL="https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64"
VINA_NAME="vina_1.2.5_linux_x86_64"

if [ ! -f "/usr/local/bin/vina" ] && [ ! -f "/usr/local/bin/${VINA_NAME}" ]; then
    echo "Installing AutoDock Vina 1.2.5..."
    wget -q "$VINA_URL"
    chmod +x "$VINA_NAME"
    sudo mv "$VINA_NAME" /usr/local/bin/vina
    rm -f "$VINA_NAME"
    echo "AutoDock Vina installed successfully."
else
    echo "AutoDock Vina is already installed."
fi

# 4) p2rank Installation
P2RANK_URL="https://github.com/rdk/p2rank/releases/download/2.4.2/p2rank_2.4.2.tar.gz"
P2RANK_DIR="p2rank_2.4.2"

if [ ! -d "$P2RANK_DIR" ]; then
    echo "Installing p2rank 2.4.2..."
    wget -q "$P2RANK_URL"
    tar -xzf "p2rank_2.4.2.tar.gz"
    rm -f "p2rank_2.4.2.tar.gz"
    echo "p2rank installed successfully."
else
    echo "p2rank is already installed."
fi

echo "Setup is complete."
