#!/bin/bash

# Update and install dependencies
sudo apt-get update -y
sudo apt-get install -y pymol openbabel wget tar openjdk-11-jdk

# Install Python libraries
pip install --upgrade pip
pip install -r requirements.txt

# Install AutoDock Vina
VINA_URL="https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64"
VINA_NAME="vina_1.2.5_linux_x86_64"

if [ ! -f "/usr/local/bin/vina" ]; then
    wget $VINA_URL
    chmod +x $VINA_NAME
    sudo mv $VINA_NAME /usr/local/bin/vina
    rm -f vina_1.2.5_linux_x86_64
    echo "AutoDock Vina installed successfully."
else
    echo "AutoDock Vina is already installed."
fi

# Install p2rank
P2RANK_URL="https://github.com/rdk/p2rank/releases/download/2.4.2/p2rank_2.4.2.tar.gz"
P2RANK_DIR="p2rank_2.4.2"

if [ ! -d "$P2RANK_DIR" ]; then
    wget $P2RANK_URL
    tar -xzf p2rank_2.4.2.tar.gz
    rm -f p2rank_2.4.2.tar.gz
    echo "p2rank installed successfully."
else
    echo "p2rank is already installed."
fi
