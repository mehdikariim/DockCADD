#!/bin/bash
# Install AutoDock Vina
wget https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v1.2.5/vina_1.2.5_linux_x86_64
chmod +x vina_1.2.5_linux_x86_64
mv vina_1.2.5_linux_x86_64 /usr/local/bin/vina

# Install P2Rank
wget https://github.com/rdk/p2rank/releases/download/2.4.2/p2rank_2.4.2.tar.gz
tar -xzf p2rank_2.4.2.tar.gz

echo "Dependencies installed successfully."
