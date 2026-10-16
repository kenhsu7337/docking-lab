#!/usr/bin/env bash
# Runs once when the codespace is created.
set -euo pipefail

VINA_VERSION=1.2.7

echo "Creating the docking environment (this takes several minutes)..."
conda env create -f environment.yml || conda env update -f environment.yml -n docking
conda run -n docking python -m ipykernel install --user --name docking --display-name "docking"

echo "Installing the Vina command-line program..."
curl -fsSL -o /tmp/vina \
  "https://github.com/ccsb-scripps/AutoDock-Vina/releases/download/v${VINA_VERSION}/vina_${VINA_VERSION}_linux_x86_64"
sudo install -m 0755 /tmp/vina /usr/local/bin/vina

# Every new terminal starts in the docking environment.
grep -q "conda activate docking" ~/.bashrc || echo "conda activate docking" >> ~/.bashrc

vina --version
conda run -n docking python -c "import vina, meeko, rdkit; print('Environment ready')"
