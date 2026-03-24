#!/bin/bash
# Run all solution notebooks in sequence to verify the full analysis flow.
# Must be executed from the repo root with the fanal conda env active.
#
# Usage:
#   source setup.sh
#   bash run_sols.sh

set -e

cd notebooks/guide

# Start fresh
rm -f collpars.py

echo "=== 1/5 fanal_selection_sols ==="
jupyter nbconvert --to notebook --execute fanal_selection_sols.ipynb

echo "=== 2/5 fanal_energy_resolution_sols ==="
jupyter nbconvert --to notebook --execute fanal_energy_resolution_sols.ipynb

echo "=== 3/5 fanal_bkg_sols ==="
jupyter nbconvert --to notebook --execute fanal_bkg_sols.ipynb

echo "=== 4/5 fanal_bkg_uncertainties_sols ==="
jupyter nbconvert --to notebook --execute fanal_bkg_uncertainties_sols.ipynb

echo "=== 5/5 fanal_data_access_sols ==="
jupyter nbconvert --to notebook --execute fanal_data_access_sols.ipynb

# Clean up generated files
rm -f *_sols.nbconvert.ipynb

echo ""
echo "All solution notebooks executed successfully."
echo "Generated collpars.py:"
cat collpars.py
