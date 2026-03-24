#!/bin/bash
# Run all solution notebooks in sequence to verify the full analysis flow.
# Must be executed from the repo root with the fanal conda env active.
#
# Usage:
#   source setup.sh
#   bash run_sols.sh [collaboration]
#
# If no collaboration is given, defaults to "new_gamma".
# Available collaborations: new_alpha, new_beta, new_gamma, new_delta, new_epsilon
#
# Output: notebooks/guide/collpars_<collaboration>.py

set -e

COLL="${1:-new_gamma}"
export FANAL_COLLABORATION="$COLL"

echo "=== Running full analysis for collaboration: $COLL ==="
echo ""

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

# Rename to collaboration-specific file
mv collpars.py "collpars_${COLL}.py"

echo ""
echo "All solution notebooks executed successfully for collaboration: $COLL"
echo "Generated collpars_${COLL}.py:"
cat "collpars_${COLL}.py"
