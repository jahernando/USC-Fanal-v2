#!/usr/bin/env bash
# deploy_to_main.sh
# =================
# Generate student notebooks from the teacher branch and push them to main.
#
# Run this script from the repo root while on the 'teacher' branch.
#
# Usage:
#   ./deploy_to_main.sh
#
# What it does:
#   1. Regenerates notebooks/student/ from notebooks/guide/ (using solution tags).
#   2. Switches to main, copies the student notebooks over guide/.
#   3. Commits and pushes main.
#   4. Returns to the teacher branch.

set -e  # exit on any error

CURRENT_BRANCH=$(git rev-parse --abbrev-ref HEAD)

if [ "$CURRENT_BRANCH" != "teacher" ]; then
    echo "ERROR: Must be on the 'teacher' branch. Currently on '$CURRENT_BRANCH'."
    exit 1
fi

# Check there are no uncommitted changes on teacher
if ! git diff --quiet HEAD; then
    echo "ERROR: You have uncommitted changes on the teacher branch. Commit them first."
    exit 1
fi

echo "==> Generating student notebooks..."
python make_student_nbs.py

echo "==> Switching to main..."
git checkout main

echo "==> Copying student notebooks to notebooks/guide/..."
cp notebooks/student/*.ipynb notebooks/guide/

echo "==> Staging student notebooks..."
git add notebooks/guide/*.ipynb

echo "==> Committing to main..."
git commit -m "student: update notebooks from teacher branch $(date +%Y-%m-%d)"

echo "==> Pushing main to origin..."
git push origin main

echo "==> Returning to teacher branch..."
git checkout teacher

echo ""
echo "Done. Students will see the updated notebooks on the main branch."
