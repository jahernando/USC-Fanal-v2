"""
make_student_nbs.py
===================
Generate student-version notebooks from the professor (guide) notebooks.

Workflow
--------
1. Read each .ipynb file from notebooks/guide/.
2. Any code cell whose metadata contains the tag ``"solution"`` is replaced
   with an empty code cell containing just ``# YOUR CODE HERE``.
3. The modified notebook is written to notebooks/student/.

Usage
-----
    python make_student_nbs.py

Adding solution tags in JupyterLab
-----------------------------------
Select a cell → right-click → "Add Cell Tag" → type ``solution`` → Enter.
Or open the Property Inspector panel (gear icon) → "Cell Metadata" → add::

    {"tags": ["solution"]}

"""

import json
import os
import shutil

GUIDE_DIR   = os.path.join('notebooks', 'guide')
STUDENT_DIR = os.path.join('notebooks', 'student')
STUB        = '# YOUR CODE HERE\n'


def is_solution(cell):
    """Return True if the cell has the 'solution' tag."""
    return 'solution' in cell.get('metadata', {}).get('tags', [])


def make_stub_cell(cell):
    """Return a copy of *cell* with source replaced by the stub comment."""
    stub = json.loads(json.dumps(cell))   # deep copy
    stub['source'] = [STUB]
    stub['outputs'] = []
    stub['execution_count'] = None
    return stub


def process_notebook(src_path, dst_path):
    with open(src_path) as f:
        nb = json.load(f)

    n_replaced = 0
    for i, cell in enumerate(nb['cells']):
        if cell['cell_type'] == 'code' and is_solution(cell):
            nb['cells'][i] = make_stub_cell(cell)
            n_replaced += 1

    with open(dst_path, 'w') as f:
        json.dump(nb, f, indent=1, ensure_ascii=False)

    return n_replaced


def main():
    os.makedirs(STUDENT_DIR, exist_ok=True)

    notebooks = sorted(
        f for f in os.listdir(GUIDE_DIR) if f.endswith('.ipynb')
    )

    if not notebooks:
        print(f'No notebooks found in {GUIDE_DIR}/')
        return

    print(f'Generating student notebooks in {STUDENT_DIR}/\n')
    for nb_name in notebooks:
        src = os.path.join(GUIDE_DIR,   nb_name)
        dst = os.path.join(STUDENT_DIR, nb_name)
        n = process_notebook(src, dst)
        print(f'  {nb_name:45s}  {n:2d} solution cell(s) replaced')

    # Copy any non-notebook files (e.g. collpars.py) that students also need
    for fname in os.listdir(GUIDE_DIR):
        if not fname.endswith('.ipynb') and not fname.startswith('__'):
            src = os.path.join(GUIDE_DIR, fname)
            dst = os.path.join(STUDENT_DIR, fname)
            if os.path.isfile(src):
                shutil.copy2(src, dst)
                print(f'  {fname:45s}  (copied)')

    print('\nDone.')


if __name__ == '__main__':
    main()
