# Teacher Workflow

## Repository structure

| Branch | Audience | Content of `notebooks/guide/` |
|---|---|---|
| `teacher` | Professor | Full solutions, cells tagged with `solution` |
| `main` | Students | Exercise cells replaced with `# YOUR CODE HERE` |

Students clone `main` and see notebooks with blank exercise cells.
The professor works on `teacher` and uses the scripts below to publish updates.

---

## Day-to-day: edit a notebook

```bash
git checkout teacher
# open JupyterLab and edit notebooks/guide/ as usual
git add notebooks/guide/<notebook>.ipynb
git commit -m "description of change"
```

---

## Publish student version

After committing your changes on `teacher`, run:

```bash
./deploy_to_main.sh
```

This script automatically:
1. Regenerates `notebooks/student/` (strips solution cells → `# YOUR CODE HERE`)
2. Copies the result to `notebooks/guide/` on `main`
3. Commits and pushes `main`
4. Returns to the `teacher` branch

---

## Mark a cell as an exercise (solution cell)

In **JupyterLab**:
1. Select the code cell that contains the solution.
2. Open the **Property Inspector** panel (⚙️ icon, right sidebar).
3. Under **Cell Metadata → Tags**, type `solution` and press Enter.
4. Save the notebook, commit, and run `./deploy_to_main.sh`.

In the **classic Jupyter Notebook** interface:
- View → Cell Toolbar → Tags → type `solution` → Add tag.

---

## Remove the exercise tag from a cell

Follow the same steps as above and delete the `solution` tag.
Then commit and run `./deploy_to_main.sh`.

---

## Notebooks and their exercise cells

The following table lists which cells are currently tagged as `solution`
in each notebook (0-based index). Update this table when you add or remove tags.

| Notebook | Tagged cells | What the student implements |
|---|---|---|
| `fanal_def.ipynb` | 16, 24, 25, 29, 30 | Variable histograms; efficiency computation and plot; energy resolution fit |
| `fanal_selection.ipynb` | 30, 31 | Efficiency computation and plot |
| `fanal_energy_resolution.ipynb` | 16, 17 | Gaussian fit to the Bi/Tl peaks; print energy resolution |
| `fanal_bkg.ipynb` | 25, 26 | Number of background events in E range and RoI |
| `fanal_bkg_uncertanties.ipynb` | 23 | Uncertainty propagation |
| `fanal_signal.ipynb` | 17, 19, 30 | `nevts_total()` function; experiment simulation; half-life |
| `fanal_signal_countexp.ipynb` | 17, 18, 20, 21, 22 | CI computation; `tau()` function; plots |

---

## Push both branches to GitHub

First time only:

```bash
git checkout teacher
git push -u origin teacher

git checkout main
git push -u origin main

git checkout teacher
```

Subsequent pushes of the `teacher` branch (student pushes are handled by `deploy_to_main.sh`):

```bash
git checkout teacher
# ... edit, commit ...
git push origin teacher
```

---

## Key files (teacher branch only)

| File | Purpose |
|---|---|
| `make_student_nbs.py` | Generates `notebooks/student/` from `notebooks/guide/` |
| `deploy_to_main.sh` | End-to-end script: generate → copy → commit → push main |
| `TEACHER_WORKFLOW.md` | This file |
