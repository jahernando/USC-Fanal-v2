# Fanal - Improvement Plan

Review of the student-facing notebooks (main branch) with proposed improvements.

## Quick wins

- [x] **Rename `fanal_bkg_uncertanties.ipynb`** → `fanal_bkg_uncertainties.ipynb`. Update `_toc.yml` and any internal links.
- [x] **Move `fanal_true.ipynb`** to `teacher` branch only (already absent from `main`).
- [x] **Replace `os.getcwd()[:-15]` path hack** in all 7 guide notebooks. Now uses `Path.cwd()` walking up to find `data/` + `ana/`, with `FANAL_ROOT` / `USCFANALDIR` env-var support. Updated `setup.sh` to export `FANAL_ROOT`.

## Pedagogical improvements

- [x] **Add a notation-to-code table** in all notebooks that write to `collpars.py` (`fanal_selection`, `fanal_energy_resolution`, `fanal_bkg`, `fanal_bkg_uncertainties`). Markdown reference table + `display_collpars()` function (in `ana/fanal_display.py`) connecting LaTeX symbols to Python variable names with computed values.
- [ ] **Add sanity checks** after exercise cells. E.g. `assert 10 < n_Bi_total < 50000` so students catch errors before they propagate through `collpars.py`.
- [ ] **Add "expected output" guidance** in the harder notebooks. Not the solution, but the expected format/order of magnitude: "you should see a plot with two curves", "n_Bi_total should be O(600)".
- [x] **Smooth the difficulty curve** in `fanal_signal_countexp.ipynb` and `fanal_signal.ipynb`: worked examples added before exercises (FC belt with nbkg=0, tau function, full generate→fit→profile pipeline with large signal).

## Structural improvements

- [x] **Flesh out `fanal_data_access.ipynb`** ("open the box"). Added: "stop and think" planning cell, progressive selection plots, counting experiment with FC CI, extended likelihood fit, profile likelihood, p-value, half-life estimation, summary table.
- [x] **Add a crib sheet notebook** `crib_pandas.ipynb` covering: HDF5 data reading, pandas DataFrame basics, boolean masks/selections, `pltext.hist` and `pltext.canvas` usage.
- [x] **Add a crib sheet notebook** `crib_fit.ipynb` covering: `hfit` (Gaussian+line fits), `efit` (extended likelihood fits), profile likelihood scans, confidence intervals. Based on patterns already used in the exercise notebooks.
- [x] **Add a presentation/article guide** (`guide_presentation.md`): deliverables, expected content, evaluation criteria (50/30/20). Plus course dynamics guide (`guide_dynamics.md`): session structure, notebook sequence, cyclic presentation scheme.
- [ ] **Consider `nbgrader`** for automated validation of student solutions.

## Already done (2026-03-23)

- [x] Fixed English spelling and grammar across all 8 notebooks (15 fixes in `fanal.ipynb`, typos in guide notebooks).
- [x] Fixed `^{108}Tl` → `^{208}Tl` in `fanal.ipynb`.
- [x] Fixed step numbering (duplicate step 5) in `fanal.ipynb`.
- [x] Added scaffolding hints to all `# YOUR CODE HERE` cells across all guide notebooks.
- [x] Reorganised branch workflow: guide NBs identical on main and teacher; solutions in separate `*_sols.ipynb` (teacher-only). `git merge main` into teacher is now conflict-free.
- [x] Added worked examples in `fanal_signal_countexp.ipynb` (FC belt nbkg=0, tau function) and `fanal_signal.ipynb` (full generate→fit→profile with large signal).

## Already done (2026-03-24)

- [x] Added notation-to-code tables (markdown + `display_collpars()` in `ana/fanal_display.py`) in all 4 notebooks that write to `collpars.py`.
- [x] Added `.fillna(0.)` to all `pd.read_hdf()` calls across 14 notebooks. NaN values in HDF5 data (blob/track columns) were generation artifacts; filling with 0 does not change analysis results.
- [x] Added `crib_pandas.ipynb` crib sheet notebook with HDF5/pandas/masks/pltext reference.
- [x] Added `crib_fit.ipynb` crib sheet notebook with hfit/efit/profile-likelihood/signal-hypothesis reference.
- [x] Added `guide_presentation.md` (deliverables, evaluation rubric) and `guide_dynamics.md` (course workflow, session structure).
- [x] Organised `_toc.yml` into three sections: Exercise notebooks, Reference material, Course structure.
- [x] Set up GitHub Pages deployment via GitHub Actions (`deploy-book.yml`). Book at https://jahernando.github.io/USC-FPII-Fanal/.
- [x] Updated README, intro.md with Claude acknowledgement. Tagged v2026.1.0.
- [x] Fixed NEXT_White_events.png not rendering in book (moved to `_static/`).
