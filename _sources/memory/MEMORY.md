# FANAL-v2 Memory

## Project overview
Physics analysis for the NEXT experiment (0νββ search in ¹³⁶Xe).
Working directory: `/Users/hernando/work/docencia/master/FPII/USC-Fanal-v2`

## Code structure
- `core/` – reusable physics/stats utilities (no NEXT-specific content)
  - `utils.py`   – array/DataFrame helpers (selection, efficiency, stats)
  - `hfit.py`    – histogram fitting (gaus, line, exp, composites)
  - `efit.py`    – extended likelihood (ComPDF, ExtComPDF, SimulExtComPDF)
  - `confint.py` – Feldman-Cousins confidence intervals
  - `pltext.py`  – matplotlib extensions
  - `fc_confint.py` – OLD draft (missing imports, typos) – superseded by confint.py
- `ana/` – NEXT/FANAL specific analysis
  - `fanal.py`     – main analysis (ELL fit, toy MC, test statistics)
  - `fanal_gen.py` – pseudo-data generation from MC templates
  - `pltfanal.py`  – analysis plots
  - `collpars.py`  – collaboration parameter constants (namedtuple Measurement)

## Refactor done (2026-03-03)
All files in core/ and ana/ were reviewed and improved:
- Added module-level docstrings to all files
- Improved all function docstrings (NumPy style)
- Fixed bugs: `mask is not None` (was `mask != None`), shadowed builtins (`range`→`val_range`, `int`→`included`), mutable default args (`[]`→`None`)
- `fanal_gen.py`: removed hardcoded path, now requires `dirpath` argument
- `pltext.py`: `str_stats` parameter renamed `formate`→`fmt`; same in `utils.py`
- Removed all large commented-out code blocks
- `fc_confint.py`: NOT modified (legacy draft with missing imports – do not use)

## Two-version notebook workflow (2026-03-03)
- `notebooks/guide/`   → professor version (source of truth, full solutions)
- `notebooks/student/` → student version (auto-generated, solution cells → `# YOUR CODE HERE`)
- Script: `make_student_nbs.py` (run from repo root) regenerates student/ from guide/
- Cell tagging: cells tagged `"solution"` in guide notebooks are stripped in the student version
  - In JupyterLab: Select cell → Property Inspector (⚙) → Tags → add `solution`
- Tagged cells per notebook (0-based indices):
  - fanal_def.ipynb: 16, 24, 25, 29, 30
  - fanal_selection.ipynb: 30, 31
  - fanal_energy_resolution.ipynb: 16, 17
  - fanal_bkg.ipynb: 25, 26
  - fanal_bkg_uncertanties.ipynb: 23
  - fanal_signal.ipynb: 17, 19, 30
  - fanal_signal_countexp.ipynb: 17, 18, 20, 21, 22

## Key API notes
- `ut.str_stats(vals, val_range=None, fmt='6.2f')` – note `val_range` and `fmt`
- `pltext.hist(...)` uses `stats_format` karg (pops it internally)
- `fanal_gen.experiment()` now requires explicit `dirpath` keyword argument
- `prepare_fit_ell` / `prepare_fit_simell`: `refnames`/`refranges` default to `None` (not `[]`)
