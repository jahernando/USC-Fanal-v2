# Student setup checklist

Use this checklist to verify that the exercise works end-to-end before distributing it to students.

---

## 1. Clone and setup

- [ ] Clone the repo: `git clone https://github.com/jahernando/USC-FPII-Fanal.git`
- [ ] `cd USC-FPII-Fanal`
- [ ] Source the environment: `source setup.sh`
- [ ] Verify: `echo $FANAL_ROOT` prints the repo root
- [ ] Verify: `python -c "import numpy, pandas, matplotlib, scipy, tables; print('OK')"` — if `tables` is missing, install it: `pip install tables`

## 2. Verify data access

- [ ] Check data files exist: `ls data/fanal_new_*.h5` (should list 5 files: new_alpha, new_beta, new_gamma, new_delta, new_epsilon)
- [ ] Open `notebooks/guide/crib_pandas.ipynb` in Jupyter and run all cells — should produce histograms without errors

## 3. Run the exercise notebooks in order

For each notebook below, open it in Jupyter, run the pre-filled cells, and complete the `# YOUR CODE HERE` exercises.

### Session 1: Selection
- [ ] `notebooks/guide/fanal_selection.ipynb`
  - [ ] Selection cuts defined, efficiency plots generated
  - [ ] `display_collpars()` table shows reasonable values
  - [ ] `collpars.py` created with selection parameters and efficiencies

### Session 2: Energy resolution
- [ ] `notebooks/guide/fanal_energy_resolution.ipynb`
  - [ ] Gaussian+line fit converges, sigma extracted
  - [ ] `display_collpars()` table shows sigma values
  - [ ] `collpars.py` updated with sigma values

### Session 3: Background estimation
- [ ] `notebooks/guide/fanal_bkg.ipynb`
  - [ ] Blind selection efficiency computed
  - [ ] Extended likelihood fit to blind data converges
  - [ ] Event counts (n_Bi, n_Tl) for total, blind, E-range, and RoI computed
  - [ ] `display_collpars()` table shows all values
  - [ ] `collpars.py` updated with background counts

### Session 4: Uncertainties
- [ ] `notebooks/guide/fanal_bkg_uncertainties.ipynb`
  - [ ] Profile likelihood scan produces parabolic curves
  - [ ] 68% CL intervals extracted
  - [ ] `display_collpars()` table shows uncertainties
  - [ ] `collpars.py` updated with un_Bi_blind, un_Tl_blind

### Session 5: Signal
- [ ] `notebooks/guide/fanal_signal_countexp.ipynb`
  - [ ] Feldman-Cousins confidence belt constructed
  - [ ] Upper limit on signal count obtained
- [ ] `notebooks/guide/fanal_signal.ipynb`
  - [ ] 3-component extended likelihood fit runs
  - [ ] Profile likelihood scan for signal parameter
  - [ ] q0 test statistic and discovery significance computed
  - [ ] Half-life limit or measurement obtained

### Final: Open the box
- [ ] `notebooks/guide/fanal_data_access.ipynb`
  - [ ] Full data (blind + RoI) analysed
  - [ ] Final result on signal (limit or evidence)

## 4. Verify reference material

- [ ] `notebooks/guide/crib_pandas.ipynb` — runs without errors
- [ ] `notebooks/guide/crib_fit.ipynb` — runs without errors (uses blind data)
- [ ] `notebooks/guide/guide_dynamics.md` — renders correctly in Jupyter / book
- [ ] `notebooks/guide/guide_presentation.md` — renders correctly, evaluation criteria visible

## 5. Verify collpars.py at the end

- [ ] `collpars.py` in `notebooks/guide/` contains all expected variables:
  ```
  collaboration, sel_ntracks, sel_eblob2, sel_erange, sel_eroi
  eff_Bi_E, eff_Bi_RoI, eff_Tl_E, eff_Tl_RoI, eff_bb_E, eff_bb_RoI
  sigma_Bi (and optionally sigma_Tl, sigma_bb)
  exposure, acc_bb
  eff_Bi_blind, eff_Tl_blind
  n_Bi_total, n_Bi_blind, n_Bi_E, n_Bi_RoI
  n_Tl_total, n_Tl_blind, n_Tl_E, n_Tl_RoI
  un_Bi_blind, un_Tl_blind
  ```
- [ ] Values are physically reasonable (not NaN, not zero, correct order of magnitude)

## 6. Clean up before distributing

- [ ] Delete `collpars.py` from `notebooks/guide/` (students must generate their own)
- [ ] Clear all notebook outputs: `jupyter nbconvert --clear-output --inplace notebooks/guide/*.ipynb`
- [ ] Verify that `TEACHER_WORKFLOW.md`, `FANAL_IMPROVEMENTS.md`, and solution notebooks (`*_sols.ipynb`) are **not** on the `main` branch

## Notes

- Students should use `main` branch (no solutions). Teacher uses `teacher` branch.
- Each collaboration gets assigned one dataset (new_alpha through new_epsilon).
- The Jupyter Book is published at https://jahernando.github.io/USC-FPII-Fanal/
