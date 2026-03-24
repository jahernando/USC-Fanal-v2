# Course dynamics

Description of how the Fanal exercise is organised: group formation, class sessions, results discussion, presentations, and deliverables.

---

## 1. Formation of collaborations

- Students form **collaborations** of 2–3 people, by their own choice.
- The instructor assigns a **dataset** (new_alpha, new_beta, new_gamma, new_delta, new_epsilon) to each collaboration. There are approximately 5 collaborations.

---

## 2. Analysis sessions (5 class sessions)

Each session follows the same structure:

1. **Introduction by the instructor** — brief presentation of the physics and analysis techniques covered in the notebook(s) of the day.
2. **Group work** — each collaboration works through 1–2 notebooks. The instructor moves between groups, answering questions and resolving doubts.
3. **Wrap-up** — the session ends when each collaboration has obtained the relevant parameters and written them into `collpars.py`.

### Notebook sequence

| Session | Notebooks | Goal |
|---------|-----------|------|
| 1 | `fanal_selection` | Define selection cuts, compute efficiencies |
| 2 | `fanal_energy_resolution` | Measure energy resolution |
| 3 | `fanal_bkg` | Estimate background event counts |
| 4 | `fanal_bkg_uncertainties` | Compute uncertainties via profile likelihood |
| 5 | `fanal_signal_countexp`, `fanal_signal` | Signal search: counting experiment, extended likelihood, hypothesis testing |

### Reference material

The following support materials are available throughout the exercise:

- `crib_pandas` — quick reference for data access, pandas, masks, and plotting.
- `crib_fit` — quick reference for histogram fits, extended likelihood, and profile likelihood scans.
- `guide_presentation` — deliverables, expected content, and evaluation criteria.

---

## 3. Inter-collaboration discussion

After all collaborations have obtained their results (end of session 5), a dedicated time slot is used for an **inter-collaboration meeting**:

- A representative from each collaboration meets with the others to share and compare results.
- Key questions to discuss: are the results from different datasets consistent? What happens when you consider them together? Can the combined sensitivity improve the limit or the evidence?
- This step mirrors real particle physics practice, where independent experiments combine results.

---

## 4. Presentations (2 class sessions)

Presentations follow a **cyclic scheme**: one collaboration presents, another asks questions, then they rotate.

| Step | Presents | Questions from |
|------|----------|----------------|
| 1 | A | B |
| 2 | B | C (A leaves) |
| 3 | C | D (B leaves) |
| 4 | D | E (C leaves) |
| 5 | E | A (D leaves) |

- Each presentation lasts **12 minutes**, followed by **8 minutes** for questions.
- Questions come from the evaluating collaboration **and** from the instructor, who is present at all presentations.
- The evaluating collaboration's questions count towards their own *Understanding* grade (see evaluation criteria in the presentation guide).

---

## 5. Deliverables

After the presentations, each collaboration submits (as PDF) to the virtual classroom:

1. **Presentation slides**
2. **Written article** (4–5 pages including figures)
3. **Methodology report**

The submission deadline is **after** the presentation sessions (single deadline for all three documents).

---

## Timeline summary

```
Sessions 1–5     Analysis work in class (1–2 notebooks per session)
After session 5   Inter-collaboration discussion
Sessions 6–7      Presentations (cyclic scheme)
After sessions    Submission of deliverables (single deadline)
```
