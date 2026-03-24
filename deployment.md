# Getting started

Step-by-step instructions to set up the Fanal exercise on your computer.
Instructions are provided for **macOS**, **Linux**, and **Windows**.

---

## 1. Install prerequisites

You need **Git**, **Python 3.11+**, and **Jupyter**.
The recommended way is to install [Miniconda](https://docs.conda.io/en/latest/miniconda.html), which works on all three platforms.

| Platform | Install Miniconda |
|----------|-------------------|
| macOS | `brew install miniconda` or download from [miniconda.io](https://docs.conda.io/en/latest/miniconda.html) |
| Linux | `wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh && bash Miniconda3-latest-Linux-x86_64.sh` |
| Windows | Download and run the installer from [miniconda.io](https://docs.conda.io/en/latest/miniconda.html). Use **Anaconda Prompt** for the commands below. |

---

## 2. Clone the repository

Open a terminal (or Anaconda Prompt on Windows) and run:

```bash
git clone https://github.com/jahernando/USC-FPII-Fanal.git
cd USC-FPII-Fanal
```

---

## 3. Create the Python environment

```bash
conda env create -f environment.yml
conda activate fanal
```

This installs Python, Jupyter, NumPy, pandas, SciPy, Matplotlib, and PyTables.

**Alternative (pip only, no conda):**

```bash
python -m venv fanal-env

# macOS / Linux:
source fanal-env/bin/activate

# Windows:
fanal-env\Scripts\activate

pip install -r requirements.txt
```

---

## 4. Configure the environment

The project needs to know where the code is. Run the setup script **every time you open a new terminal**:

**macOS / Linux:**

```bash
source setup.sh
```

**Windows (Command Prompt):**

```cmd
set FANAL_ROOT=%cd%
set PYTHONPATH=%FANAL_ROOT%;%PYTHONPATH%
```

**Windows (PowerShell):**

```powershell
$env:FANAL_ROOT = (Get-Location).Path
$env:PYTHONPATH = "$env:FANAL_ROOT;$env:PYTHONPATH"
```

> **Tip:** On Windows, you can create a `setup.bat` file with the two `set` lines above and run it with `setup.bat`.

---

## 5. Download the data

The instructor will share a OneDrive link with the data files.

1. Download the `.h5` file assigned to your collaboration (e.g. `fanal_new_gamma.h5`).
2. Place it in the `data/` directory inside the repository:

```
USC-FPII-Fanal/
├── data/
│   └── fanal_new_gamma.h5    <-- your file goes here
├── notebooks/
├── ana/
└── ...
```

If the `data/` directory does not exist, create it:

```bash
mkdir data
```

Then move or copy the downloaded file into it.

---

## 6. Verify the setup

Start Jupyter:

```bash
jupyter notebook
```

Open `notebooks/guide/crib_pandas.ipynb` and run all cells. If everything is correct you should see:

- `Fanal root: /path/to/USC-FPII-Fanal` (your actual path)
- `Loaded XXXXX events` (a number in the tens of thousands)
- Several histogram plots

If you get an error:

| Error | Solution |
|-------|----------|
| `ModuleNotFoundError: No module named 'tables'` | `pip install tables` or `conda install tables` |
| `ModuleNotFoundError: No module named 'core'` | You forgot to run `source setup.sh` (or the Windows equivalent) |
| `FileNotFoundError: ... .h5` | The data file is not in the `data/` directory, or the filename does not match |
| `No module named 'pandas'` | `pip install pandas` or `conda install pandas` |

---

## 7. Start working

The exercise notebooks are in `notebooks/guide/` and should be completed in this order:

1. `fanal_selection.ipynb`
2. `fanal_energy_resolution.ipynb`
3. `fanal_bkg.ipynb`
4. `fanal_bkg_uncertainties.ipynb`
5. `fanal_signal_countexp.ipynb`
6. `fanal_signal.ipynb`
7. `fanal_data_access.ipynb` (open the box)

Reference material (crib sheets and course guides) is available in the same directory.

---

## Quick reference

```bash
# Every time you open a new terminal:
cd USC-FPII-Fanal
conda activate fanal          # or: source fanal-env/bin/activate
source setup.sh               # macOS/Linux only; on Windows use set/env commands
jupyter notebook
```
