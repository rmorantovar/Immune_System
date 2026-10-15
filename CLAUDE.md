# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Research code for "How does the adaptive immune system learn?": stochastic simulations and analysis of naive/memory B- and T-cell repertoires responding to antigens (primary and memory responses, kinetic proofreading, epitope complexity, repertoire entropy, antigenic evolution, MHC presentation, TCR contact maps). There is no build system, package, or test suite — code is run as standalone scripts, notebooks, and Mathematica notebooks.

## Environment

- Python runs from the local venv: always use `.venv/bin/python` and `.venv/bin/pip`.
- Never install packages globally or with system Python.
- Packages the library needs: numpy, scipy, pandas, matplotlib, numba, seaborn, tqdm, joblib, and openpyxl (for `pd.read_excel`). `pickle5` is only imported on Python ≤ 3.7 and is not needed.

## Layout and how the pieces connect

- `Codes/library/` — shared Python helpers. Scripts pull in one of these modules with a star import, which also brings numpy, pandas, matplotlib, argparse, etc.:
  - `funcs.py` loads `funcs_mini.py`. Between them they hold energy models, sequence utilities, plotting colors, repertoire generation, expansions and ensembles of responses. Used by newer scripts (`from funcs import *`).
  - `functions_2.py` loads `functions_1.py`. Used by the older primary_response scripts and notebooks (`from functions_2 import *`).
  - `mf_lib.py` (mean-field engine, `run_simulation('semicomplete' | 'null' | 'complete' | 'approx', ...)`) and `mf_plotting.py` (plotting and colors).
  - `paths.py` defines `DROPBOX`, the local Dropbox folder (see Running things). `funcs_mini.py` and `functions_1.py` both import it, so it arrives through either star import. `mf_lib` users need `from paths import DROPBOX`.
- Library modules that no longer exist:
  - `functions.py` was deleted. Its replacement is `functions_2`, which has the same names, but `calculate_Q0` now normalises Q0 to 1 (the old version divided by 2× the sum, as a symmetry correction). Keep the normalisation to 1.
  - `lib_mf.py` is archived in `Codes/library/old_26_08/`. Use `mf_lib` + `mf_plotting`, which give numerically identical results. `run_simulation_semicomplete(...)` / `run_simulation_null(...)` became `run_simulation('semicomplete', ...)` / `run_simulation('null', ...)`.
  - `Immuno_models.py` (in `Codes/library/old/`) and `classes.py` (in `Codes/library/old_26_08/`) are only in archive folders. Scripts and notebooks that import them (e.g. `Jupyter_notebooks/primary_response/1_Dynamics.ipynb`, `Codes/python/primary_response/recognition.py`) don't run as-is.
- `Codes/python/<project>/` — simulation drivers, one folder per research project. `Codes/python/_general_simulations/MeRiM.py` is the general multi-epitope recurrent-infection simulator (primary + secondary/memory responses); its docstring explains `one_WT`, `secondary`, `secondary_all`, `N_inf`.
- `Codes/analysis/<project>/` and `Jupyter_notebooks/<project>/` — post-processing and figures that read simulation output.
- `Codes/cpp/primary_immune_response/` — older C++ simulators (GSL-based); `exe/` holds Python wrappers. `Codes/library/c/` and `Codes/library/randomc/` are C++ helpers/RNG.
- `Codes/in/` — input data: interaction energy matrices (`MJ.txt`, `TCRen.txt`, `BLOSUM62.txt`, alphabets), IEDB epitope tables, proteomes, contact maps. Selected via the `--energy_model` flag (default `TCRen`).
- `Codes/for_server/`, `sh/` — bash drivers that sweep parameters and call the Python scripts (e.g. `sh/memory_KPR.sh` runs MeRiM primary then secondary responses chained via `activated_repertoire.csv`).
- `mathematica/` — analytical work (`.nb`).
- Many `old*/` subfolders are archived versions; prefer the non-`old` file when editing.

## Running things

- Python env: `.venv/` at repo root (VS Code is configured to use it).
- Library imports rely on `sys.path.append(...)` with paths relative to the **current working directory**, so run scripts and notebooks from their own folder. The correct path depends on the folder's depth:
  - `'../../library/'` from `Codes/{analysis,python}/<project>/`
  - `'../../../library/'` one level deeper (e.g. `.../data/`, `.../N0/`)
  - `'../../Codes/library/'` from `Jupyter_notebooks/<project>/`

  Only `Codes/library` exists; `lib`, `my_lib` and `Codes/Python` do not. The exponential_proofreading and primary_response projects have been fixed (except archived `old*/` folders and files that need `Immuno_models`/`classes`). Other projects still have broken `'../library/'`, `'../../lib/'` and `'../../my_lib/'` lines, including `Codes/python/_general_simulations/MeRiM.py` (`'../../my_lib/'`).
- Simulations are configured via argparse flags, e.g.:
  ```
  cd Codes/python/_general_simulations
  python MeRiM.py --p 3.0 --pmem 1.0 --secondary 0 --pro memory_response --subpro Z_NA_dynamics
  ```
  Common flags: `--L0` (repertoire size), `--l` (sequence length), `--N_ens`, `--p` (proofreading steps), `--lamA`/`--lamB` (growth rates), `--E_lim`, `--t_lim`, `--chunk_size` (used when `L0 >= 1e7`), `--n_jobs`.
- Data and figures live in the author's Dropbox. Build paths from `DROPBOX` (from `Codes/library/paths.py`, which reads `~/.dropbox/info.json` and falls back to `~/Dropbox`), never from a hard-coded `/Users/<name>/...` prefix. Write `str(DROPBOX) + '/...'` or `f'{DROPBOX}/...'`. Current Dropbox layout:
  - Data and simulation output: `DROPBOX/Research_files/Immune_system/{project}/...`, e.g. `primary_response/simulations/output_N_ens_.../`, `primary_response/data/Victora_2016/`, `primary_response/data/victora_2020/` (lowercase `v`).
  - Energy matrices and alphabets (`TCRen.txt`, `Alphabet_TCRen.txt`, `MJ2.txt`, `BLOSUM62.txt`, `TCRen_potential.csv`): `DROPBOX/Research_files/Immune_system/in/`. Pass this folder to `get_motif`.
  - Figures and other documents: `DROPBOX/_Documents/Science/Projects/Immune_System/_Repository/Figures/{project}/...`.
  - Outdated layouts still appear in older scripts: `Research/Immune_system/`, `Research/immune_system/primary_immune_response/`, `Research/Evolution_Immune_System/Text_files/`, `My_Documents/...` and `_Documents/Research/Projects/...`. The exponential_proofreading and primary_response projects use the current layout, except for older primary_response scripts that read `Research/Evolution_Immune_System/Text_files/` (no longer in Dropbox; left as-is on purpose). Other projects still hard-code `/Users/robertomorantovar/...`.
- Simulation outputs are organized in parameter-named subfolders containing `activated_repertoire.csv`, `antigens.csv`, `params.json`. Existing outputs are skipped (scripts check `os.path.isfile` before running).
- C++ build (see `Codes/cpp/primary_immune_response/compiling.txt`):
  ```
  c++ --std=gnu++17 ../cpp/Primary_immune_response_local.cpp -lgsl -lgslcblas -o primary_immune_response_iMac.x
  ```
  On the cluster, add `-I/home/rmorantovar/gsl/include -L/home/rmorantovar/gsl/lib` and use the `*_server.cpp` variants.

## Gotchas

- Git tracks both `Codes/` and `codes/` (and `Jupyter_notebooks/` vs `jupyter_notebooks/`) as separate paths; on macOS's case-insensitive filesystem they collapse into one directory. Be careful when adding/moving files so they don't land under the wrong-case path in git.
- `funcs_mini.py` and `functions_1.py` set `plt.rcParams['text.usetex'] = True`, so plotting requires a LaTeX installation.
- Dropbox may keep files online-only (`ls -lO` shows `dataless`). Reading such a file blocks until Dropbox downloads it, and the download may not start from a sandboxed shell. Make the folder available offline in Finder before running scripts on it.
- Folder names are case-sensitive on Linux and the cluster even though macOS ignores case (e.g. `data/victora_2020` vs `Victora_2020`). Match the exact spelling on disk.
