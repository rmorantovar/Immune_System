# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

Research code for "How does the adaptive immune system learn?": stochastic simulations and analysis of naive/memory B- and T-cell repertoires responding to antigens (primary and memory responses, kinetic proofreading, epitope complexity, repertoire entropy, antigenic evolution, MHC presentation, TCR contact maps). There is no build system, package, or test suite — code is run as standalone scripts, notebooks, and Mathematica notebooks.

## Environment

- Python runs from the local venv: always use `.venv/bin/python` and `.venv/bin/pip`.
- Never install packages globally or with system Python.

## Layout and how the pieces connect

- `Codes/library/` — shared Python helpers. `funcs_mini.py` (energy models, sequence utilities, plotting colors/constants, imports everything via `*`) is imported by `funcs.py` (repertoire generation, expansions, ensembles of responses). Other modules: `functions_1.py`, `functions_2.py`, `mf_lib.py`/`mf_plotting.py` (mean-field). Scripts use `from funcs import *`, so numpy/pandas/matplotlib/argparse/etc. come in through this star import.
- `Codes/python/<project>/` — simulation drivers, one folder per research project. `Codes/python/_general_simulations/MeRiM.py` is the general multi-epitope recurrent-infection simulator (primary + secondary/memory responses); its docstring explains `one_WT`, `secondary`, `secondary_all`, `N_inf`.
- `Codes/analysis/<project>/` and `Jupyter_notebooks/<project>/` — post-processing and figures that read simulation output.
- `Codes/cpp/primary_immune_response/` — older C++ simulators (GSL-based); `exe/` holds Python wrappers. `Codes/library/c/` and `Codes/library/randomc/` are C++ helpers/RNG.
- `Codes/in/` — input data: interaction energy matrices (`MJ.txt`, `TCRen.txt`, `BLOSUM62.txt`, alphabets), IEDB epitope tables, proteomes, contact maps. Selected via the `--energy_model` flag (default `TCRen`).
- `Codes/for_server/`, `sh/` — bash drivers that sweep parameters and call the Python scripts (e.g. `sh/memory_KPR.sh` runs MeRiM primary then secondary responses chained via `activated_repertoire.csv`).
- `mathematica/` — analytical work (`.nb`).
- Many `old*/` subfolders are archived versions; prefer the non-`old` file when editing.

## Running things

- Python env: `.venv/` at repo root (VS Code is configured to use it).
- Library imports rely on `sys.path.append(...)` with paths **relative to the script's directory**, and these are inconsistent (`'../library/'`, `'../../lib/'`, `'../../my_lib/'`, …; `Codes/my_lib` and `lib` do not currently exist, `Codes/library` does). Run scripts from their own directory, and fix the `sys.path` line to point at `Codes/library` if imports fail.
- Simulations are configured via argparse flags, e.g.:
  ```
  cd Codes/python/_general_simulations
  python MeRiM.py --p 3.0 --pmem 1.0 --secondary 0 --pro memory_response --subpro Z_NA_dynamics
  ```
  Common flags: `--L0` (repertoire size), `--l` (sequence length), `--N_ens`, `--p` (proofreading steps), `--lamA`/`--lamB` (growth rates), `--E_lim`, `--t_lim`, `--chunk_size` (used when `L0 >= 1e7`), `--n_jobs`.
- Output paths are **hardcoded absolute paths** to the author's Dropbox (e.g. `/Users/robertomorantovar/Dropbox/Research/Immune_system/{project}/{subproject}/...`), organized in parameter-named subfolders containing `activated_repertoire.csv`, `antigens.csv`, `params.json`. Existing outputs are skipped (scripts check `os.path.isfile` before running). Note the home directory in these paths may differ from the current machine's.
- C++ build (see `Codes/cpp/primary_immune_response/compiling.txt`):
  ```
  c++ --std=gnu++17 ../cpp/Primary_immune_response_local.cpp -lgsl -lgslcblas -o primary_immune_response_iMac.x
  ```
  On the cluster, add `-I/home/rmorantovar/gsl/include -L/home/rmorantovar/gsl/lib` and use the `*_server.cpp` variants.

## Gotchas

- Git tracks both `Codes/` and `codes/` (and `Jupyter_notebooks/` vs `jupyter_notebooks/`) as separate paths; on macOS's case-insensitive filesystem they collapse into one directory. Be careful when adding/moving files so they don't land under the wrong-case path in git.
- `funcs_mini.py` sets `plt.rcParams['text.usetex'] = True`, so plotting requires a LaTeX installation.
