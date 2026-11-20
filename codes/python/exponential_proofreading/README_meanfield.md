# Running the mean-field code

Two scripts in `Codes/python/exponential_proofreading/`, run with the repo's virtual environment:

```bash
cd ~/Documents/GitHub/Immune_System/Codes/python/exponential_proofreading
PY=../../../.venv/bin/python
$PY run_meanfield.py [options] [--plot]      # simulate (and optionally plot)
$PY plot_meanfield.py <run name>             # re-plot a saved run without re-running
```

- **Data** go to `My Drive/research_files/Immune_system/exponential_proofreading/meanfield/runs/<run name>/`.
- **Figures** go to `My Drive/documents/Science/research_projects/Immune_System/figures/exponential_proofreading/meanfield/runs/<run name>/`.
- **The run name is built from the options,** e.g. `two_stage_polyclonal_recall+naive_replay_sweep_h0_N_T0`. It is printed at the start of every run.

## 1. Choices on the command line

| Option | Values | Default |
|---|---|---|
| `--model` | `continuous`: antigen loading, T-cell help and division at the same time (the semicomplete equations)<br>`two_stage`: loading for t < t_switch, then expansion limited by a finite T-cell pool | `continuous` |
| `--repertoire` | `polyclonal`: affinity grid weighted by the density of states<br>`two_clone`: two explicit clones | `polyclonal` |
| `--responses` | `naive`, or `naive recall` | `naive recall` |
| `--antigen` | recall antigen law: `feedback` (the response controls the viral load) or `replay` (the primary's antigen trajectory, no feedback) | `feedback` (becomes `replay` if you pass `--with-naive`) |
| `--with-naive` / `--no-naive` | include the naive repertoire in the recall, or recall from memory only | `--no-naive` |
| `--draw` | memory seed: `expected` (deterministic) or `random` (multinomial) | `expected` |
| `--sweep NAME[=VALUES]` | up to twice. None means a single run; one is a 1-D sweep; two is a 2-D grid | none (single run) |
| `--set PARAM=VALUE` | override any model parameter or two-clone setting (repeatable) | none |
| `--T` | end time | 12 (16 for `two_stage`) |
| `--t-obs` | times at which the end-point observables are computed | `T` |
| `--name` | your own run name | automatic |
| `--plot` | make the figures after the run | off |
| `--dry-run` | print the resolved configuration and stop | off |

`python run_meanfield.py --help` lists all the options.

**Sweep values:**
- `--sweep NAME` alone uses a default grid; or give `NAME=1,2,3`, `NAME=log:a:b:n` (10^a to 10^b, n points), or `NAME=lin:a:b:n`.
- `pi_c` is converted to h0 = b0/π_c², so the run is saved as an h0 sweep.
- Any model parameter can be swept. Default grids:

| Parameter | Default grid |
|---|---|
| `pi_c` / `h0` | 16 values, π_c = 3.2 to 10⁴ |
| `n_mem` | 10² to 10⁶, 9 values |
| `N_T0` | 10³ to 10⁷, 9 values |
| `t_switch` | 4 to 8, 9 values |
| `h_max` | 10² to 10⁴, 5 values |
| `gap` | 0, 0.5, 1, 2, 4 |
| `ratio` | 0.1 to 10, 5 values |

**Presets** go as the first argument, and options after them override them:

| Preset | What it sets |
|---|---|
| `single` | single run with per-clone time courses (Fig. 2) |
| `boost` | recall with naive, antigen replay, sweep π_c × n_mem, observed at t = 6 and 12 |
| `two_stage` | two-stage model, recall with naive and replay, sweep π_c × N_T0, T = 16 |
| `two_clone` | two-clone repertoire, sweep over the gap |

**Run names:** parameters changed with `--set` become part of the name (e.g. `…_single_N_T0=1e+04`), so runs don't overwrite each other.

## 2. Default parameters (`BASE` in `run_meanfield.py`)

| | Parameter | Default |
|---|---|---|
| Antigen | growth λ_A, clearance δ_A | 6.2, 3.5 per day |
| | innate cutoff (rate, threshold) | 2.2, 10⁴ |
| | initial N_A0 | 1 |
| pMHC | turnover δ_π | 0.1 |
| | Hill coefficient | 2 |
| | k_on | 1e2·2e5·1e6·24·3600 / N_Avogadro |
| B cells | division rate b0 | 2 per day |
| | **h0** | **b0/1000, so π_c = (b0/h0)^½ ≈ 31.6** |
| | death δ_B | 0 |
| Recall | memory h0 | `h0_mem = b0` (`BASE` `h0_mem=None`), so π_c,mem = 1; change with `--set h0_mem=...` |
| | memory capacity n_mem | 4·10⁴ |
| | GC-entry bias f in φ_naive | 1 |
| Repertoire | ΔG grid | 0 to 9, 32 points |
| | β* | 2.3 |
| | η | 1 |
| T cells | N_T0 | 10⁶ |
| | K_T | 10⁴ (continuous model) |
| | τ_eng | 0.1 day |
| Two-stage | t_switch | 6 |
| | h_max | 10³ per day |
| Potency | Z_c | 10³ |
| Two-clone | clone 1 at ΔG1, clone 2 at ΔG1 + gap, sizes N1 and N1·ratio | ΔG1 = 0, gap = 1, N1 = 1, ratio = 1 |

## 3. Response modes

| Mode | Options |
|---|---|
| Primary only | `--responses naive` |
| Recall from memory, controlling the antigen (the legacy memory phase; used for the potency figure) | default |
| Recall from memory, antigen replayed | `--antigen replay` |
| Boost: naive + memory compete for the same antigen | `--with-naive` (replay) |
| Naive + memory with feedback | `--with-naive --antigen feedback` |

**Memory seed:** formed from the primary of the same parameters, about min(n_mem, primary yield) cells drawn in proportion to the expanded clones. If no clone expanded, the recall is skipped and logged.

## 4. Outputs

**`observables.csv`:** one row per (sweep point, response, t_obs). Columns: the swept parameters, `pi_c`, `response`, `t_obs`, then:

- **Per population** (`_naive`, `_memory`, `_all`):

| Column | Meaning |
|---|---|
| `Y` | cells in expanded clones (at least doubled) |
| `L_act` | number of expanded clones |
| `N1` | largest clone |
| `Z` | potency |
| `D_KL` | specificity |
| `S2` | Rényi-2 clone entropy, −ln Σ p² over expanded clones |
| `alpha_fit` | fitted EP exponent |
| `DG_front` | affinity front of the expanded clones |
| `DG_hi`, `DG_lo` | range of clones above 2 cells |
| `f_best` | fraction of cells in the best clone |

- **Per response:**

| Column | Meaning |
|---|---|
| `N_A` | antigen load |
| `T_free_min` | minimum free T-cell fraction |
| `t_c` | time the potency reaches Z_c |

- **Naive + memory only:**

| Column | Meaning |
|---|---|
| `R` | memory-derived ÷ naive-derived cells |
| `phi_naive` | 1/(1 + f·R) |
| `R_th` | analytic κ·ν_μ·N̄_w |
| `kappa` | per-cell expansion advantage of memory |

**`dynamics.pkl`:** time courses per run: N_A, T_free, and Z, Y and L_act per population. Single, two-clone and boost runs also keep per-clone N_B and π. With `DG_dist` (default and `single`), it also keeps the final cells per clone and the memory seed, for the energy distributions.

**`config.json`:** the full resolved configuration.

## 5. Figures (`--plot`, or `plot_meanfield.py <run>`)

- **`DG_dist/`:** cell density in ΔG (cells per unit ΔG, expanded clones) at the end of the primary, in the memory seed and at the end of the recall, one figure per sweep point plus a grid (single runs and 1-D sweeps). Theory shapes scaled to each curve's total cells: dotted, EP exponentials Ω₀e^{−aΔG} (a = α, or 2α for the recall); dashed grey, the Gaussian of Fig. 4D.
- **`zipf/`:** rank-size (Zipf) plots of the expanded clones (doubled since t = 0) at the end of each response, normalised to the largest clone, one curve per sweep point (1-D sweeps); dashed: k^−ζ (primary) and k^−2ζ (recall), ζ = η(b0/λ_A + b0/min(b0, δ_A))/β*. Single runs get `manuscript/zipf.pdf` with both responses.
- **`dynamics/`:** time courses, one curve per sweep point. For 2-D sweeps, one figure each at the first, middle and last value of the second parameter.
- **`observables/`:** each observable against the swept parameter (1-D sweep), or as a heatmap (2-D sweep).
- **`manuscript/`:** produced automatically when the run matches a figure:

| Figure | Command |
|---|---|
| Fig. 2 panels (N_A, π, N_B, Z, fronts) | `$PY run_meanfield.py --plot` |
| Potency vs specificity | `$PY run_meanfield.py --sweep pi_c --plot` |
| Boost (R vs π_c, φ_naive map, r(ΔG), R vs n_mem) | `$PY run_meanfield.py boost --plot` |
| Two-clone | `$PY run_meanfield.py two_clone --plot` |
| T-cell competition (two-stage) | `$PY run_meanfield.py --model two_stage --with-naive --sweep pi_c --sweep N_T0 --plot` |

Fig. 2 panels: the three clones shown are at ΔG = 0, 2.9 and 6.39 (grid rows 0, 10, −10).

## 6. From Python

```python
import run_meanfield as rm
df, dynamics = rm.run(dict(rm.CONFIG, model='two_stage', T=16.0, t_obs=[16.0], sweep={'N_T0': [1e4, 1e6]}))
rm.main(['--model', 'two_stage', '--plot'])        # same as the command line
```

## 7. Things to keep in mind

- **Observation times snap to the time grid.** Each run has 500 time points, so `--t-obs 6` reads t = 5.99.
- **Fig. 2 for the two-stage model:**
  - the fronts panel (D) is skipped, because its formula only applies to the continuous model;
  - the axis ranges are the continuous model's, and they clip some two-stage curves.
- **h_max matters when T cells are scarce:** two-stage results with low N_T0 depend on it (up to about 40%).
- **Plotting needs LaTeX installed,** because the figures use LaTeX text rendering.
- **The previous driver is kept:** `run_meanfield_legacy.py` still works, for comparison only.

## 8. Examples

All commands assume `cd ~/Documents/GitHub/Immune_System/Codes/python/exponential_proofreading` and `PY=../../../.venv/bin/python`. Add `--dry-run` to any of them to see the resolved configuration and run name without running.

### Before you start: check a configuration

```bash
$PY run_meanfield.py --model two_stage --with-naive --sweep pi_c --sweep N_T0 --dry-run
```
This prints the run name, the model and responses, the sweep (number of values and range), T, t_obs, the time courses kept, and any changed parameters. Nothing is run.

### Single runs

```bash
# Primary response only, default parameters (pi_c ~ 31.6)
$PY run_meanfield.py --responses naive --plot

# Primary + recall (Fig. 2 panels in manuscript/naive and manuscript/recall)
$PY run_meanfield.py --plot

# Same at a different T-cell help threshold: h0 = 2e-4 gives pi_c = (b0/h0)^(1/2) = 100
$PY run_meanfield.py --set h0=2e-4 --plot

# Two-stage model with scarce T cells (T = 16 is set automatically)
$PY run_meanfield.py --model two_stage --set N_T0=1e4 --plot
```

### One-parameter sweeps

```bash
# Primary responses across the T-cell help threshold: yield, reach, specificity vs pi_c
$PY run_meanfield.py --responses naive --sweep pi_c --plot

# Potency vs specificity (manuscript/potency.pdf): recall with antigen feedback across pi_c
$PY run_meanfield.py --sweep pi_c --plot

# A custom pi_c range: 9 values from 10 to 1000 (a new --name keeps the default run intact)
$PY run_meanfield.py --sweep pi_c=log:1:3:9 --name potency_pic_10_1000 --plot

# Recall without antigen feedback, to compare with the run above
$PY run_meanfield.py --antigen replay --sweep pi_c --plot

# Two-stage primary across the T-cell pool: where does T-cell limitation start?
$PY run_meanfield.py --model two_stage --responses naive --sweep N_T0 --plot

# Two-stage primary across the loading/expansion switch time
$PY run_meanfield.py --model two_stage --responses naive --sweep t_switch --plot

# Sensitivity to the engagement cap when T cells are scarce
$PY run_meanfield.py --model two_stage --with-naive --set N_T0=1e4 --sweep h_max --plot
```

### Boost: naive and memory together

```bash
# The full boost map, pi_c x memory capacity, read at day 6 and day 12
# (manuscript/t_read_6/ and t_read_12/: R vs pi_c, phi_naive map, r(DG), R vs n_mem)
$PY run_meanfield.py --with-naive --sweep pi_c --sweep n_mem --t-obs 6 12 --plot
#   (equivalent to the preset:  $PY run_meanfield.py boost --plot)

# A quicker version at a fixed memory capacity, read at day 6
$PY run_meanfield.py --with-naive --set n_mem=1e3 --sweep pi_c --t-obs 6 --plot

# Naive/memory competition for T cells: pi_c x T-cell pool, two-stage model
$PY run_meanfield.py --model two_stage --with-naive --sweep pi_c --sweep N_T0 --plot

# The boost map for the two-stage model
$PY run_meanfield.py boost --model two_stage --plot
```

### Two clones

```bash
# Affinity gap: how much does the better clone win? (manuscript/final_frac_clone1.pdf)
$PY run_meanfield.py --repertoire two_clone --sweep gap --plot

# Initial size ratio at a fixed gap: can a larger, weaker clone win?
$PY run_meanfield.py --repertoire two_clone --set gap=2 --sweep ratio=0.01,0.1,1,10 --plot

# Gap x size ratio as heatmaps
$PY run_meanfield.py --repertoire two_clone --sweep gap --sweep ratio --plot

# Two clones competing for T cells in the two-stage boost
$PY run_meanfield.py --model two_stage --repertoire two_clone --with-naive --sweep gap --plot
```

### Re-plotting and analysing a saved run

```bash
# Re-make the figures of a saved run without re-running it (e.g. after restyling a figure)
$PY plot_meanfield.py "continuous_polyclonal_recall+naive_replay_sweep_h0_n_mem"
```

```python
# Load the end-point observables of a run into pandas
import pandas as pd
from run_meanfield import DATA_ROOT
df = pd.read_csv(f'{DATA_ROOT}/two_stage_polyclonal_recall+naive_replay_sweep_h0_N_T0/observables.csv')
recall = df[df.response == 'recall']
print(recall.pivot_table(index='N_T0', columns='pi_c', values='R'))     # R over the 2-D grid
```

### Tips

- **Run names come from the options,** so two runs that differ only in the sweep values (e.g. `--sweep pi_c` and `--sweep pi_c=log:1:3:9`) get the same name and the second overwrites the first. Use `--name` to keep both. Parameters changed with `--set` are already part of the name.
- **2-D sweeps take longer:** the boost map (16 × 9 points) takes about 2 minutes, the two-stage π_c × N_T0 map about 1 minute. Single runs and 1-D sweeps take seconds.
- **Quotes:** names containing `+` (e.g. `recall+naive`) must be quoted when passed to `plot_meanfield.py`.
