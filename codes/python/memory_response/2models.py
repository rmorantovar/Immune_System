import numpy as np
import matplotlib.pyplot as plt

rng = np.random.default_rng(3)

# ---------------------- Parameters ----------------------
T = 48.0                 # hours (observation time)
lamA = 0.25              # 1/h (signal growth rate)
lamB = 0.06              # 1/h (post-start division rate)
U0_h = 0.002 * 3600      # per hour (kon*S0 converted)
kp = 0.01                # 1/h (per-step rate)
p = 3                    # number of steps

# Choose a single clone k for panel (1)
k = 0.1                  # 1/h (off-rate; moderate affinity)

def q_of_k(k):
    return (kp/(kp + k))**p

def factor_fpam(k):
    return (U0_h / lamA) * q_of_k(k)

def t_star(k):
    return (1/lamA) * np.log(k / U0_h)

# ----------------- 1) Single-clone CSDs -----------------
N = 100000

# FPAM sampling
fac = factor_fpam(k)
U = rng.random(N)
Tact = (1/lamA) * np.log(1 + (-np.log(1-U))/fac)

def sample_yule_sizes(time_left):
    pgeom = 1 - np.exp(-lamB * time_left)
    U2 = rng.random(len(pgeom))
    n = 1 + np.floor(np.log(1-U2) / np.log(1 - pgeom + 1e-300)).astype(int)
    return n

mask_act = Tact <= T
sizes_fpam = np.ones(N, dtype=int)
if mask_act.any():
    time_left = T - Tact[mask_act]
    sizes_fpam[mask_act] = sample_yule_sizes(time_left)

# SPAM sampling
tstar = t_star(k)
sizes_spam = np.ones(N, dtype=int)
if T > tstar:
    M = lamB * q_of_k(k) * (T - tstar)
    pgeom = 1 - np.exp(-M)
    U3 = rng.random(N)
    sizes_spam = 1 + np.floor(np.log(1-U3) / np.log(1 - pgeom + 1e-300)).astype(int)

def ccdf_from_samples(x):
    x_sorted = np.sort(x)
    ccdf = 1.0 - np.arange(1, len(x_sorted)+1)/len(x_sorted)
    return x_sorted, ccdf

xF, cF = ccdf_from_samples(sizes_fpam)
xS, cS = ccdf_from_samples(sizes_spam)

plt.figure(figsize=(7,5))
plt.loglog(xF, cF, label="FPAM (single clone)")
plt.loglog(xS, cS, label="SPAM (single clone)")
plt.xlabel("Clone size C")
plt.ylabel("P(C ≥ x)")
plt.title(f"Single-clone CCDF at T={T} h (k={k} h$^{{-1}}$, p={p})")
plt.legend()
plt.tight_layout()
plt.show()

# ----------------- 2) Pooled over a repertoire Ω0(k) -----------------
kmin, kmax = 1e-3, 10.0     # 1/h
beta = 2.0
n_clones = 40000
umin, umax = np.log(kmin), np.log(kmax)
r = rng.random(n_clones)
u = np.log(r*(np.exp(beta*umax)-np.exp(beta*umin)) + np.exp(beta*umin))/beta
k_pool = np.exp(u)

fac_pool = factor_fpam(k_pool)
U = rng.random(n_clones)
Tact_pool = (1/lamA) * np.log(1 + (-np.log(1-U))/fac_pool)
sizes_fpam_pool = np.ones(n_clones, dtype=int)
mask_act_pool = Tact_pool <= T
if mask_act_pool.any():
    time_left = T - Tact_pool[mask_act_pool]
    sizes_fpam_pool[mask_act_pool] = sample_yule_sizes(time_left)

tstar_pool = t_star(k_pool)
sizes_spam_pool = np.ones(n_clones, dtype=int)
mask_on = T > tstar_pool
if mask_on.any():
    M_pool = lamB * q_of_k(k_pool[mask_on]) * (T - tstar_pool[mask_on])
    pgeom_pool = 1 - np.exp(-M_pool)
    U = rng.random(mask_on.sum())
    sizes_spam_pool[mask_on] = 1 + np.floor(np.log(1-U) / np.log(1 - pgeom_pool + 1e-300)).astype(int)

xF2, cF2 = ccdf_from_samples(sizes_fpam_pool)
xS2, cS2 = ccdf_from_samples(sizes_spam_pool)

plt.figure(figsize=(7,5))
plt.loglog(xF2, cF2, label="FPAM pooled Ω0(k)")
plt.loglog(xS2, cS2, label="SPAM pooled Ω0(k)")
plt.xlabel("Clone size C")
plt.ylabel("P(C ≥ x)")
plt.title(f"Pooled CCDFs (β={beta}, n={n_clones}, T={T} h)")
plt.legend()
plt.tight_layout()
plt.show()
