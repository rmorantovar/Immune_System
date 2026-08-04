# Synthetic CCDFs for Model B (rate advantage, constant signal)
# We'll show two figures:
# 1) Stochastic Yule mixture CCDF over clone size N
# 2) Deterministic growth CCDF over clone size x, highlighting the spike at C_max

import numpy as np
import matplotlib.pyplot as plt

# ----- Parameters -----
np.random.seed(0)
K0 = 1.0                 # k_on * S0
lambda_B = 0.8           # post-activation proliferation rate
T = 10.0                 # observation time
beta = 2.0               # exponent in log-density Omega0(k) ~ k^beta (high-affinity tail)
kmin, kmax = 1e-3, 1e2   # affinity range (k_off)
n_samples = 200_000      # number of affinity samples for the mixture

# ----- Sample k_off from density per log-scale proportional to k^beta -----
# On u=log k, density ∝ e^{beta * u} over [log kmin, log kmax] → sample by inverse transform.
umin, umax = np.log(kmin), np.log(kmax)
# CDF for U: F(u) = (e^{beta u} - e^{beta umin}) / (e^{beta umax} - e^{beta umin}) for beta != 0
u = None
if beta != 0:
    r = np.random.rand(n_samples)
    num = r*(np.exp(beta*umax) - np.exp(beta*umin)) + np.exp(beta*umin)
    u = np.log(num)/beta
else:
    # beta == 0 → uniform in log-space
    u = np.random.rand(n_samples)*(umax-umin) + umin

k_off = np.exp(u)

# ----- Model B growth rate and sizes -----
def growth_rate(k):
    return lambda_B * K0 / (K0 + k)

gk = growth_rate(k_off)
M = gk * T

# ----- 1) Stochastic Yule mixture CCDF -----
# For integer clone size N, CCDF(N) = E_k [ (1 - e^{-M(k)})^{N-1} ]
N_max = 3000
N_vals = np.arange(1, N_max+1, dtype=int)

p_k = 1 - np.exp(-M)  # geometric "success" parameter per k
# To avoid numerical underflow for large N, compute in log domain then exponentiate
log_terms = np.log(p_k + 1e-300)  # small epsilon to avoid log(0) for k with tiny M
# CCDF[N-1] = mean(exp((N-1) * log_terms))
# We'll compute for a subset of N to keep it fast and smooth
N_plot = np.unique(np.concatenate([np.arange(1,101), np.arange(100,1001,10), np.arange(1000, N_max+1, 50)]))
ccdf_stoch = np.array([np.mean(np.exp((n-1)*log_terms)) for n in N_plot])

plt.figure(figsize=(7,5))
plt.loglog(N_plot, ccdf_stoch)
plt.xlabel("Clone size N")
plt.ylabel("P(C ≥ N)  (stochastic Yule mixture)")
plt.title("Model B: CCDF over clone size (stochastic)")
plt.tight_layout()
plt.show()

# ----- 2) Deterministic growth CCDF -----
# Deterministic size for each k: C_det(k) = exp(M(k))
C_det = np.exp(M)
# Compute empirical CCDF over x by sorting
C_sorted = np.sort(C_det)
empirical_ccdf = 1.0 - np.arange(1, n_samples+1)/n_samples

plt.figure(figsize=(7,5))
plt.loglog(C_sorted, empirical_ccdf)
# Mark the saturation ceiling C_max = e^{lambda_B T}
C_max = np.exp(lambda_B * T)
plt.axvline(C_max, linestyle="--")
plt.xlabel("Clone size x")
plt.ylabel("P(C ≥ x)  (deterministic mixture)")
plt.title("Model B: CCDF over clone size (deterministic)\nNote spike at C_max (vertical dashed line)")
plt.tight_layout()
plt.show()
