import matplotlib.pyplot as plt
import numpy as np

# Time axis
T = 10
t = np.linspace(0, T, 500)

# Parameters
lambda_A = 0.5
lambda_B = 1.0
t_star = 4  # "deterministic threshold" time
T_act = 5   # stochastic activation time in FPM (example)

# Continuous-growth model: growth rate ramps up smoothly (sigmoid for illustration)
R_p = 1 / (1 + np.exp(-2*(t - t_star)))  # smooth rise around t_star
lnC_cgm = lambda_B * np.cumsum(R_p) * (t[1]-t[0])

# Step-function approximation
lnC_step = lambda_B * np.maximum(0, t - t_star)

# First-passage model (single realization): no growth until T_act, then linear
lnC_fpm = lambda_B * np.maximum(0, t - T_act)

# Plot
plt.figure(figsize=(8,5))
plt.plot(t, lnC_cgm, label="Continuous-growth (CGM)", lw=2)
plt.plot(t, lnC_step, "--", label="Step-function approx.", lw=2)
plt.plot(t, lnC_fpm, ":", label="First-passage (example path)", lw=2)

plt.axvline(t_star, color="gray", ls="--", lw=1)
plt.axvline(T_act, color="gray", ls=":", lw=1)
plt.text(t_star, 0.5, r"$t^\star$", ha="center", va="bottom")
plt.text(T_act, 0.5, r"$T_{\mathrm{act}}$", ha="center", va="bottom")

plt.xlabel("time $t$")
plt.ylabel(r"$\ln C(t)$ (clone size log)")
plt.title("Comparison of growth trajectories for one lineage")
plt.legend()
plt.grid(alpha=0.3)
plt.show()
