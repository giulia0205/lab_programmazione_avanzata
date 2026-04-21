import numpy as np
import matplotlib.pyplot as plt

mid = np.loadtxt("a_puntomedio.dat")
adp = np.loadtxt("a_adaptive_full.dat")

t_mid = mid[:, 0]
a_mid = mid[:, 1]

t_adp = adp[:, 0]
a_adp = adp[:, 1]

# Interpolazione adaptive sui tempi del midpoint
a_adp_on_mid = np.interp(t_mid, t_adp, a_adp)

err_abs = np.abs(a_mid - a_adp_on_mid)
err_rel = err_abs / np.maximum(np.abs(a_adp_on_mid), 1e-30)

one_plus_z_mid = 1.0 / a_mid
one_plus_z_adp = 1.0 / a_adp

print("max err_abs =", np.max(err_abs))
print("max err_rel =", np.max(err_rel))

# =========================
# 1) a(t)
# =========================
plt.figure(figsize=(8, 6))
plt.plot(t_mid, a_mid, label="Punto medio")
plt.plot(t_adp, a_adp, "--", label="Passo adattivo")
plt.xlabel("t [Gyr]")
plt.ylabel("a(t)")
plt.title("Evoluzione del fattore di scala")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()

# =========================
# 2) 1+z su scala log
# =========================
plt.figure(figsize=(8, 6))
plt.plot(t_mid, one_plus_z_mid, label="Punto medio")
plt.plot(t_adp, one_plus_z_adp, "--", label="Passo adattivo")
plt.yscale("log")
plt.xlabel("t [Gyr]")
plt.ylabel("1 + z = 1/a")
plt.title("Evoluzione del redshift")
plt.grid(True, which="both")
plt.legend()
plt.tight_layout()
plt.show()

# =========================
# 3) errore assoluto
# =========================
mask_abs = err_abs > 0

plt.figure(figsize=(8, 6))
plt.plot(t_mid[mask_abs], err_abs[mask_abs])
plt.yscale("log")
plt.xlabel("t [Gyr]")
plt.ylabel(r"$|a_{\rm mid} - a_{\rm adapt}|$")
plt.title("Differenza assoluta tra i due metodi")
plt.grid(True, which="both")
plt.tight_layout()
plt.show()

# =========================
# 4) errore relativo
# =========================
mask_rel = err_rel > 0

plt.figure(figsize=(8, 6))
plt.plot(t_mid[mask_rel], err_rel[mask_rel])
plt.yscale("log")
plt.xlabel("t [Gyr]")
plt.ylabel(r"$|a_{\rm mid} - a_{\rm adapt}| / a_{\rm adapt}$")
plt.title("Differenza relativa tra i due metodi")
plt.grid(True, which="both")
plt.tight_layout()
plt.show()