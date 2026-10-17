# -*- coding: utf-8 -*-
# CODE 2 : Energie interne et relation de Mayer (Etapes 4 et 5)
# Code independant : il redefinit ses propres donnees.

import numpy as np
import matplotlib.pyplot as plt

# --- Donnees du CO2 --------------------------------------------------
R = 8.314462618    # J/(mol.K)
a = 0.364          # Pa.m6/mol2
b = 4.267e-5       # m3/mol
T = 320.0          # K
v = 1.0e-3         # m3/mol


def p_vdw(T, v):
    # Pression de Van der Waals en Pa (verifie v > b)
    if np.any(np.asarray(v) <= b):
        raise ValueError("Le volume doit etre superieur a b")
    return R * T / (v - b) - a / v**2


# ===== ETAPE 4 : ENERGIE INTERNE =====================================
# (dU/dv)_T = T (dP/dT)_v - P
#   gaz parfait   : 0   (loi de Joule)
#   Van der Waals : a/v^2
v1 = 1.0e-3   # volume initial, m3/mol
v2 = 2.0e-3   # volume final, m3/mol

# Integration de a/v^2 entre v1 et v2
delta_U_vdw = a * (1 / v1 - 1 / v2)
delta_U_ideal = 0.0
print("Delta U Van der Waals (J/mol) :", delta_U_vdw)
print("Delta U gaz parfait (J/mol)   :", delta_U_ideal)
# Signe positif : le gaz vainc les forces d'attraction en se detendant,
# son energie interne augmente.

# Verification avec la relation de Maxwell
dP_dT = R / (v - b)
print("T(dP/dT)v - P =", T * dP_dT - p_vdw(T, v), "Pa")
print("a/v^2         =", a / v**2, "Pa")

# Graphique de Delta U en fonction de v2
v2_tab = np.linspace(v1, 5e-3, 200)
dU_tab = a * (1 / v1 - 1 / v2_tab)

plt.figure()
plt.plot(v2_tab, dU_tab, label="Van der Waals")
plt.plot(v2_tab, np.zeros_like(v2_tab), "--", label="Gaz parfait (0)")
plt.xlabel("Volume final v2 (m3/mol)")
plt.ylabel("Delta U (J/mol)")
plt.title("Detente isotherme a T = 320 K, v1 = 1e-3 m3/mol")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()


# ===== ETAPE 5 : RELATION DE MAYER GENERALISEE =======================
# Cp - Cv = -T [(dP/dT)_v]^2 / (dP/dv)_T      (gaz parfait : R)
def cp_moins_cv_vdw(T, v):
    dpdt = R / (v - b)                              # (dP/dT)_v
    dpdv = -R * T / (v - b)**2 + 2 * a / v**3       # (dP/dv)_T
    # Si dP/dv ~ 0, on divise presque par zero : resultat sans sens
    if np.any(np.isclose(dpdv, 0, atol=1e-12)):
        raise ValueError("Derivee en volume proche de zero")
    return -T * dpdt**2 / dpdv


print("Cp-Cv ideal :", R)
print("Cp-Cv VdW   :", cp_moins_cv_vdw(T, v))

# Graphique de Cp - Cv en fonction de v (points stables seulement)
volumes = np.linspace(5e-4, 5e-3, 250)
dpdv_tab = -R * T / (volumes - b)**2 + 2 * a / volumes**3
masque = dpdv_tab < -1e-6                 # on garde les points stables
cpcv_tab = np.full_like(volumes, np.nan)  # nan = point non trace
cpcv_tab[masque] = cp_moins_cv_vdw(T, volumes[masque])

plt.figure()
plt.plot(volumes, cpcv_tab, label="Van der Waals")
plt.axhline(R, linestyle="--", color="gray", label="Gaz parfait (R)")
plt.xlabel("Volume molaire (m3/mol)")
plt.ylabel("Cp - Cv (J/(mol.K))")
plt.title("Cp - Cv a T = 320 K")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()