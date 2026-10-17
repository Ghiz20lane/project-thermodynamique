# -*- coding: utf-8 -*-
# CODE 3 : Maxwell avec G + synthese (Etapes 6 et 7)
# Code independant : il redefinit ses propres donnees.

import numpy as np
import sympy as sp
import matplotlib.pyplot as plt

# --- Donnees du CO2 --------------------------------------------------
R = 8.314462618    # J/(mol.K)
a = 0.364          # Pa.m6/mol2
b = 4.267e-5       # m3/mol
T = 320.0          # K
v = 1.0e-3         # m3/mol


def p_ideal(T, v):
    return R * T / v


def p_vdw(T, v):
    if np.any(np.asarray(v) <= b):
        raise ValueError("Le volume doit etre superieur a b")
    return R * T / (v - b) - a / v**2


# ===== ETAPE 6 : 4e RELATION DE MAXWELL (avec G) =====================
# (dS/dP)_T = -(dv/dT)_P
# Derivation implicite : (dv/dT)_P = -(dP/dT)_v / (dP/dv)_T
Ts, vs, Rs, aa, bb = sp.symbols("T v R a b", positive=True)
P = Rs * Ts / (vs - bb) - aa / vs**2
dP_dT = sp.diff(P, Ts)
dP_dv = sp.diff(P, vs)

dvdT_P = sp.simplify(-dP_dT / dP_dv)
# subs : on remplace les lettres par les valeurs numeriques
valeur = dvdT_P.subs({Ts: T, vs: v, Rs: R, aa: a, bb: b})
print("VdW : (dv/dT)P =", float(valeur))
print("VdW : (dS/dP)T =", -float(valeur))

# Gaz parfait (a la main) : (dv/dT)_P = v/T
dvdT_ideal = v / T
print("Ideal : (dv/dT)P =", dvdT_ideal)
print("Ideal : (dS/dP)T =", -dvdT_ideal)

# ===== ETAPE 7 : SYNTHESE ============================================
v1, v2 = 1.0e-3, 2.0e-3
delta_U_vdw = a * (1 / v1 - 1 / v2)
dpdt_ideal = R / v
dpdt_vdw = R / (v - b)
dpdv_vdw = -R * T / (v - b)**2 + 2 * a / v**3
cpcv_vdw = -T * dpdt_vdw**2 / dpdv_vdw

print("\nTABLEAU DE SYNTHESE (T = 320 K, v = 1e-3 m3/mol)")
print(f"{'Grandeur':<22}{'Gaz parfait':>16}{'Van der Waals':>16}")
print(f"{'P (bar)':<22}{p_ideal(T, v)/1e5:>16.4f}{p_vdw(T, v)/1e5:>16.4f}")
print(f"{'(dP/dT)v (Pa/K)':<22}{dpdt_ideal:>16.2f}{dpdt_vdw:>16.2f}")
print(f"{'(dU/dv)T (Pa)':<22}{0.0:>16.2f}{a / v**2:>16.2f}")
print(f"{'Delta U (J/mol)':<22}{0.0:>16.2f}{delta_U_vdw:>16.2f}")
print(f"{'Cp-Cv (J/mol/K)':<22}{R:>16.4f}{cpcv_vdw:>16.4f}")
print(f"{'(dv/dT)P (m3/mol/K)':<22}{dvdT_ideal:>16.4e}{float(valeur):>16.4e}")

# Graphique final
volumes = np.linspace(5e-4, 5e-3, 250)
plt.plot(volumes, p_ideal(T, volumes) / 1e5, label="Gaz parfait")
plt.plot(volumes, p_vdw(T, volumes) / 1e5, label="Van der Waals")
plt.xlabel("Volume molaire (m3/mol)")
plt.ylabel("Pression (bar)")
plt.title("Isotherme a T = 320 K")
plt.grid(True)
plt.legend()
plt.tight_layout()
plt.show()