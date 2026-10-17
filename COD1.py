# Importation de la bibliothèque NumPy pour les calculs numériques
import numpy as np

# Constante universelle des gaz parfaits en J/(mol.K)
R = 8.314462618

# Constantes de Van der Waals
a = 0.364
b = 4.267e-5

# Fonction qui calcule la pression d'un gaz parfait
# T : température en kelvins
# V : volume molaire en m³/mol
def p_ideal(T, V):
    return R * T / V

# Fonction qui calcule la pression selon le modèle de Van der Waals
def p_vdw(T, V):

    # Vérification que le volume est supérieur à b
    if np.any(np.asarray(V) <= b):
        # Affichage d'une erreur si le volume est trop petit
        raise ValueError("Le volume doit être supérieur à b")

    # Calcul de la pression avec l'équation de Van der Waals
    return R * T / (V - b) - a / V**2

# Définition de la température en kelvins
T = 320.0

# Définition du volume molaire en m³/mol
V = 1.0e-3

# Calcul et affichage de la pression du gaz parfait en bar
print("Gaz parfait (bar) :", p_ideal(T, V) / 1e5)

# Calcul et affichage de la pression de Van der Waals en bar
print("Van der Waals (bar) :", p_vdw(T, V) / 1e5)