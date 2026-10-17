# Projet Thermodynamique : Relations de Maxwell (Python)

**Cycle Ingénieur Procédés et Ingénierie Chimique (PIC), FST Settat, Université Hassan 1er**
**Module :** Thermodynamique 1
**Enseignant :** Pr. BAKHER Zine Elabidine
**Étudiant(e) :** *Ghizlane derouich*

## Objectif

Établir les relations de Maxwell à partir des potentiels thermodynamiques (U, H, F, G), puis les appliquer à une mole de CO₂ décrite par deux modèles :

- le **gaz parfait** : P = RT/v
- le **gaz de Van der Waals** : P = RT/(v − b) − a/v²

Les calculs sont faits à la main (dans le rapport), puis vérifiés avec Python (NumPy, SymPy, Matplotlib).

## Structure du projet

```
projet-maxwell/
├── code1_equations_etat_derivees.py
├── code2_energie_interne_mayer.py
├── code3_maxwell_G_synthese.py
├── requirements.txt
└── README.md
```

| Fichier | Étapes du sujet | Contenu |
|---|---|---|
| `code1_equations_etat_derivees.py` | 2 et 3 | Pression des deux modèles, dérivée (∂P/∂T)ᵥ par SymPy et par différences finies, effet du pas h, isotherme |
| `code2_energie_interne_mayer.py` | 4 et 5 | (∂U/∂v)ₜ, ΔU pour une détente isotherme, Cₚ − Cᵥ, 2 graphiques |
| `code3_maxwell_G_synthese.py` | 6 et 7 | 4ᵉ relation de Maxwell (avec G), tableau de synthèse, graphique final |

Les trois scripts sont **indépendants** : chacun redéfinit ses constantes et ses fonctions, et peut être exécuté seul.

## Prérequis

- Python 3.8 ou plus
- Bibliothèques autorisées : NumPy, SymPy, Matplotlib

Installation :

```bash
pip install -r requirements.txt
```

## Exécution dans VS Code

1. Ouvrir le dossier du projet : **Fichier > Ouvrir le dossier**.
2. Installer l'extension **Python** (Microsoft) si ce n'est pas déjà fait.
3. Choisir l'interpréteur : `Ctrl + Maj + P`, puis **Python: Select Interpreter**.
4. Ouvrir un fichier `.py` et cliquer sur le bouton **▶ Run Python File** (en haut à droite), ou utiliser le terminal :

```bash
python code1_equations_etat_derivees.py
python code2_energie_interne_mayer.py
python code3_maxwell_G_synthese.py
```

**Remarque :** les graphiques s'ouvrent dans une fenêtre. Il faut la fermer pour que le programme continue.

## Données utilisées (CO₂, une mole)

| Symbole | Valeur | Unité |
|---|---|---|
| R | 8,314462618 | J/(mol·K) |
| a | 0,364 | Pa·m⁶/mol² |
| b | 4,267 × 10⁻⁵ | m³/mol |
| T | 320 | K |
| v (point de calcul) | 1,0 × 10⁻³ | m³/mol |

Plage de volumes des graphiques : de 5 × 10⁻⁴ à 5 × 10⁻³ m³/mol. On vérifie toujours v > b.

## Formules principales

- Maxwell (issue de F) : (∂S/∂v)ₜ = (∂P/∂T)ᵥ
- Maxwell (issue de G) : (∂S/∂P)ₜ = −(∂v/∂T)ₚ
- Énergie interne : (∂U/∂v)ₜ = T(∂P/∂T)ᵥ − P, soit 0 (gaz parfait) et a/v² (Van der Waals)
- ΔU (Van der Waals, isotherme) : ΔU = a(1/v₁ − 1/v₂)
- Mayer généralisée : Cₚ − Cᵥ = −T[(∂P/∂T)ᵥ]² / (∂P/∂v)ₜ
- Dérivation implicite : (∂v/∂T)ₚ = −(∂P/∂T)ᵥ / (∂P/∂v)ₜ

## Résultats attendus (T = 320 K, v = 1,0 × 10⁻³ m³/mol)

| Grandeur | Gaz parfait | Van der Waals |
|---|---|---|
| P (bar) | 26,606 | 24,152 |
| (∂P/∂T)ᵥ (Pa/K) | 8314,46 | 8685,05 |
| (∂U/∂v)ₜ (Pa) | 0 | 364 000 |
| ΔU, détente v₁ = 1×10⁻³ → v₂ = 2×10⁻³ (J/mol) | 0 | 182 |
| Cₚ − Cᵥ (J/(mol·K)) | 8,314 | 11,097 |
| (∂v/∂T)ₚ (m³/(mol·K)) | 3,125 × 10⁻⁶ | 3,993 × 10⁻⁶ |

## Remarques

- Le code 1 reprend le code fourni dans l'énoncé (étapes 2 et 3), avec des commentaires ajoutés.
- Pour Cₚ − Cᵥ, la formule n'a de sens que si (∂P/∂v)ₜ < 0 (phase stable). Près d'un point où cette dérivée s'annule, le résultat diverge et n'est pas interprété.
- Pour T = 320 K, supérieure à la température critique du modèle (environ 304 K), toute la plage étudiée est stable.
- Les démonstrations à la main et l'interprétation des figures sont dans le rapport remis avec ce projet.
