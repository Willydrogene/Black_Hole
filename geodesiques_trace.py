import matplotlib.pyplot as plt     # On importe matplotlib pour faire les graphiques
import numpy as np                  # On importe numpy pour gérer les tableaux et les NaN, etc.



# ------------------------------ Constantes physiques ------------------------------ #
G = 6.67430e-11                                 # Constante gravitationnelle universelle (SI)
c = 2.99792458e8                                # Vitesse de la lumière dans le vide (m/s)
M_sun = 1.98847e30                              # Masse du Soleil (kg)
masse_trou_noir = 5 * M_sun                     # Masse du trou noir : ici on prend 5 masses solaires
rs = 2 * G * masse_trou_noir/(c*c)              # Rayon de Schwarzschild : rs = 2GM / c² (ici juste pour référence, pas utilisé directement dans le plot)


# ------------------------------ Chargement des données ------------------------------ #
data = np.loadtxt("eq_diff_1.txt", skiprows=1)  # skiprows=1 pour ignorer l'en-tête texte de la première ligne


# ------------------------------ Configuration de la figure ------------------------------ #
plt.style.use("dark_background")                # On utilise un style "fond noir", plus adapté visuellement aux trous noirs
fig, ax = plt.subplots(figsize=(10,10))         # On crée une figure et un système d'axes, taille 10x10 (carrée)


# ------------------------------ Identification des séparateurs NaN ------------------------------ #
nan_indices = np.where(np.isnan(data[:,0]))[0]                  # On récupère les indices des lignes où x = NaN
indices = np.concatenate(([0], nan_indices, [len(data)]))       # On construit un tableau d'indices de "coupure" : début (0), NaN, et fin (len(data))


# ------------------------------ Boucle sur chaque segment de géodésique ------------------------------ #
for i in range(len(indices) - 1):
    start, end = indices[i], indices[i+1]       # On récupère le segment entre deux indices successifs (entre deux NaN ou entre début/fin)
    segment = data[start:end]
    
    
    segment = segment[~np.isnan(segment[:,0])]  # On supprime les lignes NaN du segment
    
    if len(segment) == 0:                       # Si, une fois les NaN enlevés, le segment est vide, on passe au suivant
        continue


    # --------------------------------- Gestion de la couleur --------------------------------- #
    val = int(segment[0,2])
    if val == 0:
        color = "y"         # uc > 1 (rayons "classiques") => on affiche en jaune
    elif val == 1:
        color = "r"         # uc = 1 (rayon critique) => en rouge
    elif val == 2:
        color = "b"         # uc < 1 (rayons qui plongent / trajectoires différentes) => en bleu
    else:
        color = "w"         # Sécurité : si jamais on a une autre valeur, on trace en blanc


    # --------------------------------- Tracé du segment --------------------------------- #
    ax.plot(segment[:,0],  segment[:,1],  color=color, alpha=0.6)       # Courbe originale
    ax.plot(segment[:,0], -segment[:,1], color=color, alpha=0.6)        # Partie symétrique par rapport à l’axe x


# ------------------------------ Mise en forme du graphique ------------------------------ #
ax.set_title("Position")                    # Titre de la figure
ax.set_xlabel("x")                          # Nom de l’axe des abscisses
ax.set_ylabel("y")                          # Nom de l’axe des ordonnées
ax.grid(True, color="grey", alpha=0.3)      # Grille légèrement visible (gris, un peu transparente)
ax.set_xlim(-150e3, 150e3)                  # Limites en x
ax.set_ylim(-150e3, 150e3)                  # Limites en y
ax.set_aspect('equal', adjustable='box')    # Rapport d’échelle identique en x et y (cercles restent des cercles)
plt.show()                                  # Affichage de la figure à l’écran

