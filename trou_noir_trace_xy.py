import matplotlib.pyplot as plt     # On importe matplotlib
import numpy as np                  # On importe numpy

data = np.loadtxt("image_detecteur.txt", skiprows=1) # Chargement des données issues du C, et on ignore la première ligne (en-tête).

plt.style.use("dark_background")    # Arrière-plan noir
plt.figure(figsize=(12,12))         # On ouvre une fenêtre de 10x10

x = data[:,0]                       # Paramètre d’impact (coordonnée radiale sur le détecteur)
y = data[:,1]                       # Angle polaire sur le détecteur
r = data[:,2]                       # Rayon du disque où le photon a été émis
Fobs = data[:,3]                    # Flux observé F_obs = F_émis / (1+z)^4
UnpZ = data[:,4]                    # Facteur (1+z)


plt.scatter(
    x, y,
    c=Fobs,                         # Couleur = flux observé (ce que verrait réellement l'observateur)
    cmap="inferno",                 # essai : plasma, inferno, viridis...
    s=4,                            # Taille des points (petits = plus lisse)
    alpha=0.7,                      # Opacité totale
    marker ='.'                     # Petits points légers
)

plt.colorbar(label= r"Flux Observé ($F_{\text{obs}}$)")     # Barre de couleur (échelle des intensités)

plt.title(r"Lumière observée par un disque d'accrétion")    # Titre du plot
plt.xlabel("x (coordonnée détecteur)", fontsize=12)         # Nom de l'axe des x
plt.ylabel("y (coordonnée détecteur)", fontsize=12)         # Nom de l'axe des y

plt.xlim(-25, 25)                                           # On fixe les limites de l'axe des x
plt.ylim(-25, 25)                                           # On fixe les limites de l'axe des y
plt.grid(False)                                             # On précise qu'on ne veut pas de grille sur le plot
plt.gca().set_aspect("equal")                               # On précise qu'on veut un repère orthonormé
plt.show()                                                  # On affiche notre plot

