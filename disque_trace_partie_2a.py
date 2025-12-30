import matplotlib.pyplot as plt     # On importe matplotlib
import numpy as np                  # On importe numpy

data = np.loadtxt("image_detecteur.txt", skiprows=1)  # On lit le fichier texte en ignorant la première ligne (entête)

plt.style.use("dark_background")    # Arrière-plan noir
plt.figure(figsize=(10,10))         # On ouvre une fenêtre de 10x10

b = data[:,0]           # On initialise b comme étant la première colonne de notre fichier
alpha = data[:,1]       # On initialise alpha comme étant la deuxième colonne de notre fichier
r = data[:,2]           # couleur !  # On intialise r comme étant la troisième colonne de notre fichier

x = b * np.cos(alpha)   # On calcule x pour mettre en coordonnées cartésiennes
y = b * np.sin(alpha)   # On calcule y

plt.scatter(
    x, y,
    c=r,                # couleur = r
    cmap="plasma",      # On peut mettre : plasma, inferno, viridis...
    s=5,                # taille des points
    alpha=0.8,          # Transparence
)

plt.colorbar(label="r = b / sin(theta)")  # barre de couleur

plt.title("Projection du disque d'accrétion sur le détecteur", fontsize=16)     #Titre du plot
plt.xlabel("x (coordonnée détecteur)", fontsize=12)                             # Nom de l'axe des x
plt.ylabel("y (coordonnée détecteur)", fontsize=12)                             # Nom de l'axe des y

plt.xlim(-10, 10)               # On fixe les limites de l'axe des x
plt.ylim(-10, 10)               # On fixe les limites de l'axe des y
plt.grid(False)                 # On précise qu'on ne veut pas de grille sur le plot
plt.gca().set_aspect("equal")   # On précise qu'on veut un repère orthonormé
plt.show()                      # On affiche notre plot