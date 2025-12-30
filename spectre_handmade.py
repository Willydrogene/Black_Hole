import matplotlib.pyplot as plt
import numpy as np
import tkinter as tk
from tkinter import simpledialog

data = np.loadtxt("image_detecteur_30_deg.txt", skiprows=1)


b = data[:,0]
alpha = data[:,1]
r = data[:,2] 
Fobs = data[:,3]
UnpZ = data[:,4]

nu_0 = 1.0
I_const = 1.0

#####################################détermination des grandeurs#####################################
I_0_r = I_const/(r**2)

#intégration sur bdalphadb

# Élément de surface sur le détecteur
# dS = b db dα
db = np.mean(np.diff(np.sort(np.unique(b))))                             # résolution en b
dalpha = np.mean(np.diff(np.sort(np.unique(alpha))))                     # résolution en alpha
dS = b * db * dalpha

#calcul grandeur Intensité et flux observé

I_obs = I_0_r/(UnpZ**3)                                                  #en un point du détecteur alpha b qui a UnpZ associé
Flux_obs = I_obs*dS                                                      #flux sur la totalité du détecteur

#calcul fréquence émise qui sera décalé par redshift = fréquence observé
nu_obs = nu_0/(UnpZ)

#####################################initialisation de l'histogramme des fréquences#####################################

nu_min = np.min(nu_obs)                         #bornes de l'intervalle des fréquences observé
nu_max = np.max(nu_obs)

range_bins = (nu_min, nu_max)     
bin_number = 200                               #nombre de sous intervalles
h = np.abs(nu_max - nu_min)/bin_number          #largueur des sous intervalles

nu = nu_min + h                                 #initialisation de la première valeurs de fréquences

F_tot_tab = np.array([])                        #initialisation des axes
nu_bin_center = np.array([])

for i in range(1,bin_number):
    
    F_tot = 0
    nu_center = nu-h/2                                       #on centre la classe (moyenne)
    nu_bin_center = np.append(nu_bin_center,nu_center)       #on stock la valeur

    mask = (nu_obs >= nu-h) & (nu_obs < nu)                 #on applique un masque sur les données pour extraire les valeurs ciblé et les indices
    selection = nu_obs[mask]
    indices = np.where(mask)[0]
    F_tot = np.sum(Flux_obs[mask])
    F_tot_tab = np.append(F_tot_tab,F_tot)
    nu += h                                                  #on passe à l'étape suivante

plt.plot(nu_bin_center, F_tot_tab)
plt.xlabel(r"Fréquence observée $\nu$")
plt.ylabel(r"Flux F($\nu$)")
plt.axvline(x=nu_0, color='red', linestyle='--', alpha=0.6, label=r'Fréquence émise $\nu_0$')   


plt.style.use("dark_background")    # Arrière-plan noir

x = b * np.cos(alpha)               # Coordonnées cartésiennes x
y = b * np.sin(alpha)               # Coordonnées cartésiennes y

z = UnpZ-1                          #Redshift z
limit1= abs(np.max(z))              #on détermine les bornes de la fluctuation du redshift pour que le blanc soit centré en 0       
limit2= abs(np.min(z))              #blanc = pas de distorsion des redshift
limit = 0   
if limit1>limit2:
    limit = limit1
else:
    limit = limit2

plt.figure(figsize=(10,10)) 

plt.scatter(
    x, y,
    c=UnpZ-1,                       # Couleur = flux observé (ce que verrait réellement l'observateur)
    cmap="seismic",                 # essai : plasma, inferno, viridis...    RdBu_r pour la carte des Redshifts,  hot, afmhot
    s=4,                            # Taille des points (petits = plus lisse)
    alpha=1,                        # Opacité totale
    marker ='.',                    # Petits points légers
    vmin=-limit,
    vmax=limit
)

plt.colorbar(label= r"")                                    # Barre de couleur (échelle des intensités)

plt.xlabel("x (coordonnée détecteur)", fontsize=12)         # Nom de l'axe des x
plt.ylabel("y (coordonnée détecteur)", fontsize=12)         # Nom de l'axe des y

plt.xlim(-25, 25)                                           # On fixe les limites de l'axe des x
plt.ylim(-25, 25)                                           # On fixe les limites de l'axe des y
plt.grid(False)                                             # On précise qu'on ne veut pas de grille sur le plot
plt.gca().set_aspect("equal")                               # On précise qu'on veut un repère orthonormé


plt.figure(figsize=(10,10)) 

plt.scatter(
    x, y,
    c=Fobs,                         # Couleur = flux observé (ce que verrait réellement l'observateur)
    cmap="afmhot",                  # essai : plasma, inferno, viridis...    RdBu_r pour la carte des Redshifts,  hot, afmhot
    s=4,                            # Taille des points (petits = plus lisse)
    alpha=0.7,                      # Opacité totale
    marker ='.',                    # Petits points légers
)

plt.colorbar(label= r"")                                    # Barre de couleur (échelle des intensités)

plt.xlabel("x (coordonnée détecteur)", fontsize=12)         # Nom de l'axe des x
plt.ylabel("y (coordonnée détecteur)", fontsize=12)         # Nom de l'axe des y

plt.xlim(-25, 25)                                           # On fixe les limites de l'axe des x
plt.ylim(-25, 25)                                           # On fixe les limites de l'axe des y
plt.grid(False)                                             # On précise qu'on ne veut pas de grille sur le plot
plt.gca().set_aspect("equal")                               # On précise qu'on veut un repère orthonormé
plt.show()
